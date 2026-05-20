// sinGRAAL — distributed DP coordinator
//
// Architecture:
//   coordinator machine:  kangaroo --serve [--bind 0.0.0.0:5135]
//   worker  machines:     kangaroo --coordinator <host:port> [--all-gpus]
//
// Protocol (binary, little-endian):
//   Handshake — worker→coord: [MAGIC "SGR2"]
//               coord→worker: [MAGIC "SGR2"]
//
//   Worker→Coordinator per batch:
//     [n_dps: u32]
//     [n_dps × DPWire (68 bytes each)]
//
//   Coordinator→Worker response:
//     [0x00]               = OK, keep going
//     [0x01][key: 32 bytes] = FOUND — here is the secret key
//
// The coordinator maintains the global DP table.  Workers just send DPs
// and act on FOUND signals.  Horizontal scaling is free: spin up N more
// GPU machines, point them at the same coordinator, and throughput scales
// linearly with zero coordination overhead.

use std::collections::HashMap;
use std::io::{BufReader, BufWriter, Read, Write};
use std::net::{TcpListener, TcpStream};
use std::sync::{
    atomic::{AtomicBool, AtomicU64, Ordering},
    Arc, Mutex,
};
use std::thread;
use std::time::Instant;

use crate::secp::{fe_to_hex, Pt};
use crate::{process_dp, DpRecord};

// ─── Wire constants ───────────────────────────────────────────────────────────

const MAGIC: [u8; 4] = *b"SGR2";

/// On-wire size of one DP entry (no padding bytes sent)
/// canon_x: 4×8 = 32, scalar: 4×8 = 32, is_wild: 4  → 68 bytes
pub const DP_WIRE: usize = 68;

const MSG_OK:    u8 = 0x00;
const MSG_FOUND: u8 = 0x01;

// ─── Shared state ─────────────────────────────────────────────────────────────

type Table  = Arc<Mutex<HashMap<[u64; 4], DpRecord>>>;
type Result = Arc<Mutex<Option<[u64; 4]>>>;

fn fe_from_wire(b: &[u8]) -> [u64; 4] {
    [
        u64::from_le_bytes(b[ 0.. 8].try_into().unwrap()),
        u64::from_le_bytes(b[ 8..16].try_into().unwrap()),
        u64::from_le_bytes(b[16..24].try_into().unwrap()),
        u64::from_le_bytes(b[24..32].try_into().unwrap()),
    ]
}

fn fe_to_wire(f: [u64; 4]) -> [u8; 32] {
    let mut b = [0u8; 32];
    for i in 0..4 { b[i*8..(i+1)*8].copy_from_slice(&f[i].to_le_bytes()); }
    b
}

// ─── Coordinator server ───────────────────────────────────────────────────────

/// Run a coordinator server.  Blocks until a solution is found or the process
/// is killed.  Workers connect, send DPs, and receive FOUND signals.
pub fn serve(bind: &str, target: Pt) -> Option<[u64; 4]> {
    let listener = TcpListener::bind(bind)
        .unwrap_or_else(|e| panic!("coordinator: cannot bind {bind}: {e}"));
    eprintln!(
        "[coord] listening on {bind} | target x=0x{}",
        fe_to_hex(target.x)
    );

    let table:    Table  = Arc::new(Mutex::new(HashMap::new()));
    let result:   Result = Arc::new(Mutex::new(None));
    let found:    Arc<AtomicBool> = Arc::new(AtomicBool::new(false));
    let dp_total: Arc<AtomicU64>  = Arc::new(AtomicU64::new(0));
    let workers:  Arc<AtomicU64>  = Arc::new(AtomicU64::new(0));
    let t0 = Instant::now();

    // Progress printer thread
    {
        let found  = found.clone();
        let dpt    = dp_total.clone();
        let wk     = workers.clone();
        let table  = table.clone();
        thread::spawn(move || {
            let mut last = 0u64;
            loop {
                thread::sleep(std::time::Duration::from_secs(10));
                if found.load(Ordering::Relaxed) { break; }
                let now   = dpt.load(Ordering::Relaxed);
                let rate  = (now - last) as f64 / 10.0;
                last = now;
                eprintln!(
                    "[coord] {:.1}M DPs total | {:.1}k DP/s | {} workers | table={}",
                    now as f64 / 1e6,
                    rate / 1e3,
                    wk.load(Ordering::Relaxed),
                    table.lock().unwrap().len(),
                );
            }
        });
    }

    for stream in listener.incoming() {
        if found.load(Ordering::Relaxed) { break; }
        let Ok(s) = stream else { continue };

        let table   = table.clone();
        let result  = result.clone();
        let found   = found.clone();
        let dp_cnt  = dp_total.clone();
        let wk_cnt  = workers.clone();

        thread::spawn(move || {
            wk_cnt.fetch_add(1, Ordering::Relaxed);
            worker_session(s, table, result, found, dp_cnt, target);
            wk_cnt.fetch_sub(1, Ordering::Relaxed);
        });
    }

    let _ = t0;
    let k = result.lock().unwrap().take();
    k
}

fn worker_session(
    stream:   TcpStream,
    table:    Table,
    result:   Result,
    found:    Arc<AtomicBool>,
    dp_total: Arc<AtomicU64>,
    target:   Pt,
) {
    let peer = stream.peer_addr()
        .map(|a| a.to_string())
        .unwrap_or_else(|_| "?".into());

    let mut r = BufReader::new(stream.try_clone().unwrap());
    let mut w = BufWriter::new(stream);

    // Handshake
    let mut m = [0u8; 4];
    if r.read_exact(&mut m).is_err() || m != MAGIC { return; }
    if w.write_all(&MAGIC).is_err() || w.flush().is_err() { return; }
    eprintln!("[coord] +{peer}");

    let mut dp_buf = Vec::new();

    loop {
        // Propagate FOUND to this worker if another thread solved it
        if found.load(Ordering::Acquire) {
            if let Some(k) = *result.lock().unwrap() {
                let mut msg = [0u8; 33];
                msg[0] = MSG_FOUND;
                msg[1..].copy_from_slice(&fe_to_wire(k));
                let _ = w.write_all(&msg);
                let _ = w.flush();
            }
            break;
        }

        // Read batch header
        let mut nb = [0u8; 4];
        if r.read_exact(&mut nb).is_err() { break; }
        let n = u32::from_le_bytes(nb) as usize;

        // Read DPs
        if n > 0 {
            dp_buf.resize(n * DP_WIRE, 0);
            if r.read_exact(&mut dp_buf).is_err() { break; }
            dp_total.fetch_add(n as u64, Ordering::Relaxed);
        }

        // Skip processing if already found by another thread
        if found.load(Ordering::Relaxed) { continue; }

        // Process against DP table
        let mut solved: Option<[u64; 4]> = None;
        if n > 0 {
            let mut tbl = table.lock().unwrap();
            for i in 0..n {
                let off = i * DP_WIRE;
                let cx  = fe_from_wire(&dp_buf[off..]);
                let sc  = fe_from_wire(&dp_buf[off + 32..]);
                let iw  = u32::from_le_bytes(
                    dp_buf[off + 64..off + 68].try_into().unwrap()
                ) != 0;
                if let Some(k) = process_dp(cx, sc, iw, &mut tbl, target) {
                    solved = Some(k);
                    break;
                }
            }
        }

        if let Some(k) = solved {
            found.store(true, Ordering::Release);
            *result.lock().unwrap() = Some(k);
            eprintln!("[coord] *** SOLVED *** key=0x{}", fe_to_hex(k));
            let mut msg = [0u8; 33];
            msg[0] = MSG_FOUND;
            msg[1..].copy_from_slice(&fe_to_wire(k));
            let _ = w.write_all(&msg);
            let _ = w.flush();
            break;
        }

        // ACK
        if w.write_all(&[MSG_OK]).is_err() || w.flush().is_err() { break; }
    }
    eprintln!("[coord] -{peer}");
}

// ─── Worker-side coordinator connection ───────────────────────────────────────

/// Persistent TCP connection from a GPU worker to the coordinator.
#[allow(dead_code)]
pub struct CoordConn {
    r: BufReader<TcpStream>,
    w: BufWriter<TcpStream>,
    addr: String,
}

#[allow(dead_code)]
impl CoordConn {
    /// Connect to coordinator and perform handshake.
    pub fn connect(addr: &str) -> std::io::Result<Self> {
        let s = TcpStream::connect(addr)?;
        s.set_nodelay(true)?;
        let mut w = BufWriter::new(s.try_clone()?);
        let mut r = BufReader::new(s);

        // Send magic, receive magic
        w.write_all(&MAGIC)?;
        w.flush()?;
        let mut m = [0u8; 4];
        r.read_exact(&mut m)?;
        if m != MAGIC {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "coordinator handshake failed",
            ));
        }
        eprintln!("[worker] connected to coordinator {addr}");
        Ok(CoordConn { r, w, addr: addr.to_owned() })
    }

    /// Send a batch of DPs to the coordinator.
    ///
    /// Each entry is (canon_x, scalar, is_wild).
    /// Returns `Some(key)` if the coordinator has found the solution.
    pub fn send_batch(
        &mut self,
        dps: &[([u64; 4], [u64; 4], bool)],
    ) -> std::io::Result<Option<[u64; 4]>> {
        let n = dps.len() as u32;
        self.w.write_all(&n.to_le_bytes())?;
        for &(cx, sc, iw) in dps {
            for &limb in &cx { self.w.write_all(&limb.to_le_bytes())?; }
            for &limb in &sc { self.w.write_all(&limb.to_le_bytes())?; }
            self.w.write_all(&(iw as u32).to_le_bytes())?;
        }
        self.w.flush()?;

        let mut resp = [0u8; 1];
        self.r.read_exact(&mut resp)?;
        if resp[0] == MSG_FOUND {
            let mut kb = [0u8; 32];
            self.r.read_exact(&mut kb)?;
            return Ok(Some(fe_from_wire(&kb)));
        }
        Ok(None)
    }

    pub fn addr(&self) -> &str { &self.addr }
}
