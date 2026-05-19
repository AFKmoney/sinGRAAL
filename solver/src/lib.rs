use wasm_bindgen::prelude::*;

mod secp;
mod kangaroo;
use secp::*;
use kangaroo::*;

// ─── Literal helpers ─────────────────────────────────────────────────────────

// Literal: positive int = var true, negative = var false.
// Internal index for literal l: lit_idx(l) = 2*(|l|-1) + (if l<0 {1} else {0})
#[inline] fn lit_idx(l: i32) -> usize {
    if l > 0 { 2 * (l as usize - 1) } else { 2 * ((-l) as usize - 1) + 1 }
}
#[inline] fn var_of(l: i32) -> usize { (l.unsigned_abs() - 1) as usize }
#[inline] fn neg(l: i32) -> i32 { -l }

// ─── CDCL Solver with Two-Watched Literals ───────────────────────────────────

/// Status codes returned by step()
#[derive(PartialEq, Clone, Copy)]
enum Status { Running, Sat, Unsat }

pub struct CDCLInner {
    n: usize,                    // number of variables
    clauses: Vec<Vec<i32>>,      // all clauses (original + learned)
    watch: Vec<Vec<usize>>,      // watch[lit_idx] = clause indices watching this lit
    val: Vec<i8>,                // val[var] = 0 unset, 1 true, -1 false
    trail: Vec<i32>,             // assignment trail (literals in order)
    trail_lim: Vec<usize>,       // trail index at start of each decision level
    reason: Vec<i32>,            // reason[var] = clause_idx+1, or 0 for decisions
    level: Vec<u32>,             // decision level of each variable
    prop_q: usize,               // propagation queue head (index into trail)
    // VSIDS
    activity: Vec<f64>,
    act_bump: f64,
    // Stats
    pub decisions: u64,
    pub propagations: u64,
    pub conflicts: u64,
    pub learned: u64,
    pub backjumps: u64,
    pub max_backjump: u64,
    pub assigned: u32,
    status: Status,
}

impl CDCLInner {
    fn new(n: usize) -> Self {
        CDCLInner {
            n,
            clauses: Vec::new(),
            watch: vec![Vec::new(); 2 * n],
            val: vec![0i8; n],
            trail: Vec::new(),
            trail_lim: Vec::new(),
            reason: vec![0i32; n],
            level: vec![0u32; n],
            prop_q: 0,
            activity: vec![0.0f64; n],
            act_bump: 1.0,
            decisions: 0, propagations: 0, conflicts: 0,
            learned: 0, backjumps: 0, max_backjump: 0, assigned: 0,
            status: Status::Running,
        }
    }

    fn add_clause(&mut self, lits: &[i32]) {
        if lits.is_empty() {
            self.status = Status::Unsat;
            return;
        }
        let ci = self.clauses.len();
        let mut cl = lits.to_vec();
        if cl.len() == 1 {
            // Unit clause: force immediately
            let l = cl[0];
            let v = var_of(l);
            if self.val[v] == 0 {
                self.enqueue(l, 0);
            } else if (l > 0 && self.val[v] == -1) || (l < 0 && self.val[v] == 1) {
                self.status = Status::Unsat;
            }
            self.clauses.push(cl);
            return;
        }
        // Two-watched literals: watch first two lits
        self.watch[lit_idx(neg(cl[0]))].push(ci);
        self.watch[lit_idx(neg(cl[1]))].push(ci);
        self.clauses.push(cl);
    }

    fn enqueue(&mut self, l: i32, reason_ci: usize) {
        let v = var_of(l);
        let val = if l > 0 { 1i8 } else { -1i8 };
        self.val[v] = val;
        self.level[v] = self.trail_lim.len() as u32;
        self.reason[v] = reason_ci as i32;
        self.trail.push(l);
        self.assigned += 1;
    }

    fn lit_val(&self, l: i32) -> i8 {
        let v = var_of(l);
        let a = self.val[v];
        if a == 0 { 0 } else if l > 0 { a } else { -a }
    }

    // Two-watched literal propagation. Returns conflict clause index or None.
    fn propagate(&mut self) -> Option<usize> {
        while self.prop_q < self.trail.len() {
            let p = self.trail[self.prop_q];
            self.prop_q += 1;
            self.propagations += 1;

            // p was set true, so neg(p) is false. Update watchers of neg(p).
            let false_lit = neg(p);
            let fli = lit_idx(false_lit);

            let mut j = 0usize;
            let n_watch = self.watch[fli].len();
            let mut conflict = None;

            let mut i = 0usize;
            while i < n_watch {
                let ci = self.watch[fli][i];

                // Ensure false_lit is at position 1
                if self.clauses[ci][0] == false_lit {
                    self.clauses[ci].swap(0, 1);
                }
                // clauses[ci][0] = other watched literal
                let w0 = self.clauses[ci][0];

                // If other watched lit is true, clause satisfied
                if self.lit_val(w0) == 1 {
                    self.watch[fli][j] = ci;
                    j += 1;
                    i += 1;
                    continue;
                }

                // Try to find a new non-false literal to watch (from index 2 onward)
                let clen = self.clauses[ci].len();
                let mut new_watch: Option<(usize, i32)> = None;
                for k in 2..clen {
                    if self.lit_val(self.clauses[ci][k]) != -1 {
                        new_watch = Some((k, self.clauses[ci][k]));
                        break;
                    }
                }

                if let Some((k, nw)) = new_watch {
                    // Swap new watch into position 1
                    self.clauses[ci][1] = nw;
                    self.clauses[ci][k] = false_lit;
                    // Remove from fli watch, add to nw watch
                    self.watch[lit_idx(neg(nw))].push(ci);
                    // Don't copy to j → effectively removes ci from fli list
                    i += 1;
                } else if self.lit_val(w0) == -1 {
                    // All literals false → conflict
                    // Copy remaining watches
                    while i < n_watch {
                        self.watch[fli][j] = self.watch[fli][i];
                        j += 1;
                        i += 1;
                    }
                    self.watch[fli].truncate(j);
                    conflict = Some(ci);
                    break;
                } else {
                    // w0 is unset → unit propagation
                    self.watch[fli][j] = ci;
                    j += 1;
                    i += 1;
                    self.enqueue(w0, ci + 1); // reason = ci+1 (0 reserved for decisions)
                }
            }

            if conflict.is_none() {
                self.watch[fli].truncate(j);
            } else {
                return conflict;
            }
        }
        None
    }

    // First-UIP conflict analysis. Returns learned clause (UIP lit at index 0)
    // and the backjump level.
    fn analyze(&mut self, conflict_ci: usize) -> (Vec<i32>, u32) {
        let dl = self.trail_lim.len() as u32;
        let mut seen = vec![false; self.n];
        let mut learnt: Vec<i32> = vec![0]; // placeholder for UIP lit
        let mut counter = 0i32;             // literals at current DL to resolve
        let mut trail_pos = self.trail.len() as i32 - 1;
        let mut ci = conflict_ci;

        loop {
            // Collect vars from clause first, then bump (avoids borrow conflict)
            let clause_vars: Vec<i32> = self.clauses[ci].clone();
            for &l in &clause_vars {
                let v = var_of(l);
                if !seen[v] {
                    seen[v] = true;
                    self.bump_activity(v);
                    if self.level[v] == dl {
                        counter += 1;
                    } else if self.level[v] > 0 {
                        learnt.push(neg(l));
                    }
                }
            }

            // Walk back trail to find next seen literal at current DL
            while !seen[var_of(self.trail[trail_pos as usize])] {
                trail_pos -= 1;
            }
            let uip_lit = self.trail[trail_pos as usize];
            let v = var_of(uip_lit);
            seen[v] = false;
            trail_pos -= 1;
            counter -= 1;

            if counter == 0 {
                // Found first UIP
                learnt[0] = neg(uip_lit);
                break;
            }

            // Resolve with reason clause
            let r = self.reason[v];
            if r == 0 {
                // Decision variable — shouldn't happen before counter reaches 0, but handle
                learnt[0] = neg(uip_lit);
                break;
            }
            ci = (r - 1) as usize;
        }

        // Compute backjump level = max DL of literals in learnt (except UIP)
        let backjump_lvl = if learnt.len() == 1 {
            0
        } else {
            learnt[1..].iter()
                .map(|&l| self.level[var_of(l)])
                .max()
                .unwrap_or(0)
        };

        // Put highest-level literal at index 1 (for efficient watching)
        if learnt.len() > 1 {
            let best_pos = learnt[1..].iter().enumerate()
                .max_by_key(|(_, &l)| self.level[var_of(l)])
                .map(|(i, _)| i + 1)
                .unwrap_or(1);
            learnt.swap(1, best_pos);
        }

        (learnt, backjump_lvl)
    }

    // Backjump to given decision level
    fn backtrack(&mut self, target_dl: u32) {
        let target_idx = if target_dl as usize >= self.trail_lim.len() {
            self.trail.len()
        } else {
            self.trail_lim[target_dl as usize]
        };
        for i in target_idx..self.trail.len() {
            let v = var_of(self.trail[i]);
            self.val[v] = 0;
            self.assigned -= 1;
        }
        self.trail.truncate(target_idx);
        self.trail_lim.truncate(target_dl as usize);
        self.prop_q = self.trail.len();
    }

    fn bump_activity(&mut self, v: usize) {
        self.activity[v] += self.act_bump;
        if self.activity[v] > 1e100 {
            // Rescale
            for a in &mut self.activity { *a *= 1e-100; }
            self.act_bump *= 1e-100;
        }
    }

    fn decay_activity(&mut self) {
        self.act_bump /= 0.95; // equivalent to multiplying all activities by 0.95
    }

    // VSIDS: pick unassigned variable with highest activity
    fn pick_var(&self) -> Option<usize> {
        let mut best_v = None;
        let mut best_a = -1.0f64;
        for v in 0..self.n {
            if self.val[v] == 0 && self.activity[v] > best_a {
                best_a = self.activity[v];
                best_v = Some(v);
            }
        }
        best_v
    }

    // Run up to max_conflicts CDCL steps
    fn solve_steps(&mut self, max_conflicts: u32) -> Status {
        if self.status != Status::Running { return self.status; }

        // Initial propagation
        if self.conflicts == 0 && self.trail.len() > self.prop_q {
            if let Some(_) = self.propagate() {
                self.status = Status::Unsat;
                return self.status;
            }
        }

        let mut conflict_budget = max_conflicts;

        loop {
            if let Some(conflict_ci) = self.propagate() {
                self.conflicts += 1;
                conflict_budget = conflict_budget.saturating_sub(1);

                if self.trail_lim.is_empty() {
                    self.status = Status::Unsat;
                    return self.status;
                }

                let dl = self.trail_lim.len() as u32;
                let (learnt, backjump_lvl) = self.analyze(conflict_ci);
                let jump_dist = dl.saturating_sub(backjump_lvl);
                if u64::from(jump_dist) > self.max_backjump { self.max_backjump = u64::from(jump_dist); }
                if jump_dist > 0 { self.backjumps += 1; }

                self.backtrack(backjump_lvl);
                self.decay_activity();

                // Add learned clause and force its UIP literal
                let uip = learnt[0];
                let lci = self.clauses.len();
                if learnt.len() == 1 {
                    self.clauses.push(learnt.clone());
                    self.enqueue(uip, 0);
                } else {
                    self.watch[lit_idx(neg(learnt[0]))].push(lci);
                    self.watch[lit_idx(neg(learnt[1]))].push(lci);
                    self.clauses.push(learnt);
                    self.learned += 1;
                    self.enqueue(uip, lci + 1);
                }

                if conflict_budget == 0 { return Status::Running; }
            } else {
                // No conflict — pick next decision
                match self.pick_var() {
                    None => {
                        self.status = Status::Sat;
                        return Status::Sat;
                    }
                    Some(v) => {
                        self.decisions += 1;
                        self.trail_lim.push(self.trail.len());
                        // Decide: prefer true if activity is tied, else use activity sign
                        let l = if self.activity[v] >= 0.0 { (v + 1) as i32 } else { -((v + 1) as i32) };
                        self.enqueue(l, 0);
                    }
                }
            }
        }
    }
}

// ─── WASM Public API ─────────────────────────────────────────────────────────

#[wasm_bindgen]
pub struct WasmSolver {
    inner: CDCLInner,
    finalized: bool,
}

#[wasm_bindgen]
impl WasmSolver {
    #[wasm_bindgen(constructor)]
    pub fn new(num_vars: u32) -> WasmSolver {
        WasmSolver {
            inner: CDCLInner::new(num_vars as usize),
            finalized: false,
        }
    }

    /// Add a clause. Pass literals as a flat JS array, e.g. [1, -2, 3]
    pub fn add_clause(&mut self, lits: &[i32]) {
        self.inner.add_clause(lits);
    }

    /// Add multiple clauses packed as flat array with 0 separators.
    /// e.g. [1, -2, 0, 3, 4, 0] = two clauses [1,-2] and [3,4]
    pub fn add_clauses_packed(&mut self, data: &[i32]) {
        let mut start = 0;
        for i in 0..=data.len() {
            if i == data.len() || data[i] == 0 {
                if i > start {
                    self.inner.add_clause(&data[start..i]);
                }
                start = i + 1;
            }
        }
    }

    /// Run up to max_conflicts CDCL steps. Returns JSON status string.
    pub fn step(&mut self, max_conflicts: u32) -> String {
        let s = self.inner.solve_steps(max_conflicts);
        let status = match s {
            Status::Running => "running",
            Status::Sat     => "sat",
            Status::Unsat   => "unsat",
        };
        format!(
            r#"{{"status":"{status}","decisions":{d},"propagations":{p},"conflicts":{c},"learned":{l},"backjumps":{b},"maxBackjump":{mb},"assigned":{a},"numVars":{n}}}"#,
            status = status,
            d = self.inner.decisions,
            p = self.inner.propagations,
            c = self.inner.conflicts,
            l = self.inner.learned,
            b = self.inner.backjumps,
            mb = self.inner.max_backjump,
            a = self.inner.assigned,
            n = self.inner.n,
        )
    }

    /// After SAT: get assignment as flat i8 array indexed by var-1.
    /// val[i] = 1 (true), -1 (false), 0 (unset)
    pub fn get_assignment(&self) -> Vec<i8> {
        self.inner.val.clone()
    }

    /// Number of clauses (including learned)
    pub fn clause_count(&self) -> u32 {
        self.inner.clauses.len() as u32
    }

    /// Number of assigned variables
    pub fn assigned_count(&self) -> u32 {
        self.inner.assigned
    }

    /// Return all unit-forced literals as JSON array [{var, val}].
    /// These are bits of k that CDCL has proven must have a specific value.
    pub fn get_forced_bits(&self, k_lits: &[i32]) -> String {
        let mut forced = Vec::new();
        for (bit_idx, &lit) in k_lits.iter().enumerate() {
            let v = (lit.unsigned_abs() - 1) as usize;
            if v < self.inner.n && self.inner.val[v] != 0 {
                // Check this was forced at level 0 (invariant, not guessed)
                if self.inner.level[v] == 0 {
                    let val = if lit > 0 { self.inner.val[v] == 1 }
                              else       { self.inner.val[v] == -1 };
                    forced.push(format!(r#"{{"bit":{bit},"val":{val}}}"#,
                        bit = bit_idx,
                        val = if val { 1 } else { 0 }));
                }
            }
        }
        format!("[{}]", forced.join(","))
    }

    /// Number of bits proven at level-0 (invariant constraints from CDCL)
    pub fn proven_bit_count(&self, k_lits: &[i32]) -> u32 {
        k_lits.iter().filter(|&&lit| {
            let v = (lit.unsigned_abs() - 1) as usize;
            v < self.inner.n && self.inner.val[v] != 0 && self.inner.level[v] == 0
        }).count() as u32
    }
}

// ─── WASM Kangaroo Engine ─────────────────────────────────────────────────────

#[wasm_bindgen]
pub struct WasmKangaroo {
    engine: KangarooEngine,
}

#[wasm_bindgen]
impl WasmKangaroo {
    /// Create a new Kangaroo for ECDLP: find k such that k·G = (target_x, target_y)
    /// range_low, target_x, target_y: 64-char hex strings (no 0x prefix)
    /// range_bits: 2^range_bits is the search range size
    /// dp_bits: number of leading zero bits for distinguished points (higher = fewer DPs, faster per DP)
    #[wasm_bindgen(constructor)]
    pub fn new(
        target_x: &str,
        target_y: &str,
        range_low: &str,
        range_bits: u32,
        num_tames: u32,
        dp_bits: u32,
    ) -> Result<WasmKangaroo, JsValue> {
        let tx = fe_from_hex(target_x).ok_or("bad target_x")?;
        let ty = fe_from_hex(target_y).ok_or("bad target_y")?;
        let rl = fe_from_hex(range_low).ok_or("bad range_low")?;
        let target = Pt { x: tx, y: ty, inf: false };
        let engine = KangarooEngine::new(target, rl, range_bits, num_tames as usize, dp_bits);
        Ok(WasmKangaroo { engine })
    }

    /// Step the engine by `iterations` per animal per call.
    /// Returns JSON: {status, steps, dps, key, wildX, tameX}
    pub fn step(&mut self, iterations: u32) -> String {
        let found = self.engine.step(iterations);
        let status = if found { "found" } else { "running" };
        let key = self.engine.solution
            .map(|k| fe_to_hex(k))
            .unwrap_or_default();
        format!(
            r#"{{"status":"{s}","steps":{st},"dps":{dp},"key":"{k}","wildX":"{wx}","tameX":"{tx}"}}"#,
            s  = status,
            st = self.engine.steps,
            dp = self.engine.dps_found,
            k  = key,
            wx = fe_to_hex(self.engine.last_wild_x),
            tx = fe_to_hex(self.engine.last_tame_x),
        )
    }

    /// Show the 6-orbit for a given point (educational / visualization)
    /// Returns JSON array: [{name, x, y}, ...]
    pub fn orbit_demo(x: &str, y: &str) -> Result<String, JsValue> {
        let fx = fe_from_hex(x).ok_or("bad x")?;
        let fy = fe_from_hex(y).ok_or("bad y")?;
        let p = Pt { x: fx, y: fy, inf: false };
        let orb = orbit_demo(p);
        let items: Vec<String> = orb.iter().map(|(name, ox, oy)| {
            format!(r#"{{"name":"{}","x":"{}","y":"{}"}}"#, name, ox, oy)
        }).collect();
        Ok(format!("[{}]", items.join(",")))
    }

    /// Verify that k·G = target (sanity check after solve)
    pub fn verify(k_hex: &str, target_x: &str, target_y: &str) -> bool {
        let k  = match fe_from_hex(k_hex)   { Some(v) => v, None => return false };
        let tx = match fe_from_hex(target_x) { Some(v) => v, None => return false };
        let ty = match fe_from_hex(target_y) { Some(v) => v, None => return false };
        let pt = scalar_mul(gen(), k);
        !pt.inf && pt.x == tx && pt.y == ty
    }

    /// Compute k·G and return hex (x, y) – useful for testing
    pub fn scalar_mul_hex(k_hex: &str) -> Result<String, JsValue> {
        let k = fe_from_hex(k_hex).ok_or("bad k")?;
        let pt = scalar_mul(gen(), k);
        if pt.inf { return Ok(r#"{"inf":true}"#.to_string()); }
        Ok(format!(r#"{{"x":"{}","y":"{}"}}"#, fe_to_hex(pt.x), fe_to_hex(pt.y)))
    }

    pub fn steps(&self) -> f64 { self.engine.steps as f64 }
    pub fn dps(&self) -> u32   { self.engine.dps_found as u32 }
    pub fn solved(&self) -> bool { self.engine.solved }
}
