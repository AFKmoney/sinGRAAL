// Build script: compile the Toom-Cook-6 CUDA kernel with AUTOMATIC GPU detection.
//
// Architecture selection (in priority order):
//   1. CUDA_ARCH env var set (e.g. CUDA_ARCH=sm_89) → single-arch build (fastest
//      compile, smallest binary). Use this when you know the exact target GPU.
//   2. No CUDA_ARCH, but nvidia-smi present on the build host → detect the local
//      GPU's compute capability and build natively for it (+ PTX for JIT).
//   3. No CUDA_ARCH and no GPU on the build host (typical cloud/CI image build)
//      → build a MULTI-ARCH FATBINARY covering every common datacenter/consumer
//      GPU. The CUDA driver picks the matching cubin at load time, and the
//      embedded PTX JIT-compiles for any future architecture. One image runs
//      everywhere: T4, V100, A100, A10, RTX 30xx, RTX 40xx, L4/L40, H100, …
//
// Requirements: nvcc on PATH, CUDA toolkit installed.
//
// Toom-6 is compiled by default. Add -DUSE_SCHOOLBOOK to compare with the
// standard 16-MAD schoolbook baseline.

use std::env;
use std::path::PathBuf;
use std::process::Command;

// Common architectures covering Volta → Hopper.
// compute_XX,code=sm_XX → native cubin;  the final compute_90 also emits PTX.
const FATBIN_ARCHS: &[u32] = &[70, 75, 80, 86, 89, 90];

/// Query the local GPU compute capability via nvidia-smi, e.g. "8.9" → 89.
fn detect_local_cc() -> Option<u32> {
    let out = Command::new("nvidia-smi")
        .args(["--query-gpu=compute_cap", "--format=csv,noheader,nounits"])
        .output()
        .ok()?;
    if !out.status.success() {
        return None;
    }
    let text = String::from_utf8_lossy(&out.stdout);
    // Take the first GPU's capability (e.g. "8.9").
    let first = text.lines().next()?.trim();
    let digits: String = first.chars().filter(|c| c.is_ascii_digit()).collect();
    digits.parse::<u32>().ok()
}

/// Build the nvcc -gencode flags for the chosen strategy.
fn gencode_flags() -> (Vec<String>, String) {
    // 1. Explicit override (ignore empty string, e.g. ENV CUDA_ARCH="").
    if let Ok(arch) = env::var("CUDA_ARCH").map(|s| s.trim().to_string()) {
      if !arch.is_empty() {
        let cc: String = arch.trim_start_matches("sm_").to_string();
        let flags = vec![
            "-gencode".to_string(),
            format!("arch=compute_{cc},code=sm_{cc}"),
            // Embed PTX of the same arch for forward-compat JIT.
            "-gencode".to_string(),
            format!("arch=compute_{cc},code=compute_{cc}"),
        ];
        return (flags, format!("explicit CUDA_ARCH=sm_{cc} (+PTX)"));
      }
    }

    // 2. Detect the build host's GPU.
    if let Some(cc) = detect_local_cc() {
        let flags = vec![
            "-gencode".to_string(),
            format!("arch=compute_{cc},code=sm_{cc}"),
            "-gencode".to_string(),
            format!("arch=compute_{cc},code=compute_{cc}"),
        ];
        return (flags, format!("detected local GPU sm_{cc} (+PTX)"));
    }

    // 3. No GPU on the build host → portable fatbinary.
    let mut flags = Vec::new();
    for &cc in FATBIN_ARCHS {
        flags.push("-gencode".to_string());
        flags.push(format!("arch=compute_{cc},code=sm_{cc}"));
    }
    // Embed PTX of the highest arch so newer GPUs JIT-compile at load time.
    let top = FATBIN_ARCHS.last().copied().unwrap_or(90);
    flags.push("-gencode".to_string());
    flags.push(format!("arch=compute_{top},code=compute_{top}"));
    let archs: Vec<String> = FATBIN_ARCHS.iter().map(|c| format!("sm_{c}")).collect();
    (flags, format!("portable fatbinary [{}] (+compute_{top} PTX)", archs.join(",")))
}

fn main() {
    println!("cargo:rerun-if-changed=cuda/kangaroo_toom6.cu");
    println!("cargo:rerun-if-changed=cuda/secp256k1_toom6.cuh");
    println!("cargo:rerun-if-env-changed=CUDA_ARCH");

    let cuda_available = Command::new("nvcc")
        .arg("--version")
        .status()
        .map(|s| s.success())
        .unwrap_or(false);

    if !cuda_available {
        println!("cargo:warning=nvcc not found — CUDA disabled, CPU-only build");
        return;
    }

    let (gencode, strategy) = gencode_flags();
    println!("cargo:warning=CUDA build strategy: {strategy}");

    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let obj = out.join("kangaroo_toom6.o");
    let lib = out.join("libkangaroo_toom6.a");

    let mut nvcc = Command::new("nvcc");
    nvcc.arg("-O3");
    for f in &gencode {
        nvcc.arg(f);
    }
    nvcc.args([
        "--compiler-options",
        "-fPIC",
        "-c",
        "cuda/kangaroo_toom6.cu",
        "-o",
        obj.to_str().unwrap(),
    ]);

    let status = nvcc.status().expect("nvcc compilation failed to launch");
    assert!(status.success(), "nvcc exited with non-zero status");

    let ar_status = Command::new("ar")
        .args(["rcs", lib.to_str().unwrap(), obj.to_str().unwrap()])
        .status()
        .expect("ar failed to launch");
    assert!(ar_status.success(), "ar exited with non-zero status");

    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=kangaroo_toom6");
    println!("cargo:rustc-link-lib=cudart");
    println!("cargo:rustc-cfg=feature=\"cuda\"");

    // Link the CUDA runtime from the detected toolkit path.
    if let Ok(cuda_path) = env::var("CUDA_PATH") {
        println!("cargo:rustc-link-search=native={cuda_path}/lib64");
    } else {
        println!("cargo:rustc-link-search=native=/usr/local/cuda/lib64");
    }
}
