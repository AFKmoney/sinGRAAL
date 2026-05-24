// Build script: compile CUDA kernel with Toom-Cook-6 field multiplication.
//
// Requirements: nvcc on PATH, CUDA toolkit installed.
// Arch: sm_80 (Ampere A100/A10), sm_86 (RTX 3090), sm_89 (Ada L4/L40).
//       Override via CUDA_ARCH env var, e.g. CUDA_ARCH=sm_86 cargo build
//
// Toom-6 is compiled by default.  Add -DUSE_SCHOOLBOOK to compare with
// the standard 16-MAD schoolbook baseline.

use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=cuda/kangaroo_toom6.cu");
    println!("cargo:rerun-if-changed=cuda/secp256k1_toom6.cuh");

    let cuda_available = Command::new("nvcc").arg("--version").status()
        .map(|s| s.success()).unwrap_or(false);

    if !cuda_available {
        eprintln!("cargo:warning=nvcc not found — CUDA disabled, CPU-only build");
        return;
    }

    let arch  = env::var("CUDA_ARCH").unwrap_or_else(|_| "sm_80".to_string());
    let out   = PathBuf::from(env::var("OUT_DIR").unwrap());
    let obj   = out.join("kangaroo_toom6.o");
    let lib   = out.join("libkangaroo_toom6.a");

    let status = Command::new("nvcc")
        .args([
            "-O3",
            &format!("-arch={arch}"),
            "--compiler-options", "-fPIC",
            "-c", "cuda/kangaroo_toom6.cu",
            "-o", obj.to_str().unwrap(),
        ])
        .status()
        .expect("nvcc compilation failed");

    assert!(status.success(), "nvcc exited with non-zero status");

    let ar_status = Command::new("ar")
        .args(["rcs", lib.to_str().unwrap(), obj.to_str().unwrap()])
        .status()
        .expect("ar failed");

    assert!(ar_status.success(), "ar exited with non-zero status");

    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=kangaroo_toom6");
    println!("cargo:rustc-link-lib=cudart");
    println!("cargo:rustc-cfg=feature=\"cuda\"");

    // Link CUDA runtime from detected CUDA path
    if let Ok(cuda_path) = env::var("CUDA_PATH") {
        println!("cargo:rustc-link-search=native={}/lib64", cuda_path);
    } else {
        println!("cargo:rustc-link-search=native=/usr/local/cuda/lib64");
    }
}
