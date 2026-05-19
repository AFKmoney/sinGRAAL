use std::process::Command;
use std::path::PathBuf;

fn main() {
    let cuda_dir = PathBuf::from("cuda");
    let out_dir  = PathBuf::from(std::env::var("OUT_DIR").unwrap());

    // Only build CUDA if nvcc is available
    let has_nvcc = Command::new("nvcc").arg("--version").status().is_ok();

    if has_nvcc {
        let cu_src  = cuda_dir.join("kangaroo.cu");
        let obj_out = out_dir.join("kangaroo.o");
        let lib_out = out_dir.join("libkangaroo_cuda.a");

        // Detect sm version; default sm_80 (A100 / RTX 30xx)
        let arch = std::env::var("CUDA_ARCH").unwrap_or_else(|_| "sm_80".into());

        let status = Command::new("nvcc")
            .args([
                "-O3",
                &format!("-arch={}", arch),
                "--compiler-options", "-fPIC",
                "-c", cu_src.to_str().unwrap(),
                "-o", obj_out.to_str().unwrap(),
            ])
            .status()
            .expect("nvcc failed");

        assert!(status.success(), "CUDA compilation failed");

        Command::new("ar")
            .args(["rcs", lib_out.to_str().unwrap(), obj_out.to_str().unwrap()])
            .status()
            .expect("ar failed");

        println!("cargo:rustc-link-search=native={}", out_dir.display());
        println!("cargo:rustc-link-lib=static=kangaroo_cuda");
        println!("cargo:rustc-link-lib=dylib=cudart");
        println!("cargo:rustc-link-lib=dylib=stdc++");
        println!("cargo:rustc-cfg=feature=\"cuda\"");
        println!("cargo:rustc-cfg=feature_cuda");
        println!("cargo:rerun-if-changed=cuda/kangaroo.cu");
        println!("cargo:rerun-if-changed=cuda/secp256k1.cuh");
    } else {
        eprintln!("[build.rs] nvcc not found — building CPU-only fallback");
        println!("cargo:rerun-if-changed=cuda/kangaroo.cu");
    }
}
