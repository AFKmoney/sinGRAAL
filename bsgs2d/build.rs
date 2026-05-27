// sinGRAAL bsgs2d — build.rs
// Compile CUDA filter if the "cuda" feature is enabled.

fn main() {
    #[cfg(feature = "cuda")]
    compile_cuda();
}

#[cfg(feature = "cuda")]
fn compile_cuda() {
    use std::process::Command;
    use std::path::PathBuf;

    let cuda_arch = std::env::var("CUDA_ARCH").unwrap_or_else(|_| "sm_80".into());
    let out_dir   = std::env::var("OUT_DIR").expect("OUT_DIR not set");

    let cu_src = "cuda/direct_filter.cu";
    let obj    = PathBuf::from(&out_dir).join("direct_filter.o");
    let lib    = PathBuf::from(&out_dir).join("libs3filter.a");

    println!("cargo:rerun-if-changed={cu_src}");

    // Compile .cu → .o
    let status = Command::new("nvcc")
        .args([
            "-O3",
            "--generate-code", &format!("arch=compute_{},code={}", &cuda_arch[3..], cuda_arch),
            "-Xcompiler", "-fPIC",
            "--std=c++17",
            "-c", cu_src,
            "-o", obj.to_str().unwrap(),
        ])
        .status()
        .expect("nvcc not found — install CUDA toolkit or remove --features cuda");

    if !status.success() {
        panic!("nvcc compilation failed for {cu_src}");
    }

    // Archive .o → static lib
    Command::new("ar")
        .args(["rcs", lib.to_str().unwrap(), obj.to_str().unwrap()])
        .status()
        .expect("ar failed");

    println!("cargo:rustc-link-search=native={out_dir}");
    println!("cargo:rustc-link-lib=static=s3filter");
    println!("cargo:rustc-link-lib=dylib=cuda");
    println!("cargo:rustc-link-lib=dylib=cudart");
}
