use std::process::Command;
use std::path::PathBuf;
use std::io::Write;

fn main() {
    let cuda_dir = PathBuf::from("cuda");
    let out_dir  = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    let cuda_requested = std::env::var("CARGO_FEATURE_CUDA").is_ok();

    let has_nvcc = Command::new("nvcc").arg("--version").status().is_ok();

    if has_nvcc {
        let cu_src  = cuda_dir.join("kangaroo.cu");
        let obj_out = out_dir.join("kangaroo.o");
        let lib_out = out_dir.join("libkangaroo_cuda.a");

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

    } else if cuda_requested {
        // nvcc not found but --features cuda was requested.
        // Build a stub library so the linker doesn't fail with cryptic errors.
        // GPU functions will panic at runtime with a clear message.
        let stub_c = out_dir.join("kangaroo_stub.c");
        let stub_o = out_dir.join("kangaroo_stub.o");
        let lib_out = out_dir.join("libkangaroo_cuda.a");

        let stub_src = r#"
#include <stdio.h>
#include <stdlib.h>
static void cuda_not_available(void) {
    fprintf(stderr, "\n[sinGRAAL] ERROR: Built with --features cuda but nvcc was not found at compile time.\n");
    fprintf(stderr, "[sinGRAAL] To use GPU acceleration, install the CUDA toolkit and rebuild:\n");
    fprintf(stderr, "[sinGRAAL]   cargo build --release --features cuda\n\n");
    abort();
}
int  cuda_device_count(void){cuda_not_available();return 0;}
void cuda_set_device(int a){cuda_not_available();}
void cuda_device_name(int a,char*b,int c){cuda_not_available();}
void cuda_device_memory(int a,unsigned long long*b,unsigned long long*c){cuda_not_available();}
void*kangaroo_init(int a,int b){cuda_not_available();return 0;}
void kangaroo_set_jumps(void*a,void*b,int c){cuda_not_available();}
int  kangaroo_num_jumps(void){cuda_not_available();return 0;}
void kangaroo_launch_persistent(void*a,void*b,int c,int d,int e,unsigned long long f,unsigned int g){cuda_not_available();}
void kangaroo_terminate(void){cuda_not_available();}
void kangaroo_free(void*a){cuda_not_available();}
int  kangaroo_read_dps_live(void*a,void*b,unsigned int c){cuda_not_available();return 0;}
unsigned long long kangaroo_read_step_count(void){cuda_not_available();return 0;}
void kangaroo_update_dp_threshold(unsigned int a){cuda_not_available();}
void kangaroo_read_animals(void*a,int b,void*c){cuda_not_available();}
void kangaroo_write_animals(void*a,int b,void*c){cuda_not_available();}
"#;

        std::fs::write(&stub_c, stub_src).expect("write stub");

        let cc_ok = Command::new("cc")
            .args(["-c", stub_c.to_str().unwrap(), "-o", stub_o.to_str().unwrap()])
            .status()
            .map(|s| s.success())
            .unwrap_or(false);

        if cc_ok {
            let _ = Command::new("ar")
                .args(["rcs", lib_out.to_str().unwrap(), stub_o.to_str().unwrap()])
                .status();
            println!("cargo:rustc-link-search=native={}", out_dir.display());
            println!("cargo:rustc-link-lib=static=kangaroo_cuda");
        }

        let _ = writeln!(std::io::stderr(),
            "[build.rs] WARNING: nvcc not found — GPU functions replaced with stubs. \
             Install CUDA toolkit for real GPU support.");
        println!("cargo:rerun-if-changed=cuda/kangaroo.cu");

    } else {
        eprintln!("[build.rs] CPU-only build (no --features cuda)");
        println!("cargo:rerun-if-changed=cuda/kangaroo.cu");
    }
}
