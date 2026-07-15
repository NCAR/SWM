use std::{env, path::PathBuf, process::Command};

fn main() {
    // Tell cargo to invalidate the built crate whenever files of interest changes.
    println!("cargo:rerun-if-changed={}", "cuda");

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());

    // Specify the desired architecture version.
    let arch = "compute_80"; // For example, using SM 8.6 (Ampere architecture).
    let code = "sm_80"; // For the same SM 8.6 (Ampere architecture).

    // build the cuda kernels
    let cuda_src = PathBuf::from("src/cuda/kernels/my_struct_kernel.cu");
    let ptx_file = out_dir.join("my_struct_kernel.ptx");

    let nvcc_status = Command::new("nvcc")
        .arg("-ptx")
        .arg("-o")
        .arg(&ptx_file)
        .arg(&cuda_src)
        .arg(format!("-arch={}", arch))
        .arg(format!("-code={}", code))
        .status()
        .unwrap();

    assert!(
        nvcc_status.success(),
        "Failed to compile CUDA source to PTX."
    );
}
