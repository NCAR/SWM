#![allow(non_snake_case)]
//! This file outlines a typical build process which can be used for more complex CUDA projects utilising this crate.
//! It does the following:
//!     1. Use a `build.rs` file to compile your CUDA code/project into a PTX file. Your CUDA code/project can be as complicated as you need them to be, including multiple files, with headers for your struct definitions, each kernel in it's own file, etc.
//!     2. The build process compiles the kernels into a PTX file, which is written to the output directory
//!     3. The build process then uses the `bindgen` crate to generate Rust bindings for the structs defined in your CUDA code
//!     4. In the `main.rs` code, the PTX code is included as a string via the `!include_str` macro, which is then compiled using the functions in this crate (detailed in previous examples)
//!
//! The advantages of having this build process for more complex CUDA projects:
//!     - You only need to define your structs once, in your CUDA code, and the Rust bindings are generated automatically
//!     - You have full intellisense for your CUDA code since they can be stored under a separate folder or even as part of a separate project
//!
//! There are two files in this example: `main.rs` and `build.rs`. You can reference them and add to your project accordingly. The `cuda` folder in this example gives a simple example of defining structs in a separate header, including creating a `wrapper.h` header for `bindgen`

use std::time::Instant;
use cudarc::driver::*;
use cudarc::nvrtc::Ptx;

mod consts;
use consts::*;

// include the compiled PTX code as string
const CUDA_KERNEL_MY_STRUCT: &str = include_str!(concat!(env!("OUT_DIR"), "/my_struct_kernel.ptx"));

fn main() -> Result<(), DriverError> {
    // -----------------------------------------------------------------------
    // GPU Setup
    // -----------------------------------------------------------------------
    
    // setup GPU device
    let now = Instant::now();

    let ctx = CudaContext::new(0)?;
    let stream = ctx.default_stream();

    println!("Time taken to initialise CUDA: {:.2?}", now.elapsed());

    // compile ptx
    let now = Instant::now();

    // loads kernels
    let my_module = ctx.load_module(Ptx::from_src(CUDA_KERNEL_MY_STRUCT))?;
    let my_function = my_module.load_function("my_struct_kernel")?;
    let init_conds = my_module.load_function("init_conds");

    println!("Time taken to compile and load PTX: {:.2?}", now.elapsed());

    // create data
    let now = Instant::now();

    let n = TOT_LEN;
    // let my_structs = vec![MyStruct { data: [1.0; 4] }; n];
    let my_structs: Vec<f64> = vec![1.0; TOT_LEN];

    // copy to GPU
    let mut gpu_my_structs = stream.clone_htod(&my_structs)?;

    println!("Time taken to initialise data: {:.2?}", now.elapsed());

    let now = Instant::now();
    let mut launch_args = stream.launch_builder(&my_function);
    launch_args.arg(&mut gpu_my_structs);
    launch_args.arg(&n);
    let cfg = LaunchConfig::for_num_elems(n as u32);
    unsafe { launch_args.launch(cfg) }?;

    println!("Time taken to call kernel: {:.2?}", now.elapsed());

    let my_structs = stream.clone_dtoh(&gpu_my_structs)?;

    assert!(my_structs.iter().all(|&x| x == 2.0));

    // -----------------------------------------------------------------------
    // Simulation Parameters and Constants Setup
    // -----------------------------------------------------------------------

    // parameters
    let dt: f64 = 90.;
    let mut tdt: f64 = dt;

    let dx: f64 = 100000.;
    let dy: f64 = 100000.;
    let fsdx: f64 = 4. / dx;
    let fsdy: f64 = 4. / dy;

    let a: f64 = 1000000.;
    let alpha: f64 = 0.001;

    // -----------------------------------------------------------------------
    // Define Solution Arrays directly to GPU
    // -----------------------------------------------------------------------

    let mut gpu_u = stream.alloc_zeros::<f64>(TOT_LEN)?;
    let mut gpu_v = stream.alloc_zeros::<f64>(TOT_LEN)?;
    let mut gpu_p = stream.alloc_zeros::<f64>(TOT_LEN)?;
    let mut gpu_psi = stream.alloc_zeros::<f64>(TOT_LEN)?;
    let mut gpu_unew = stream.alloc_zeros::<f64>(TOT_LEN)?;
    let mut gpu_vnew = stream.alloc_zeros::<f64>(TOT_LEN)?;
    let mut gpu_pnew = stream.alloc_zeros::<f64>(TOT_LEN)?;
    let mut gpu_uold = stream.alloc_zeros::<f64>(TOT_LEN)?;
    let mut gpu_vold = stream.alloc_zeros::<f64>(TOT_LEN)?;
    let mut gpu_pold = stream.alloc_zeros::<f64>(TOT_LEN)?;
    let mut gpu_cu = stream.alloc_zeros::<f64>(TOT_LEN)?;
    let mut gpu_cv = stream.alloc_zeros::<f64>(TOT_LEN)?;
    let mut gpu_z = stream.alloc_zeros::<f64>(TOT_LEN)?;
    let mut gpu_h = stream.alloc_zeros::<f64>(TOT_LEN)?;

    // -----------------------------------------------------------------------
    // Initialize Data
    // -----------------------------------------------------------------------

    // initialize velocities u and v, pressure p
    let mut launch_kern = stream.launch_builder(&init_conds);
    launch_kern.arg(&mut gpu_u);
    launch_kern.arg(&mut gpu_v);
    launch_kern.arg(&mut gpu_p);
    launch_kern.arg(&mut gpu_psi);
    launch_kern.arg(&dx);
    launch_kern.arg(&dy);
    launch_kern.arg(&a);
    launch_kern.arg(&M_LEN);
    launch_kern.arg(&N_LEN);
    launch_kern.arg(&di);
    launch_kern.arg(&dj);
    launch_kern.arg(&pcf);
    let cfg = LaunchConfig::for_num_elems(TOT_LEN as u32);
    unsafe { launch_kern.launch(cfg) }?;

    // periodic boundary conditions
    // apply_uv_bcs(&mut u, &mut v);

    // initialize old arrays
    // for i in 0..M_LEN {
    //     for j in 0..N_LEN {
    //         uold[idx(i,j)] = u[idx(i,j)];
    //         vold[idx(i,j)] = v[idx(i,j)];
    //         pold[idx(i,j)] = p[idx(i,j)];
    //     }
    // }

    Ok(())
}
