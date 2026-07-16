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
use std::f64::consts;

mod constants;
use constants::*;

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
    let init_conds = my_module.load_function("init_conds")?;
    let apply_uv_bcs = my_module.load_function("apply_uv_bcs")?;
    let init_olds = my_module.load_function("init_olds")?;
    let update_intermed_vars = my_module.load_function("update_intermed_vars")?;
    let apply_intermed_bcs = my_module.load_function("apply_intermed_bcs")?;
    let time_update_new_vars = my_module.load_function("time_update_new_vars")?;
    let apply_uvp_bcs = my_module.load_function("apply_uvp_bcs")?;
    let smooth_update_old_vars = my_module.load_function("smooth_update_old_vars")?;

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

    let mut u: Vec<f64> = vec![0.0; TOT_LEN];
    let mut v: Vec<f64> = vec![0.0; TOT_LEN];
    let mut p: Vec<f64> = vec![0.0; TOT_LEN];

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

    // set params
    let el: f64 = N as f64 * dx;
    let pi = consts::PI;
    let tpi: f64 = pi + pi;
    let di: f64 = tpi / M as f64;
    let dj: f64 = tpi / N as f64;
    let pcf: f64 = pi * pi * a * a / (el * el);

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
    launch_kern.arg(&TOT_LEN);
    launch_kern.arg(&di);
    launch_kern.arg(&dj);
    launch_kern.arg(&pcf);
    let cfg = LaunchConfig::for_num_elems(TOT_LEN as u32);
    unsafe { launch_kern.launch(cfg) }?;

    // periodic boundary conditions
    let mut launch_kern = stream.launch_builder(&apply_uv_bcs);
    launch_kern.arg(&mut gpu_u);
    launch_kern.arg(&mut gpu_v);
    launch_kern.arg(&M);
    launch_kern.arg(&N);
    launch_kern.arg(&N_LEN);
    let mnmin = M.min(N);
    let cfg = LaunchConfig::for_num_elems(mnmin as u32);
    unsafe { launch_kern.launch(cfg) }?;

    // initialize old arrays
    let mut launch_kern = stream.launch_builder(&init_olds);
    launch_kern.arg(&mut gpu_uold);
    launch_kern.arg(&mut gpu_vold);
    launch_kern.arg(&mut gpu_pold);
    launch_kern.arg(&gpu_u);
    launch_kern.arg(&gpu_v);
    launch_kern.arg(&gpu_p);
    launch_kern.arg(&TOT_LEN);
    let cfg = LaunchConfig::for_num_elems(TOT_LEN as u32);
    unsafe { launch_kern.launch(cfg) }?;

    let mut time = 0.;

    // -----------------------------------------------------------------------
    // Time Marching Loop
    // -----------------------------------------------------------------------
    
    for ncycle in 1..=ITMAX {

        // let mut c1 = tstart.elapsed().as_secs_f64();

        // compute intermediate variables cu, cv, z, and h using u, v, and p
        let mut launch_kern = stream.launch_builder(&update_intermed_vars);
        launch_kern.arg(&gpu_u);
        launch_kern.arg(&gpu_v);
        launch_kern.arg(&gpu_p);
        launch_kern.arg(&fsdx);
        launch_kern.arg(&fsdy);
        launch_kern.arg(&mut gpu_cu);
        launch_kern.arg(&mut gpu_cv);
        launch_kern.arg(&mut gpu_z);
        launch_kern.arg(&mut gpu_h);
        launch_kern.arg(&TOT_LEN);
        let cfg = LaunchConfig::for_num_elems(TOT_LEN as u32);
        unsafe { launch_kern.launch(cfg) }?;

        // let mut c2 = tstart.elapsed().as_secs_f64();
        // t100 = t100 + (c2 - c1);

        // apply periodic boundary conditions to intermediate variables
        let mut launch_kern = stream.launch_builder(&apply_intermed_bcs);
        launch_kern.arg(&mut gpu_cu);
        launch_kern.arg(&mut gpu_cv);
        launch_kern.arg(&mut gpu_z);
        launch_kern.arg(&mut gpu_h);
        launch_kern.arg(&M);
        launch_kern.arg(&N);
        launch_kern.arg(&N_LEN);
        let cfg = LaunchConfig::for_num_elems(mnmin as u32);
        unsafe { launch_kern.launch(cfg) }?;

        // time update to new variables
        let tdts8 = tdt / 8.0;
        let tdtsdx = tdt / dx;
        let tdtsdy = tdt / dy;

        // c1 = tstart.elapsed().as_secs_f64();

        let mut launch_kern = stream.launch_builder(&time_update_new_vars);
        launch_kern.arg(&gpu_uold);
        launch_kern.arg(&gpu_vold);
        launch_kern.arg(&gpu_pold);
        launch_kern.arg(&gpu_cu);
        launch_kern.arg(&gpu_cv);
        launch_kern.arg(&gpu_z);
        launch_kern.arg(&gpu_h);
        launch_kern.arg(&tdts8);
        launch_kern.arg(&tdtsdx);
        launch_kern.arg(&tdtsdy);
        launch_kern.arg(&mut gpu_unew);
        launch_kern.arg(&mut gpu_vnew);
        launch_kern.arg(&mut gpu_pnew);
        launch_kern.arg(&N_LEN);
        launch_kern.arg(&TOT_LEN);
        let cfg = LaunchConfig::for_num_elems(TOT_LEN as u32);
        unsafe { launch_kern.launch(cfg) }?;

        // c2 = tstart.elapsed().as_secs_f64();
        // t200 = t200 + (c2 - c1);
        
        // apply periodic boundary conitions to new variables
        let mut launch_kern = stream.launch_builder(&apply_uvp_bcs);
        launch_kern.arg(&mut gpu_u);
        launch_kern.arg(&mut gpu_v);
        launch_kern.arg(&mut gpu_p);
        launch_kern.arg(&M);
        launch_kern.arg(&N);
        launch_kern.arg(&N_LEN);
        let cfg = LaunchConfig::for_num_elems(mnmin as u32);
        unsafe { launch_kern.launch(cfg) }?;

        // update time
        time = time + dt;

        // update update old vars and solution
        if ncycle > 1 {

            // c1 = tstart.elapsed().as_secs_f64();

            // smooth old vars using time filter
            let mut launch_kern = stream.launch_builder(&smooth_update_old_vars);
            launch_kern.arg(&gpu_u);
            launch_kern.arg(&gpu_v);
            launch_kern.arg(&gpu_p);
            launch_kern.arg(&gpu_unew);
            launch_kern.arg(&gpu_vnew);
            launch_kern.arg(&gpu_pnew);
            launch_kern.arg(&mut gpu_uold);
            launch_kern.arg(&mut gpu_vold);
            launch_kern.arg(&mut gpu_pold);
            launch_kern.arg(&alpha);
            launch_kern.arg(&N_LEN);
            launch_kern.arg(&TOT_LEN);
            let cfg = LaunchConfig::for_num_elems(TOT_LEN as u32);
            unsafe { launch_kern.launch(cfg) }?;

            // update u, v, and p to new solution
            mem::swap(&mut gpu_u, &mut gpu_unew);
            mem::swap(&mut gpu_v, &mut gpu_vnew);
            mem::swap(&mut gpu_p, &mut gpu_pnew);

            // c2 = tstart.elapsed().as_secs_f64(); 
            // t300 = t300 + (c2 - c1);
        } else {
            // update tdt for subsequent timesteps
            tdt = tdt + tdt;

            // no smoothing for first timestep
            // this might be redundant
            mem::swap(&mut gpu_uold, &mut gpu_u);
            mem::swap(&mut gpu_vold, &mut gpu_v);
            mem::swap(&mut gpu_pold, &mut gpu_p);

            // update u, v, and p to new solution
            // might be able to take out of if statement
            mem::swap(&mut gpu_u, &mut gpu_unew);
            mem::swap(&mut gpu_v, &mut gpu_vnew);
            mem::swap(&mut gpu_p, &mut gpu_pnew);
        }
    }

    let u = stream.clone_dtoh(&gpu_u)?;
    let v = stream.clone_dtoh(&gpu_v)?;
    let p = stream.clone_dtoh(&gpu_p)?;

    // save solutions to txt files
    if VAL_OUT {
        print_data_to_file("u_rust.txt", &u);
        print_data_to_file("v_rust.txt", &v);
        print_data_to_file("p_rust.txt", &p);
    }

    Ok(())
}
