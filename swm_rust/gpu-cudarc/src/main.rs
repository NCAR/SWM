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

use std::mem;
use std::fs::File;
use std::path::Path;
use std::io::prelude::*;

mod constants;
use constants::*;

// include the compiled PTX code as string
const CUDA_KERNEL: &str = include_str!(concat!(env!("OUT_DIR"), "/kernels.ptx"));

fn main() -> Result<(), DriverError> {
    // -----------------------------------------------------------------------
    // GPU Setup
    // -----------------------------------------------------------------------

    let ctx = CudaContext::new(0)?;
    let stream = ctx.default_stream();

    // loads kernels
    let my_module = ctx.load_module(Ptx::from_src(CUDA_KERNEL))?;
    let init_conds = my_module.load_function("init_conds")?;
    let apply_uv_bcs = my_module.load_function("apply_uv_bcs")?;
    let init_olds = my_module.load_function("init_olds")?;
    let update_intermed_vars = my_module.load_function("update_intermed_vars")?;
    let apply_intermed_bcs = my_module.load_function("apply_intermed_bcs")?;
    let time_update_new_vars = my_module.load_function("time_update_new_vars")?;
    let apply_uvp_bcs = my_module.load_function("apply_uvp_bcs")?;
    let smooth_update_old_vars = my_module.load_function("smooth_update_old_vars")?;

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
    // Define Arrays directly to GPU
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

    // solution arrays on CPU
    let mut u: Vec<f64> = vec![0.0; TOT_LEN];
    let mut v: Vec<f64> = vec![0.0; TOT_LEN];
    let mut p: Vec<f64> = vec![0.0; TOT_LEN];


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
    let mnmax = M.max(N);
    let cfg = LaunchConfig::for_num_elems(mnmax as u32);
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

    // Start timer
    let tstart = Instant::now(); 
    let mut time = 0.;
    let mut t100 = 0.;
    let mut t200 = 0.;
    let mut t300 = 0.;

    // -----------------------------------------------------------------------
    // Time Marching Loop
    // -----------------------------------------------------------------------
    
    for ncycle in 1..=ITMAX {

        let mut c1 = tstart.elapsed().as_secs_f64();

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
        launch_kern.arg(&N_LEN);
        launch_kern.arg(&TOT_LEN);
        let cfg = LaunchConfig::for_num_elems(TOT_LEN as u32);
        unsafe { launch_kern.launch(cfg) }?;

        let mut c2 = tstart.elapsed().as_secs_f64();
        t100 = t100 + (c2 - c1);

        // apply periodic boundary conditions to intermediate variables
        let mut launch_kern = stream.launch_builder(&apply_intermed_bcs);
        launch_kern.arg(&mut gpu_cu);
        launch_kern.arg(&mut gpu_cv);
        launch_kern.arg(&mut gpu_z);
        launch_kern.arg(&mut gpu_h);
        launch_kern.arg(&M);
        launch_kern.arg(&N);
        launch_kern.arg(&N_LEN);
        let cfg = LaunchConfig::for_num_elems(mnmax as u32);
        unsafe { launch_kern.launch(cfg) }?;

        // time update to new variables
        let tdts8 = tdt / 8.0;
        let tdtsdx = tdt / dx;
        let tdtsdy = tdt / dy;

        c1 = tstart.elapsed().as_secs_f64();

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

        c2 = tstart.elapsed().as_secs_f64();
        t200 = t200 + (c2 - c1);
        
        // apply periodic boundary conitions to new variables
        let mut launch_kern = stream.launch_builder(&apply_uvp_bcs);
        launch_kern.arg(&mut gpu_unew);
        launch_kern.arg(&mut gpu_vnew);
        launch_kern.arg(&mut gpu_pnew);
        launch_kern.arg(&M);
        launch_kern.arg(&N);
        launch_kern.arg(&N_LEN);
        let cfg = LaunchConfig::for_num_elems(mnmax as u32);
        unsafe { launch_kern.launch(cfg) }?;

        // update time
        time = time + dt;

        // update update old vars and solution
        if ncycle > 1 {

            c1 = tstart.elapsed().as_secs_f64();

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

            c2 = tstart.elapsed().as_secs_f64(); 
            t300 = t300 + (c2 - c1);
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

    // copy solutions over to CPU
    let mut c1 = tstart.elapsed().as_secs_f64();
    stream.memcpy_dtoh(&gpu_u, &mut u)?;
    stream.memcpy_dtoh(&gpu_v, &mut v)?;
    stream.memcpy_dtoh(&gpu_p, &mut p)?;
    let mut c2 = tstart.elapsed().as_secs_f64();
    let cpy_time = c2 - c1;
    println!("memcpy time: {:?}",cpy_time);

    c1 = tstart.elapsed().as_secs_f64();
    let u = stream.clone_dtoh(&gpu_u)?;
    let v = stream.clone_dtoh(&gpu_v)?;
    let p = stream.clone_dtoh(&gpu_p)?;
    c2 = tstart.elapsed().as_secs_f64();
    let clone_time = c2 - c1;
    println!("clone time: {:?}",clone_time);

    // End time
    let ctime = tstart.elapsed().as_secs_f64();

    let ptime = time / 3600.;

    // print out timings
    if TIMING {
        println!(" cycle number {:?} model time in hours {:?}\n", ITMAX, ptime);

        let mut mfs100 = 0.0;
        let mut mfs200 = 0.0;
        let mut mfs300 = 0.0;
        // gdr t100 etc. now an accumulation of all l100 time
        if t100 > 0. { mfs100 = ITMAX as f64 * 24. * M as f64 * N as f64 / t100 / 1000000.; }
        if t200 > 0. { mfs200 = ITMAX as f64 * 26. * M as f64 * N as f64 / t200 / 1000000.; }
        if t300 > 0. { mfs300 = ITMAX as f64 * 15. * M as f64 * N as f64 / t300 / 1000000.; }

        let tcyc = ctime / ITMAX as f64;

        println!(" cycle number {:?} total computer time {:?} time per cycle {:?}", ITMAX, ctime, tcyc);
        println!(" time and megaflops for loop 100 {:?} {:?}", t100, mfs100);
        println!(" time and megaflops for loop 200 {:?} {:?}", t200, mfs200);
        println!(" time and megaflops for loop 300 {:?} {:?}", t300, mfs300);
    }

    // save solutions to txt files
    if VAL_OUT {
        print_data_to_file("u_rust.txt", &u);
        print_data_to_file("v_rust.txt", &v);
        print_data_to_file("p_rust.txt", &p);
    }

    if SUCCINCT {
        println!("Version: GPU cudarc");
        println!("Grid Size: {:?}x{:?}, Number of iterations: {:?}, Total computer time: {:.2}", M, N, ITMAX, ctime);
    }

    Ok(())
}

fn print_data_to_file(pathname: &str, data: &Vec<f64>) {
    // define path and display
    let path = Path::new(pathname);
    let display = path.display();

    // Open a file in write-only mode, returns `io::Result<File>`
    let mut file = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", display, why),
        Ok(file) => file,
    };

    // Create string from data
    let mut s = String::from("");
    for i in 0..M_LEN {
        for j in 0..N_LEN {
            s += &format!("{:.6} ", data[i * N_LEN + j]);
        }
        s += "\n";
    }

    // Write string to `file`, returns `io::Result<()>`
    match file.write_all(s.as_bytes()) {
        Err(why) => panic!("couldn't write to {}: {}", display, why),
        Ok(_) => println!("successfully wrote to {}", display),
    }
}