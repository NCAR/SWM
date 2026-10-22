// Timing and environment stuff 
use std::time::Instant;
// use std::env;

use std::mem;

// Utils is a helper module that contains some utility functions in src/utils.rs
mod consts;
mod types;
mod utils;
use consts::*;
use types::{Arr,idx,make_arr};
use utils::*;

fn main() {

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
    // Define Solution Arrays
    // -----------------------------------------------------------------------

    let mut u: Arr = make_arr();
    let mut v: Arr = make_arr();
    let mut p: Arr = make_arr();
    let mut unew: Arr = make_arr();
    let mut vnew: Arr = make_arr();
    let mut pnew: Arr = make_arr();
    let mut uold: Arr = make_arr();
    let mut vold: Arr = make_arr();
    let mut pold: Arr = make_arr();
    let mut cu: Arr = make_arr();
    let mut cv: Arr = make_arr();
    let mut z: Arr = make_arr();
    let mut h: Arr = make_arr();

    // -----------------------------------------------------------------------
    // Initialize Data
    // -----------------------------------------------------------------------

    // initialize velocities u and v, pressure p
    init_conds(&mut u, &mut v, &mut p, dx, dy, a);

    // periodic boundary conditions
    apply_uv_bcs(&mut u, &mut v);

    // initialize old arrays
    for i in 0..M_LEN {
        for j in 0..N_LEN {
            uold[idx(i,j)] = u[idx(i,j)];
            vold[idx(i,j)] = v[idx(i,j)];
            pold[idx(i,j)] = p[idx(i,j)];
        }
    }

    // print out initial values
    if VERBOSE {
        println!(" number of points in the x direction {:}", N); 
        println!(" number of points in the y direction {:}", M); 
        println!(" grid spacing in the x direction     {:}", dx); 
        println!(" grid spacing in the y direction     {:}", dy); 
        println!(" time step                           {:}", dt); 
        println!(" time filter parameter               {:}", alpha);

        let mnmin = M.min(N);
        println!(" initial diagonal elements of p");
        for i in 0..mnmin {
          print!("{:.6} ",p[idx(i,i)]);
        }
        println!("\n initial diagonal elements of u");
        for i in 0..mnmin {
          print!("{:.6} ",u[idx(i,i)]);
        }
        println!("\n initial diagonal elements of v");
        for i in 0..mnmin {
          print!("{:.6} ",v[idx(i,i)]);
        }
        print!("\n");
    }

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
        update_intermed_vars(&u, &v, &p, fsdx, fsdy, &mut cu, &mut cv, &mut z, &mut h);

        let mut c2 = tstart.elapsed().as_secs_f64();
        t100 = t100 + (c2 - c1);

        // apply periodic boundary conditions to intermediate variables
        apply_intermed_bcs(&mut cu, &mut cv, &mut z, &mut h);

        // time update to new variables
        let tdts8 = tdt / 8.0;
        let tdtsdx = tdt / dx;
        let tdtsdy = tdt / dy;

        c1 = tstart.elapsed().as_secs_f64();

        time_update_new_vars(&uold, &vold, &pold, &cu, &cv, &z, &h, tdts8, tdtsdx, tdtsdy, &mut unew, &mut vnew, &mut pnew);

        c2 = tstart.elapsed().as_secs_f64();
        t200 = t200 + (c2 - c1);
        
        // apply periodic boundary conitions to new variables
        apply_uvp_bcs(&mut unew, &mut vnew, &mut pnew);

        // update time
        time = time + dt;

        // update update old vars and solution
        if ncycle > 1 {

            c1 = tstart.elapsed().as_secs_f64();

            // smooth old vars using time filter
            smooth_update_old_vars(&u, &v, &p, &unew, &vnew, &pnew, &mut uold, &mut vold, &mut pold, alpha);

            // update u, v, and p to new solution
            mem::swap(&mut u, &mut unew);
            mem::swap(&mut v, &mut vnew);
            mem::swap(&mut p, &mut pnew);

            c2 = tstart.elapsed().as_secs_f64(); 
            t300 = t300 + (c2 - c1);
        } else {
            // update tdt for subsequent timesteps
            tdt = tdt + tdt;

            // no smoothing for first timestep
            // this might be redundant
            mem::swap(&mut uold, &mut u);
            mem::swap(&mut vold, &mut v);
            mem::swap(&mut pold, &mut p);

            // update u, v, and p to new solution
            // might be able to take out of if statement
            mem::swap(&mut u, &mut unew);
            mem::swap(&mut v, &mut vnew);
            mem::swap(&mut p, &mut pnew);
        }
    }

    // End time
    let ctime = tstart.elapsed().as_secs_f64();

    let ptime = time / 3600.;

    // print out final values
    if VERBOSE {
        let mnmin = M.min(N);
        println!(" diagonal elements of p");
        for i in 0..mnmin {
        print!("{:?} ",pnew[idx(i,i)]);
        }
        println!("\n diagonal elements of u");
        for i in 0..mnmin {
        print!("{:?} ",unew[idx(i,i)]);
        }
        println!("\n diagonal elements of v");
        for i in 0..mnmin {
        print!("{:?} ",vnew[idx(i,i)]);
        }
        print!("\n");
    }
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

    if CSV_OUT {
        #[cfg(feature = "box")]
        let filename = "rust_box_times.csv";

        #[cfg(feature = "vec")]
        let filename = "rust_vec_times.csv";

        #[cfg(feature = "ndarray")]
        let filename = "rust_ndarray_times.csv";

        #[cfg(feature = "mdarray")]
        let filename = "rust_mdarray_times.csv";

        let _ = write_int_float_to_csv(filename,M,ctime);
    }

    if SUCCINCT {
        #[cfg(feature = "box")]
        println!("Version: box");

        #[cfg(feature = "vec")]
        println!("Version: vec");

        #[cfg(feature = "ndarray")]
        println!("Version: ndarray");

        #[cfg(feature = "mdarray")]
        println!("Version: mdarray");

        println!("Grid Size: {:?}x{:?}, Number of iterations: {:?}, Total computer time: {:.2}", M, N, ITMAX, ctime);
    }

    // -----------------------------------------------------------------------
    // End
    // -----------------------------------------------------------------------

}
