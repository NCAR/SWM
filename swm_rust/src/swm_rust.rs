// Timing and environment stuff 
// use std::time::Instant;
// use std::env;

use std::f64::consts;
use std::mem;

// Utils is a helper module that contains some utility functions in src/utils.rs
mod utils;
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

    let el: f64 = N as f64 * dx;
    let pi = consts::PI;
    let tpi: f64 = pi + pi;
    let di: f64 = tpi / M as f64;
    let dj: f64 = tpi / N as f64;
    let pcf: f64 = pi * pi * a * a / (el * el);

    // indices
    let mut idx: usize;
    let mut idx01: usize;
    let mut idx10: usize;
    let mut idx11: usize;

    // -----------------------------------------------------------------------
    // Define Solution Arrays
    // -----------------------------------------------------------------------

    let mut u: Box<[f64]> = Box::from(vec![0.0; TOT_LEN]);
    let mut v: Box<[f64]> = Box::from(vec![0.0; TOT_LEN]);
    let mut p: Box<[f64]> = Box::from(vec![0.0; TOT_LEN]);
    let mut unew: Box<[f64]> = Box::from(vec![0.0; TOT_LEN]);
    let mut vnew: Box<[f64]> = Box::from(vec![0.0; TOT_LEN]);
    let mut pnew: Box<[f64]> = Box::from(vec![0.0; TOT_LEN]);
    let mut uold: Box<[f64]> = Box::from(vec![0.0; TOT_LEN]);
    let mut vold: Box<[f64]> = Box::from(vec![0.0; TOT_LEN]);
    let mut pold: Box<[f64]> = Box::from(vec![0.0; TOT_LEN]);
    let mut cu: Box<[f64]> = Box::from(vec![0.0; TOT_LEN]);
    let mut cv: Box<[f64]> = Box::from(vec![0.0; TOT_LEN]);
    let mut z: Box<[f64]> = Box::from(vec![0.0; TOT_LEN]);
    let mut h: Box<[f64]> = Box::from(vec![0.0; TOT_LEN]);
    let mut psi: Box<[f64]> = Box::from(vec![0.0; TOT_LEN]);

    // -----------------------------------------------------------------------
    // Initialize Data
    // -----------------------------------------------------------------------

    // initialize stream function psi and pressure p
    for i in 0..M_LEN {
        for j in 0..N_LEN {
            idx = ij_to_idx(i,j); // [i][j]
            psi[idx] = a * ( ( (i as f64) + 0.5 ) * di ).sin() * ( ( (j as f64) + 0.5 ) * dj ).sin();
            p[idx] = pcf * ( ( 2.0 * (i as f64) * di ).cos() + ( 2.0 * (j as f64) * dj ).cos() ) + 50000.;
        }
    }

    // initialize velocities u and v
    for i in 0..M {
        for j in 0..N {
            idx01 = ij_to_idx(i,j+1); // [i][j+1]
            idx10 = ij_to_idx(i+1,j); // [i+1][j]
            idx11 = ij_to_idx(i+1,j+1); // [i+1][j+1]
            u[idx10] = -(psi[idx11] - psi[idx10]) / dy;
            v[idx01] = (psi[idx11] - psi[idx01]) / dx;
        }
    }

    // periodic boundary conditions
    apply_uv_bcs(&mut u, &mut v);

    // initialize old arrays
    for i in 0..M_LEN {
        for j in 0..N_LEN {
            idx = ij_to_idx(i,j);
            uold[idx] = u[idx];
            vold[idx] = v[idx];
            pold[idx] = p[idx];
        }
    }

    // print out initial values
    // println!(" number of points in the x direction {:}", N); 
    // println!(" number of points in the y direction {:}", M); 
    // println!(" grid spacing in the x direction     {:}", dx); 
    // println!(" grid spacing in the y direction     {:}", dy); 
    // println!(" time step                           {:}", dt); 
    // println!(" time filter parameter               {:}", alpha);

    // let mnmin = M.min(N);
    // println!(" initial diagonal elements of p");
    // for i in 0..mnmin {
    //   print!("{:.6} ",p[ij_to_idx(i,i)]);
    // }
    // println!("\n initial diagonal elements of u");
    // for i in 0..mnmin {
    //   print!("{:.6} ",u[ij_to_idx(i,i)]);
    // }
    // println!("\n initial diagonal elements of v");
    // for i in 0..mnmin {
    //   print!("{:.6} ",v[ij_to_idx(i,i)]);
    // }
    // print!("\n");

    // print initial values to txt files
    // print_data_to_file("uinit_rust.txt", &u);
    // print_data_to_file("vinit_rust.txt", &v);
    // print_data_to_file("pinit_rust.txt", &p);

    // -----------------------------------------------------------------------
    // Time Marching Loop
    // -----------------------------------------------------------------------

    for ncycle in 1..=ITMAX { // fix
        // compute intermediate variables cu, cv, z, and h using u, v, and p
        update_intermed_vars(&u, &v, &p, fsdx, fsdy, &mut cu, &mut cv, &mut z, &mut h);

        // apply periodic boundary conditions to intermediate variables
        apply_intermed_bcs(&mut cu, &mut cv, &mut z, &mut h);

        // time update to new variables
        let tdts8 = tdt / 8.0;
        let tdtsdx = tdt / dx;
        let tdtsdy = tdt / dy;
        time_update_new_vars(&uold, &vold, &pold, &cu, &cv, &z, &h, tdts8, tdtsdx, tdtsdy, &mut unew, &mut vnew, &mut pnew);
        
        // apply periodic boundary conitions to new variables
        apply_uvp_bcs(&mut unew, &mut vnew, &mut pnew);

        // update update old vars and solution
        if ncycle > 1 {
            // smooth old vars using time filter
            smooth_update_old_vars(&u, &v, &p, &unew, &vnew, &pnew, &mut uold, &mut vold, &mut pold, alpha);

            // update u, v, and p to new solution
            mem::swap(&mut u, &mut unew);
            mem::swap(&mut v, &mut vnew);
            mem::swap(&mut p, &mut pnew);
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

    // -----------------------------------------------------------------------
    // End
    // -----------------------------------------------------------------------
    
    // print data to txt files
    // print_data_to_file("cu_rust.txt", &cu);
    // print_data_to_file("cv_rust.txt", &cv);
    // print_data_to_file("z_rust.txt", &z);
    // print_data_to_file("h_rust.txt", &h);
    // print_data_to_file("u_rust.txt", &u);
    // print_data_to_file("v_rust.txt", &v);
    // print_data_to_file("p_rust.txt", &p);


    // code for arg parsing and timing

    // let args: Vec<String> = env::args().collect();

    // let in_val: i32 = args[1].parse().expect("Input not an integer");

    // let start = Instant::now();
    // println!("Hello, world!");
    // let result = rand_util(in_val);
    // println!("Result: {}", result);

    // let elapsed_time = start.elapsed();
    // println!("Elapsed time: {:?}", elapsed_time.as_secs_f64());
}