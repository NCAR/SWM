use std::fs::File;
use std::path::Path;
use std::io::prelude::*;
use std::f64::consts;

use crate::consts::*;
use crate::types::{Arr,idx,make_arr};

pub fn init_conds(u: &mut Arr, v: &mut Arr, p: &mut Arr, dx: f64, dy: f64, a: f64) { 
    // init psi
    let mut psi: Arr = make_arr();

    // set params
    let el: f64 = N as f64 * dx;
    let pi = consts::PI;
    let tpi: f64 = pi + pi;
    let di: f64 = tpi / M as f64;
    let dj: f64 = tpi / N as f64;
    let pcf: f64 = pi * pi * a * a / (el * el);
    
    // initialize stream function psi and pressure p
    for i in 0..M_LEN {
        for j in 0..N_LEN {
            psi[idx(i,j)] = a * ( ( (i as f64) + 0.5 ) * di ).sin() * ( ( (j as f64) + 0.5 ) * dj ).sin();
            p[idx(i,j)] = pcf * ( ( 2.0 * (i as f64) * di ).cos() + ( 2.0 * (j as f64) * dj ).cos() ) + 50000.;
        }
    }

    // initialize velocities u and v
    for i in 0..M {
        for j in 0..N {
            u[idx(i+1,j)] = -(psi[idx(i+1,j+1)] - psi[idx(i+1,j)]) / dy;
            v[idx(i,j+1)] = (psi[idx(i+1,j+1)] - psi[idx(i,j+1)]) / dx;
        }
    }
}

pub fn apply_uv_bcs(u: &mut Arr, v: &mut Arr) {
    for j in 0..N {
        u[idx(0,j)] = u[idx(M,j)];
        v[idx(M,j+1)] = v[idx(0,j+1)];
    }

    for i in 0..M {
        u[idx(i+1,N)] = u[idx(i+1,0)];
        v[idx(i,0)] = v[idx(i,N)];
    }

    u[idx(0,N)] = u[idx(M,0)];
    v[idx(M,0)] = v[idx(0,N)];
}

pub fn update_intermed_vars(u: &Arr, v: &Arr, p: &Arr, fsdx: f64, fsdy: f64, cu: &mut Arr, cv: &mut Arr, z: &mut Arr, h: &mut Arr) {
    for i in 0..M {
        for j in 0..N {
            cu[idx(i+1,j)] = 0.5 * (p[idx(i+1,j)] + p[idx(i,j)]) * u[idx(i+1,j)];
            cv[idx(i,j+1)] = 0.5 * (p[idx(i,j+1)] + p[idx(i,j)]) * v[idx(i,j+1)];
            z[idx(i+1,j+1)] = (fsdx * (v[idx(i+1,j+1)] - v[idx(i,j+1)]) - fsdy * (u[idx(i+1,j+1)] - u[idx(i+1,j)])) / (p[idx(i,j)] + p[idx(i+1,j)] + p[idx(i+1,j+1)] + p[idx(i,j+1)]);
            h[idx(i,j)] = p[idx(i,j)] + 0.25 * (u[idx(i+1,j)] * u[idx(i+1,j)] + u[idx(i,j)] * u[idx(i,j)] + v[idx(i,j+1)] * v[idx(i,j+1)] + v[idx(i,j)] * v[idx(i,j)]);
        }
    }
}

pub fn apply_intermed_bcs(cu: &mut Arr, cv: &mut Arr, z: &mut Arr, h: &mut Arr) {
    for j in 0..N {
        cu[idx(0,j)] = cu[idx(M,j)];
        cv[idx(M,j+1)] = cv[idx(0,j+1)];
        z[idx(0,j+1)] = z[idx(M,j+1)];
        h[idx(M,j)] = h[idx(0,j)];
    }
    
    for i in 0..M {
        cu[idx(i+1,N)] = cu[idx(i+1,0)];
        cv[idx(i,0)] = cv[idx(i,N)];
        z[idx(i+1,0)] = z[idx(i+1,N)];
        h[idx(i,N)] = h[idx(i,0)];
    }

    cu[idx(0,N)] = cu[idx(M,0)];
    cv[idx(M,0)] = cv[idx(0,N)];
    z[idx(0,0)] = z[idx(M,N)];
    h[idx(M,N)] = h[idx(0,0)];
}

pub fn time_update_new_vars(uold: &Arr, vold: &Arr, pold: &Arr, cu: &Arr, cv: &Arr, z: &Arr, h: &Arr, tdts8: f64, tdtsdx: f64, tdtsdy: f64, unew: &mut Arr, vnew: &mut Arr, pnew: &mut Arr) {
    for i in 0..M {
        for j in 0..N {
            unew[idx(i+1,j)] = uold[idx(i+1,j)] + tdts8 * (z[idx(i+1,j+1)] + z[idx(i+1,j)]) * (cv[idx(i+1,j+1)] + cv[idx(i,j+1)] + cv[idx(i,j)] + cv[idx(i+1,j)]) - tdtsdx * (h[idx(i+1,j)] - h[idx(i,j)]);
            vnew[idx(i,j+1)] = vold[idx(i,j+1)] - tdts8 * (z[idx(i+1,j+1)] + z[idx(i,j+1)]) * (cu[idx(i+1,j+1)] + cu[idx(i,j+1)] + cu[idx(i,j)] + cu[idx(i+1,j)]) - tdtsdy * (h[idx(i,j+1)] - h[idx(i,j)]);
            pnew[idx(i,j)] = pold[idx(i,j)] - tdtsdx * (cu[idx(i+1,j)] - cu[idx(i,j)]) - tdtsdy * (cv[idx(i,j+1)] - cv[idx(i,j)]);
        }
    }
}

pub fn apply_uvp_bcs(u: &mut Arr, v: &mut Arr, p: &mut Arr) {
    for j in 0..N {
        u[idx(0,j)] = u[idx(M,j)];
        v[idx(M,j+1)] = v[idx(0,j+1)];
        p[idx(M,j)] = p[idx(0,j)];
    }

    for i in 0..M {
        u[idx(i+1,N)] = u[idx(i+1,0)];
        v[idx(i,0)] = v[idx(i,N)];
        p[idx(i,N)] = p[idx(i,0)];
    }

    u[idx(0,N)] = u[idx(M,0)];
    v[idx(M,0)] = v[idx(0,N)];
    p[idx(M,N)] = p[idx(0,0)];
}

pub fn smooth_update_old_vars(u: &Arr, v: &Arr, p: &Arr, unew: &Arr, vnew: &Arr, pnew: &Arr, uold: &mut Arr, vold: &mut Arr, pold: &mut Arr, alpha: f64) {
    for i in 0..M_LEN {
        for j in 0..N_LEN {
            uold[idx(i,j)] = u[idx(i,j)] + alpha * (unew[idx(i,j)] - 2. * u[idx(i,j)] + uold[idx(i,j)]);
            vold[idx(i,j)] = v[idx(i,j)] + alpha * (vnew[idx(i,j)] - 2. * v[idx(i,j)] + vold[idx(i,j)]);
            pold[idx(i,j)] = p[idx(i,j)] + alpha * (pnew[idx(i,j)] - 2. * p[idx(i,j)] + pold[idx(i,j)]);
        }
    }
}

pub fn print_data_to_file(pathname: &str, data: &Arr) {
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
            s += &format!("{:.6} ", data[idx(i,j)]);
        }
        s += "\n";
    }

    // Write string to `file`, returns `io::Result<()>`
    match file.write_all(s.as_bytes()) {
        Err(why) => panic!("couldn't write to {}: {}", display, why),
        Ok(_) => println!("successfully wrote to {}", display),
    }
}
