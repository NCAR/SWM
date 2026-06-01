use std::fs::File;
use std::path::Path;
use std::io::prelude::*;

// declare constants
pub const M: usize = 256;
pub const N: usize = 256;
pub const M_LEN: usize = M+1;
pub const N_LEN: usize = N+1;
pub const TOT_LEN: usize = (M_LEN)*(N_LEN);
pub const ITMAX: usize = 1000;

// Below is a dummy variable to test stuff. 
// pub fn rand_util(in_val: i32) -> f64 {
//     let mut x = 0.0;
//     for i in 0..in_val {
//         x += (i as f64).sin() * (i as f64).cos();
//     }
//     // Default we return the last line that has no semicolon, but we can also use the return keyword if needed.
//     x
// }

pub fn ij_to_idx(i: usize, j: usize) -> usize {
    i*N_LEN + j
}

pub fn apply_uv_bcs(u: &mut Box<[f64]>, v: &mut Box<[f64]>) {
    for j in 0..N {
        u[ij_to_idx(0,j)] = u[ij_to_idx(M,j)];
        v[ij_to_idx(M,j+1)] = v[ij_to_idx(0,j+1)];
    }

    for i in 0..M {
        u[ij_to_idx(i+1,N)] = u[ij_to_idx(i+1,0)];
        v[ij_to_idx(i,0)] = v[ij_to_idx(i,N)];
    }

    u[ij_to_idx(0,N)] = u[ij_to_idx(M,0)];
    v[ij_to_idx(M,0)] = v[ij_to_idx(0,N)];
}

pub fn update_intermed_vars(u: &Box<[f64]>, v: &Box<[f64]>, p: &Box<[f64]>, fsdx: f64, fsdy: f64, cu: &mut Box<[f64]>, cv: &mut Box<[f64]>, z: &mut Box<[f64]>, h: &mut Box<[f64]>) {
    for i in 0..M {
        for j in 0..N {
            let idx00 = ij_to_idx(i,j);
            let idx01 = ij_to_idx(i,j+1);
            let idx10 = ij_to_idx(i+1,j);
            let idx11 = ij_to_idx(i+1,j+1);
            cu[idx10] = 0.5 * (p[idx10] + p[idx00]) * u[idx10];
            cv[idx01] = 0.5 * (p[idx01] + p[idx00]) * v[idx01];
            z[idx11] = (fsdx * (v[idx11] - v[idx01]) - fsdy * (u[idx11] - u[idx10])) / (p[idx00] + p[idx10] + p[idx11] + p[idx01]);
            h[idx00] = p[idx00] + 0.25 * (u[idx10] * u[idx10] + u[idx00] * u[idx00] + v[idx01] * v[idx01] + v[idx00] * v[idx00]);
        }
    }
}

pub fn apply_intermed_bcs(cu: &mut Box<[f64]>, cv: &mut Box<[f64]>, z: &mut Box<[f64]>, h: &mut Box<[f64]>) {
    for j in 0..N {
        cu[ij_to_idx(0,j)] = cu[ij_to_idx(M,j)];
        cv[ij_to_idx(M,j+1)] = cv[ij_to_idx(0,j+1)];
        z[ij_to_idx(0,j+1)] = z[ij_to_idx(M,j+1)];
        h[ij_to_idx(M,j)] = h[ij_to_idx(0,j)];
    }
    
    for i in 0..M {
        cu[ij_to_idx(i+1,N)] = cu[ij_to_idx(i+1,0)];
        cv[ij_to_idx(i,0)] = cv[ij_to_idx(i,N)];
        z[ij_to_idx(i+1,0)] = z[ij_to_idx(i+1,N)];
        h[ij_to_idx(i,N)] = h[ij_to_idx(i,0)];
    }

    cu[N] = cu[M*N_LEN];
    cv[M*N_LEN] = cv[N];
    z[0] = z[M*N_LEN+N];
    h[M*N_LEN+N] = h[0];
}

pub fn time_update_new_vars(uold: &Box<[f64]>, vold: &Box<[f64]>, pold: &Box<[f64]>, cu: &Box<[f64]>, cv: &Box<[f64]>, z: &Box<[f64]>, h: &Box<[f64]>, tdts8: f64, tdtsdx: f64, tdtsdy: f64, unew: &mut Box<[f64]>, vnew: &mut Box<[f64]>, pnew: &mut Box<[f64]>) {
    for i in 0..M {
        for j in 0..N {
            let idx00 = ij_to_idx(i,j);
            let idx01 = ij_to_idx(i,j+1);
            let idx10 = ij_to_idx(i+1,j);
            let idx11 = ij_to_idx(i+1,j+1);
            unew[idx10] = uold[idx10] + tdts8 * (z[idx11] + z[idx10]) * (cv[idx11] + cv[idx01] + cv[idx00] + cv[idx10]) - tdtsdx * (h[idx10] - h[idx00]);
            vnew[idx01] = vold[idx01] - tdts8 * (z[idx11] + z[idx01]) * (cu[idx11] + cu[idx01] + cu[idx00] + cu[idx10]) - tdtsdy * (h[idx01] - h[idx00]);
            pnew[idx00] = pold[idx00] - tdtsdx * (cu[idx10] - cu[idx00]) - tdtsdy * (cv[idx01] - cv[idx00]);
        }
    }
}

pub fn apply_uvp_bcs(u: &mut Box<[f64]>, v: &mut Box<[f64]>, p: &mut Box<[f64]>) {
    for j in 0..N {
        u[ij_to_idx(0,j)] = u[ij_to_idx(M,j)];
        v[ij_to_idx(M,j+1)] = v[ij_to_idx(0,j+1)];
        p[ij_to_idx(M,j)] = p[ij_to_idx(0,j)];
    }

    for i in 0..M {
        u[ij_to_idx(i+1,N)] = u[ij_to_idx(i+1,0)];
        v[ij_to_idx(i,0)] = v[ij_to_idx(i,N)];
        p[ij_to_idx(i,N)] = p[ij_to_idx(i,0)];
    }

    u[ij_to_idx(0,N)] = u[ij_to_idx(M,0)];
    v[ij_to_idx(M,0)] = v[ij_to_idx(0,N)];
    p[ij_to_idx(M,N)] = p[ij_to_idx(0,0)];
}

pub fn smooth_update_old_vars(u: &Box<[f64]>, v: &Box<[f64]>, p: &Box<[f64]>, unew: &Box<[f64]>, vnew: &Box<[f64]>, pnew: &Box<[f64]>, uold: &mut Box<[f64]>, vold: &mut Box<[f64]>, pold: &mut Box<[f64]>, alpha: f64) {
    for i in 0..M_LEN {
        for j in 0..N_LEN {
            let idx = ij_to_idx(i,j);
            uold[idx] = u[idx] + alpha * (unew[idx] - 2. * u[idx] + uold[idx]);
            vold[idx] = v[idx] + alpha * (vnew[idx] - 2. * v[idx] + vold[idx]);
            pold[idx] = p[idx] + alpha * (pnew[idx] - 2. * p[idx] + pold[idx]);
        }
    }
}

pub fn print_data_to_file(pathname: &str, data: &Box<[f64]>) {
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
            s += &format!("{:.6} ", data[ij_to_idx(i,j)]);
        }
        s += "\n";
    }

    // Write string to `file`, returns `io::Result<()>`
    match file.write_all(s.as_bytes()) {
        Err(why) => panic!("couldn't write to {}: {}", display, why),
        Ok(_) => println!("successfully wrote to {}", display),
    }
}