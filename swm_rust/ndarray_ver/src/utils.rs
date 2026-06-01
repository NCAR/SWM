use std::fs::File;
use std::path::Path;
use std::io::prelude::*;
use ndarray::Array2;

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

// pub fn ij_to_idx(i: usize, j: usize) -> usize {
//     i*N_LEN + j
// }

pub fn apply_uv_bcs(u: &mut Array2<f64>, v: &mut Array2<f64>) {
    for j in 0..N {
        u[[0,j]] = u[[M,j]];
        v[[M,j+1]] = v[[0,j+1]];
    }

    for i in 0..M {
        u[[i+1,N]] = u[[i+1,0]];
        v[[i,0]] = v[[i,N]];
    }

    u[[0,N]] = u[[M,0]];
    v[[M,0]] = v[[0,N]];
}

pub fn update_intermed_vars(u: &Array2<f64>, v: &Array2<f64>, p: &Array2<f64>, fsdx: f64, fsdy: f64, cu: &mut Array2<f64>, cv: &mut Array2<f64>, z: &mut Array2<f64>, h: &mut Array2<f64>) {
    for i in 0..M {
        for j in 0..N {
            // let idx00 = ij_to_idx(i,j);
            // let idx01 = ij_to_idx(i,j+1);
            // let idx10 = ij_to_idx(i+1,j);
            // let idx11 = ij_to_idx(i+1,j+1);
            cu[[i+1,j]] = 0.5 * (p[[i+1,j]] + p[[i,j]]) * u[[i+1,j]];
            cv[[i,j+1]] = 0.5 * (p[[i,j+1]] + p[[i,j]]) * v[[i,j+1]];
            z[[i+1,j+1]] = (fsdx * (v[[i+1,j+1]] - v[[i,j+1]]) - fsdy * (u[[i+1,j+1]] - u[[i+1,j]])) / (p[[i,j]] + p[[i+1,j]] + p[[i+1,j+1]] + p[[i,j+1]]);
            h[[i,j]] = p[[i,j]] + 0.25 * (u[[i+1,j]] * u[[i+1,j]] + u[[i,j]] * u[[i,j]] + v[[i,j+1]] * v[[i,j+1]] + v[[i,j]] * v[[i,j]]);
        }
    }
}

pub fn apply_intermed_bcs(cu: &mut Array2<f64>, cv: &mut Array2<f64>, z: &mut Array2<f64>, h: &mut Array2<f64>) {
    for j in 0..N {
        cu[[0,j]] = cu[[M,j]];
        cv[[M,j+1]] = cv[[0,j+1]];
        z[[0,j+1]] = z[[M,j+1]];
        h[[M,j]] = h[[0,j]];
    }
    
    for i in 0..M {
        cu[[i+1,N]] = cu[[i+1,0]];
        cv[[i,0]] = cv[[i,N]];
        z[[i+1,0]] = z[[i+1,N]];
        h[[i,N]] = h[[i,0]];
    }

    // cu[N] = cu[M*N_LEN];
    // cv[M*N_LEN] = cv[N];
    // z[0] = z[M*N_LEN+N];
    // h[M*N_LEN+N] = h[0];
    cu[[0,N]] = cu[[M,0]];
    cv[[M,0]] = cv[[0,N]];
    z[[0,0]] = z[[M,N]];
    h[[M,N]] = h[[0,0]];
}

pub fn time_update_new_vars(uold: &Array2<f64>, vold: &Array2<f64>, pold: &Array2<f64>, cu: &Array2<f64>, cv: &Array2<f64>, z: &Array2<f64>, h: &Array2<f64>, tdts8: f64, tdtsdx: f64, tdtsdy: f64, unew: &mut Array2<f64>, vnew: &mut Array2<f64>, pnew: &mut Array2<f64>) {
    for i in 0..M {
        for j in 0..N {
            // let idx00 = ij_to_idx(i,j);
            // let idx01 = ij_to_idx(i,j+1);
            // let idx10 = ij_to_idx(i+1,j);
            // let idx11 = ij_to_idx(i+1,j+1);
            unew[[i+1,j]] = uold[[i+1,j]] + tdts8 * (z[[i+1,j+1]] + z[[i+1,j]]) * (cv[[i+1,j+1]] + cv[[i,j+1]] + cv[[i,j]] + cv[[i+1,j]]) - tdtsdx * (h[[i+1,j]] - h[[i,j]]);
            vnew[[i,j+1]] = vold[[i,j+1]] - tdts8 * (z[[i+1,j+1]] + z[[i,j+1]]) * (cu[[i+1,j+1]] + cu[[i,j+1]] + cu[[i,j]] + cu[[i+1,j]]) - tdtsdy * (h[[i,j+1]] - h[[i,j]]);
            pnew[[i,j]] = pold[[i,j]] - tdtsdx * (cu[[i+1,j]] - cu[[i,j]]) - tdtsdy * (cv[[i,j+1]] - cv[[i,j]]);
        }
    }
}

pub fn apply_uvp_bcs(u: &mut Array2<f64>, v: &mut Array2<f64>, p: &mut Array2<f64>) {
    for j in 0..N {
        u[[0,j]] = u[[M,j]];
        v[[M,j+1]] = v[[0,j+1]];
        p[[M,j]] = p[[0,j]];
    }

    for i in 0..M {
        u[[i+1,N]] = u[[i+1,0]];
        v[[i,0]] = v[[i,N]];
        p[[i,N]] = p[[i,0]];
    }

    u[[0,N]] = u[[M,0]];
    v[[M,0]] = v[[0,N]];
    p[[M,N]] = p[[0,0]];
}

pub fn smooth_update_old_vars(u: &Array2<f64>, v: &Array2<f64>, p: &Array2<f64>, unew: &Array2<f64>, vnew: &Array2<f64>, pnew: &Array2<f64>, uold: &mut Array2<f64>, vold: &mut Array2<f64>, pold: &mut Array2<f64>, alpha: f64) {
    for i in 0..M_LEN {
        for j in 0..N_LEN {
            // let idx = ij_to_idx(i,j);
            uold[[i,j]] = u[[i,j]] + alpha * (unew[[i,j]] - 2. * u[[i,j]] + uold[[i,j]]);
            vold[[i,j]] = v[[i,j]] + alpha * (vnew[[i,j]] - 2. * v[[i,j]] + vold[[i,j]]);
            pold[[i,j]] = p[[i,j]] + alpha * (pnew[[i,j]] - 2. * p[[i,j]] + pold[[i,j]]);
        }
    }
}

pub fn print_data_to_file(pathname: &str, data: &Array2<f64>) {
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
            s += &format!("{:.6} ", data[[i,j]]);
        }
        s += "\n";
    }

    // Write string to `file`, returns `io::Result<()>`
    match file.write_all(s.as_bytes()) {
        Err(why) => panic!("couldn't write to {}: {}", display, why),
        Ok(_) => println!("successfully wrote to {}", display),
    }
}