use std::fs::File;
use std::path::Path;
use std::io::prelude::*;
use std::f64::consts;
use ndarray::{Array2,Zip,s};

use csv::WriterBuilder;
use std::error::Error;
use std::fs::OpenOptions;

// declare constants
pub const M: usize = const_env::env_lit!("M", 256);
pub const N: usize = const_env::env_lit!("N", 256);
pub const M_LEN: usize = M+1;
pub const N_LEN: usize = N+1;
pub const ITMAX: usize = 4000;
pub const VERBOSE: bool = false;  // print out initial and final values
pub const TIMING: bool = false;   // print out timings
pub const VAL_OUT: bool = false;  // save final solutions to txt files
pub const SUCCINCT: bool = true;  // print out grid size, itmax, and final time
pub const CSV_OUT: bool = true;  // save time to csv file

pub fn init_conds(u: &mut Array2<f64>, v: &mut Array2<f64>, p: &mut Array2<f64>, dx: f64, dy: f64, a: f64) { 
    // init psi
    let mut psi = Array2::<f64>::zeros((M_LEN,N_LEN));

    // set params
    let el: f64 = N as f64 * dx;
    let pi = consts::PI;
    let tpi: f64 = pi + pi;
    let di: f64 = tpi / M as f64;
    let dj: f64 = tpi / N as f64;
    let pcf: f64 = pi * pi * a * a / (el * el);

    // initialize stream function psi and pressure p
    // psi[[i,j]] = a * ( ( (i as f64) + 0.5 ) * di ).sin() * ( ( (j as f64) + 0.5 ) * dj ).sin()
    Zip::indexed(psi.view_mut())
        .par_for_each(|(i, j), val| {
            *val = a * ( ( (i as f64) + 0.5 ) * di ).sin() * ( ( (j as f64) + 0.5 ) * dj ).sin();
        });
    // p[[i,j]] = pcf * ( ( 2.0 * (i as f64) * di ).cos() + ( 2.0 * (j as f64) * dj ).cos() ) + 50000.
    Zip::indexed(p.view_mut())
        .par_for_each(|(i, j), val| {
            *val = pcf * ( ( 2.0 * (i as f64) * di ).cos() + ( 2.0 * (j as f64) * dj ).cos() ) + 50000.;
        });

    // initialize velocities u and v
    // u[[i+1,j]] = -(psi[[i+1,j+1]] - psi[[i+1,j]]) / dy;
    Zip::from(u.slice_mut(s![1.., ..N]))
        .and(psi.slice(s![1.., 1..]))
        .and(psi.slice(s![1.., ..N]))
        .par_for_each(|u_val, &psi_up, &psi_down| {
            *u_val = -(psi_up - psi_down) / dy;
        });
    // v[[i,j+1]] = (psi[[i+1,j+1]] - psi[[i,j+1]]) / dx;
    Zip::from(v.slice_mut(s![..M, 1..]))
        .and(psi.slice(s![1.., 1..]))
        .and(psi.slice(s![..M, 1..]))
        .par_for_each(|v_val, &psi_right, &psi_left| {
            *v_val = (psi_right - psi_left) / dx;
        });
}

pub fn apply_uv_bcs(u: &mut Array2<f64>, v: &mut Array2<f64>) {
    // u[[0,j]] = u[[M,j]]  for j in 0..N
    let (mut row0, rowm) = u.multi_slice_mut((s![0, ..N], s![M, ..N]));
    row0.assign(&rowm);

    // v[[M,j+1]] = v[[0,j+1]]  for j in 0..N
    let (mut rowm, row0) = v.multi_slice_mut((s![M, 1..], s![0, 1..]));
    rowm.assign(&row0);

    // u[[i+1,N]] = u[[i+1,0]]  for i in 0..M
    let (mut coln, col0) = u.multi_slice_mut((s![1.., N], s![1.., 0]));
    coln.assign(&col0);

    // v[[i,0]] = v[[i,N]]  for i in 0..M
    let (mut col0, coln) = v.multi_slice_mut((s![..M, 0], s![..M, N]));
    col0.assign(&coln);

    u[[0,N]] = u[[M,0]];
    v[[M,0]] = v[[0,N]];
}

pub fn update_intermed_vars(u: &Array2<f64>, v: &Array2<f64>, p: &Array2<f64>, fsdx: f64, fsdy: f64, cu: &mut Array2<f64>, cv: &mut Array2<f64>, z: &mut Array2<f64>, h: &mut Array2<f64>, term1: &mut Array2<f64>, term2: &mut Array2<f64>) {
    // cu[[i+1,j]] = 0.5 * (p[[i+1,j]] + p[[i,j]]) * u[[i+1,j]];
    Zip::from(cu.slice_mut(s![1.., ..N]))
        .and(p.slice(s![1.., ..N]))
        .and(p.slice(s![..M, ..N]))
        .and(u.slice(s![1.., ..N]))
        .par_for_each(|cu10, &p10, &p00, &u10| {
            *cu10 = 0.5 * (p10 + p00) * u10;
        });
    
    // cv[[i,j+1]] = 0.5 * (p[[i,j+1]] + p[[i,j]]) * v[[i,j+1]];
    Zip::from(cv.slice_mut(s![..M, 1..]))
        .and(p.slice(s![..M, 1..]))
        .and(p.slice(s![..M, ..N]))
        .and(v.slice(s![..M, 1..]))
        .par_for_each(|cv01, &p01, &p00, &v01| {
            *cv01 = 0.5 * (p01 + p00) * v01;
        });
    
    // z[[i+1,j+1]] = (fsdx * (v[[i+1,j+1]] - v[[i,j+1]]) - fsdy * (u[[i+1,j+1]] - u[[i+1,j]])) / (p[[i,j]] + p[[i+1,j]] + p[[i+1,j+1]] + p[[i,j+1]]);
    Zip::from(term1.view_mut())
        .and(v.slice(s![1.., 1..]))
        .and(v.slice(s![..M, 1..]))
        .and(u.slice(s![1.., 1..]))
        .and(u.slice(s![1.., ..N]))
        .par_for_each(|val, &v11, &v01, &u11, &u10| {
            *val = fsdx * (v11 - v01) - fsdy * (u11 - u10);
        });
    four_node_sum(p, term2);
    Zip::from(z.slice_mut(s![1.., 1..]))
        .and(term1.view())
        .and(term2.view())
        .par_for_each(|z11, &val1, &val2| {
            *z11 = val1 / val2;
        });
    
    // h[[i,j]] = p[[i,j]] + 0.25 * (u[[i+1,j]] * u[[i+1,j]] + u[[i,j]] * u[[i,j]] + v[[i,j+1]] * v[[i,j+1]] + v[[i,j]] * v[[i,j]]);
    Zip::from(h.slice_mut(s![..M, ..N]))
        .and(u.slice(s![1.., ..N]))
        .and(u.slice(s![..M, ..N]))
        .and(v.slice(s![..M, 1..]))
        .and(v.slice(s![..M, ..N]))
        .and(p.slice(s![..M, ..N]))
        .par_for_each(|h00, &u10, &u00, &v01, &v00, &p00| {
            *h00 = &p00 + 0.25 * (&u10 * &u10 + &u00 * &u00 + &v01 * &v01 + &v00 * &v00);
        });
}

pub fn apply_intermed_bcs(cu: &mut Array2<f64>, cv: &mut Array2<f64>, z: &mut Array2<f64>, h: &mut Array2<f64>) {
    // cu[[0,j]] = cu[[M,j]]  for j in 0..N
    let (mut row0, rowm) = cu.multi_slice_mut((s![0, ..N], s![M, ..N]));
    row0.assign(&rowm);

    // cv[[M,j+1]] = cv[[0,j+1]]  for j in 0..N
    let (mut rowm, row0) = cv.multi_slice_mut((s![M, 1..], s![0, 1..]));
    rowm.assign(&row0);

    // z[[0,j+1]] = z[[M,j+1]]  for j in 0..N
    let (mut row0, rowm) = z.multi_slice_mut((s![0, 1..], s![M, 1..]));
    row0.assign(&rowm);

    // h[[M,j]] = h[[0,j]]  for j in 0..N
    let (mut rowm, row0) = h.multi_slice_mut((s![M, ..N], s![0, ..N]));
    rowm.assign(&row0);

    // cu[[i+1,N]] = cu[[i+1,0]]  for i in 0..M
    let (mut coln, col0) = cu.multi_slice_mut((s![1.., N], s![1.., 0]));
    coln.assign(&col0);

    // cv[[i,0]] = cv[[i,N]]  for i in 0..M
    let (mut col0, coln) = cv.multi_slice_mut((s![..M, 0], s![..M, N]));
    col0.assign(&coln);

    // z[[i+1,0]] = z[[i+1,N]]  for i in 0..M
    let (mut coln, col0) = z.multi_slice_mut((s![1.., 0], s![1.., N]));
    coln.assign(&col0);

    // h[[i,N]] = h[[i,0]]  for i in 0..M
    let (mut col0, coln) = h.multi_slice_mut((s![..M, N], s![..M, 0]));
    col0.assign(&coln);

    cu[[0,N]] = cu[[M,0]];
    cv[[M,0]] = cv[[0,N]];
    z[[0,0]] = z[[M,N]];
    h[[M,N]] = h[[0,0]];
}

pub fn time_update_new_vars(uold: &Array2<f64>, vold: &Array2<f64>, pold: &Array2<f64>, cu: &Array2<f64>, cv: &Array2<f64>, z: &Array2<f64>, h: &Array2<f64>, tdts8: f64, tdtsdx: f64, tdtsdy: f64, unew: &mut Array2<f64>, vnew: &mut Array2<f64>, pnew: &mut Array2<f64>, term1: &mut Array2<f64>, term2: &mut Array2<f64>) {
    // unew[[i+1,j]] = uold[[i+1,j]] + tdts8 * (z[[i+1,j+1]] + z[[i+1,j]]) * (cv[[i+1,j+1]] + cv[[i,j+1]] + cv[[i,j]] + cv[[i+1,j]]) - tdtsdx * (h[[i+1,j]] - h[[i,j]]);
    four_node_sum(cv, term1);
    Zip::from(term2.view_mut())
        .and(z.slice(s![1.., 1..]))
        .and(z.slice(s![1.., ..N]))
        .and(term1.view())
        .par_for_each(|val2, &z11, &z10, &val1| {
            *val2 = tdts8 * (z11 + z10) * val1;
        });
    Zip::from(unew.slice_mut(s![1.., ..N]))
        .and(term2.view())
        .and(uold.slice(s![1.., ..N]))
        .and(h.slice(s![1.., ..N]))
        .and(h.slice(s![..M, ..N]))
        .par_for_each(|unew10, &val2, &uold10, &h10, &h00| {
            *unew10 = uold10 + val2 - tdtsdx * (h10 - h00);
        });

    // vnew[[i,j+1]] = vold[[i,j+1]] - tdts8 * (z[[i+1,j+1]] + z[[i,j+1]]) * (cu[[i+1,j+1]] + cu[[i,j+1]] + cu[[i,j]] + cu[[i+1,j]]) - tdtsdy * (h[[i,j+1]] - h[[i,j]]);
    four_node_sum(cu, term1);
    Zip::from(term2.view_mut())
        .and(z.slice(s![1.., 1..]))
        .and(z.slice(s![..M, 1..]))
        .and(term1.view())
        .par_for_each(|val2, &z11, &z01, &val1| {
            *val2 = - tdts8 * (z11 + z01) * val1;
        });
    Zip::from(vnew.slice_mut(s![..M, 1..]))
        .and(term2.view())
        .and(vold.slice(s![0..M, 1..]))
        .and(h.slice(s![..M, 1..]))
        .and(h.slice(s![..M, ..N]))
        .par_for_each(|vnew01, &val2, &vold01, &h01, &h00| {
            *vnew01 = vold01 + val2 - tdtsdy * (h01 - h00);
        });

    // pnew[[i,j]] = pold[[i,j]] - tdtsdx * (cu[[i+1,j]] - cu[[i,j]]) - tdtsdy * (cv[[i,j+1]] - cv[[i,j]]);
    Zip::from(pnew.slice_mut(s![..M, ..N]))
        .and(pold.slice(s![..M, ..N]))
        .and(cu.slice(s![1.., ..N]))
        .and(cu.slice(s![..M, ..N]))
        .and(cv.slice(s![..M, 1..]))
        .and(cv.slice(s![..M, ..N]))
        .par_for_each(|pnew00, &pold00, &cu10, &cu00, &cv01, &cv00| {
            *pnew00 = pold00 - tdtsdx * (cu10 - cu00) - tdtsdy * (cv01 - cv00);
        });
}

pub fn apply_uvp_bcs(u: &mut Array2<f64>, v: &mut Array2<f64>, p: &mut Array2<f64>) {
    // u[[0,j]] = u[[M,j]]  for j in 0..N
    let (mut row0, rowm) = u.multi_slice_mut((s![0, ..N], s![M, ..N]));
    row0.assign(&rowm);

    // v[[M,j+1]] = v[[0,j+1]]  for j in 0..N
    let (mut rowm, row0) = v.multi_slice_mut((s![M, 1..], s![0, 1..]));
    rowm.assign(&row0);

    // p[[M,j]] = p[[0,j]]  for j in 0..N
    let (mut row0, rowm) = p.multi_slice_mut((s![M, ..N], s![0, ..N]));
    row0.assign(&rowm);

    // u[[i+1,N]] = u[[i+1,0]]  for i in 0..M
    let (mut coln, col0) = u.multi_slice_mut((s![1.., N], s![1.., 0]));
    coln.assign(&col0);

    // v[[i,0]] = v[[i,N]]  for i in 0..M
    let (mut col0, coln) = v.multi_slice_mut((s![..M, 0], s![..M, N]));
    col0.assign(&coln);

    // p[[i,N]] = p[[i,0]]  for i in 0..M
    let (mut coln, col0) = p.multi_slice_mut((s![..M, N], s![..M, 0]));
    coln.assign(&col0);

    u[[0,N]] = u[[M,0]];
    v[[M,0]] = v[[0,N]];
    p[[M,N]] = p[[0,0]];
}

pub fn smooth_update_old_vars(u: &Array2<f64>, v: &Array2<f64>, p: &Array2<f64>, unew: &Array2<f64>, vnew: &Array2<f64>, pnew: &Array2<f64>, uold: &mut Array2<f64>, vold: &mut Array2<f64>, pold: &mut Array2<f64>, alpha: f64) {
    // uold[[i,j]] = u[[i,j]] + alpha * (unew[[i,j]] - 2. * u[[i,j]] + uold[[i,j]]);
    Zip::from(uold.view_mut())
        .and(u)
        .and(unew)
        .par_for_each(|uold00, &u00, &unew00| {
            *uold00 = u00 + alpha * (unew00 - 2. * u00 + *uold00);
        });

    // vold[[i,j]] = v[[i,j]] + alpha * (vnew[[i,j]] - 2. * v[[i,j]] + vold[[i,j]]);
    Zip::from(vold.view_mut())
        .and(v)
        .and(vnew)
        .par_for_each(|vold00, &v00, &vnew00| {
            *vold00 = v00 + alpha * (vnew00 - 2. * v00 + *vold00);
        });

    // pold[[i,j]] = p[[i,j]] + alpha * (pnew[[i,j]] - 2. * p[[i,j]] + pold[[i,j]]);
    Zip::from(pold.view_mut())
        .and(p)
        .and(pnew)
        .par_for_each(|pold00, &p00, &pnew00| {
            *pold00 = p00 + alpha * (pnew00 - 2. * p00 + *pold00);
        });
}

pub fn four_node_sum(arr: &Array2<f64>, sum: &mut Array2<f64>) {
    // sum[[i,j]] = arr[[i+1,j+1]] + arr[[i,j+1]] + arr[[i,j]] + arr[[i+1,j]];
    Zip::from(sum.view_mut())
        .and(arr.slice(s![1.., 1..]))
        .and(arr.slice(s![..M, 1..]))
        .and(arr.slice(s![..M, ..N]))
        .and(arr.slice(s![1.., ..N]))
        .par_for_each(|val, &arr11, &arr01, &arr00, &arr10| {
            *val = arr11 + arr01 + arr00 + arr10;
        });
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

pub fn write_int_float_to_csv(
    filename: &str,
    int_val: usize,
    float_val: f64,
) -> Result<(), Box<dyn Error>> {
    let file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(filename)?;

    let mut wtr = WriterBuilder::new()
        .has_headers(false)
        .from_writer(file);

    wtr.write_record(&[int_val.to_string(), float_val.to_string()])?;
    wtr.flush()?;

    Ok(())
}