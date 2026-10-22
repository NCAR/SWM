use crate::consts::*;

// --- Crates ---
#[cfg(feature = "ndarray")]
use ndarray::Array2;

#[cfg(feature = "mdarray")]
use mdarray::DArray;

// --- Arr type alias ---
#[cfg(feature = "box")]
pub type Arr = Box<[f64]>;

#[cfg(feature = "vec")]
pub type Arr = Vec<f64>;

#[cfg(feature = "ndarray")]
pub type Arr = Array2<f64>;

#[cfg(feature = "mdarray")]
pub type Arr = DArray<f64, 2>;

// --- Index type alias ---
#[cfg(any(feature = "box", feature = "vec"))]
pub type Idx = usize;

#[cfg(any(feature = "ndarray", feature = "mdarray"))]
pub type Idx = [usize; 2];

// --- Index function ---
#[cfg(any(feature = "box", feature = "vec"))]
#[inline(always)]
pub fn idx(i: usize, j: usize) -> Idx {
    i * N_LEN + j
}

#[cfg(any(feature = "ndarray", feature = "mdarray"))]
#[inline(always)]
pub fn idx(i: usize, j: usize) -> Idx {
    [i, j]
}

// --- Initialize Arr ---
pub fn make_arr() -> Arr {
    #[cfg(feature = "box")]
    return Box::from(vec![0.0; TOT_LEN]);

    #[cfg(feature = "vec")]
    return vec![0.0; TOT_LEN];

    #[cfg(feature = "ndarray")]
    return Array2::zeros((M_LEN,N_LEN));

    #[cfg(feature = "mdarray")]
    return DArray::<f64, 2>::zeros([M_LEN,N_LEN]);
}