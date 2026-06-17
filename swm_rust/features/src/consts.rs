// declare constants
pub const M: usize = 256;
pub const N: usize = 256;
pub const M_LEN: usize = M+1;
pub const N_LEN: usize = N+1;
pub const TOT_LEN: usize = (M_LEN)*(N_LEN);
pub const ITMAX: usize = 4000;
pub const VERBOSE: bool = false;  // print out initial and final values
pub const TIMING: bool = true;    // print out timings
pub const VAL_OUT: bool = false;  // save final solutions to txt files