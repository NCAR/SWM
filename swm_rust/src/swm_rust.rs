// Timing and environment stuff 
use std::time::Instant;
use std::env;

// Utils is a helper module that contains some utility functions in src/utils.rs
mod utils;
use utils::*;

fn main() {

    let args: Vec<String> = env::args().collect();

    let in_val: i32 = args[1].parse().expect("Input not an integer");

    let start = Instant::now();
    println!("Hello, world!");
    let result = rand_util(in_val);
    println!("Result: {}", result);

    let elapsed_time = start.elapsed();
    println!("Elapsed time: {:?}", elapsed_time.as_secs_f64());
}