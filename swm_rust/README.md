# SWM - Rust implementation

This folder houses source code and project information for the Rust implementation of the SWM mini-app.
General information about Rust can be found [here](https://doc.rust-lang.org/book/ch01-00-getting-started.html). 

The Rust-lang book has an extensive programming guide [online](https://doc.rust-lang.org/book/ch03-00-common-programming-concepts.html).

## Install, compiling, running

### Install locally
Install on Linux/MacOS with 

`curl --proto '=https' --tlsv1.2 https://sh.rustup.rs -sSf | sh` 

and restart your shell. 
Check the installation with 

`rustc --version` 

or 

`cargo --version` 

Additional information and troubleshooting installs can be found [here](https://doc.rust-lang.org/book/ch01-01-installation.html).

### Install on Derecho

Same as local installation, no additional modules needed. May need to run

`module --force purge`

before compiling Rust code on Derecho.

### Compile and run

#### Cargo (recommended)

To build only with cargo use : `cargo build`

To build (if files have changed) and run : `cargo run {optional arguments}`

Default is to run debugging, flag for no debugging / best performance : `--release`

Flag to optimize for specific CPU architecture : `RUSTFLAGS="-C target-cpu=native"`

To find errors without actually building : `cargo check`
Cargo also offers a more indepth explaination of errors with either `rustc --explain {CODE}` or `cargo --explain {CODE}` where `CODE` is the error code provided by the compiler (ex. E0308).

Run `cargo --help` for additional information on using cargo. 

#### C-like version

Compile with the Rust compiler:

`rustc main.rs`

and then run the executable with:

`./main {optional arguments}`

### Features

Features version is currently housed in `features`. This version uses nested for loops to updates arrays in the chosen type. Use features flags `box`, `vec`, `ndarray`, and `mdarray` to change array type. Default is to run with Box. For example, to run with Vec : 

`cargo run --features vec --release`

### Zip ndarray

In `zip_ndarray`, ndarray array types are used and filled using the `zip` function in place of nested for loops. This has shown to be the fastest method in serial.

### Zip ndarray rayon

In `zip_ndarray_rayon`, the code from `zip_ndarray` is parellized with shared memory using the `rayon` crate. This is done by replacing `.for_each` with `.par_for_each` in the zip function calls. The number of threads can by specified at compile time with

`RAYON_NUM_THREADS={num}`

### Grid size

Grid size is a compile time argument. Default is 256x256. For example, to run on a 2048x2048 grid : 

`M=2048 N=2048 cargo run --release`

## Crates - packages for Rust

Rust utilizes packages, similar to how Python uses packages. These are known as crates.
An example crate is the `clap` crate for command line argument parsing. 

`clap` can be added to our project by using

`cargo add clap`. This automatically adds it to the `Cargo.toml` file as a dependency with its version number. 

To use functions from `clap`, we add something like 

`use clap::Parser;`

to the top of our `.rs` file and then use as we would any other Rust function. 

## Performance related notes

Cargo does support a test suite if needed, however I don't envision using that for performance testing work. 

`swm_rust.rs` currently contains an example for timing. This should be comparable to how walltime is computed in the C examples.

We will want to look for CUDA/GPU support crates eventually. See [NERSC docs](https://docs.nersc.gov/development/languages/rust/) for additional information on how to setup Rust to run efficiently on an HPC system. 
