# SWM - Rust implementation

This folder houses source code and project informatin for the Rust implementation of the SWM mini-app.
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

Probably the same as locally but need to double check. 

### Compile and run

#### Cargo (recommended)

To build only with cargo use : `cargo build`

To build (if files have changed) and run : `cargo run {optional arguments}`

To find errors without actually building : `cargo check`
Cargo also offers a more indepth explaination of errors with either `rustc --explain {CODE}` or `cargo --explain {CODE}` where `CODE` is the error code provided by the compiler (ex. E0308).

Run `cargo --help` for additional information on using cargo. 

#### C-like version

Compile with the Rust compiler:

`rustc main.rs`

and then run the executable with:

`./main {optional arguments}`

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

We will want to look for CUDA/GPU support crates eventually. 
