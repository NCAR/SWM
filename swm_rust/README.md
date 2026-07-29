# SWM - Rust implementation

This folder houses source code and project information for the Rust implementation of the SWM mini-app.
General information about Rust can be found [here](https://doc.rust-lang.org/book/ch01-00-getting-started.html). 

The [Rust-lang book](https://doc.rust-lang.org/book/) has an extensive programming guide.

The [Rust Cookbook](https://rust-lang-nursery.github.io/rust-cookbook/) is another helpful piece of documentation for developing Rust code.

[Rust by Example](https://doc.rust-lang.org/rust-by-example/) contains example code snippets for popular Rust concepts.

The [Rust Performance Book](https://nnethercote.github.io/perf-book/) has useful tips for maximizing performance in Rust.

See [NERSC docs](https://docs.nersc.gov/development/languages/rust/) for additional information on how to setup Rust to run efficiently on an HPC system. 

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

Executables are housed in `target` folder.

Default is to run debugging, flag for no debugging / best performance : `--release`

Flag to optimize for specific CPU architecture : `RUSTFLAGS="-C target-cpu=native"`

To find errors without actually making executable : `cargo check`
Cargo also offers a more indepth explaination of errors with either `rustc --explain {CODE}` or `cargo --explain {CODE}` where `CODE` is the error code provided by the compiler (ex. E0308).

Run `cargo --help` for additional information on using cargo. 

#### C-like version

This is for general use, current implementations use cargo.

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

[crates.io](https://crates.io/) is Rust's crate registry where you can search for crates and find links to documentation.

## Cargo - Rust's build system and package manager

Package dependencies are managed by Cargo. You simply list the packages with the version you want in a `Cargo.toml`. You can also manage paths to the source code, features (conditional compilation), and much more. A custom build script can be used for more complicated compile setups. An example is shown in `features/build.rs`. More information can be found in the [Cargo book](https://doc.rust-lang.org/cargo/).

## Versions

This folders houses four versions of SWM in Rust. `features` and `zip_ndarray` house serial versions of the code. `zip_ndarray_rayon` houses a shared memory version. `gpu-cudarc` houses a GPU version.

### Features

Features version is currently housed in `features`. This serial version uses nested for loops to updates arrays in the chosen type. Use features flags `box`, `vec`, `ndarray`, and `mdarray` to change array type. Default is to run with Box. For example, to run with Vec: 

`cargo run --features vec --release`

Information on the various array options:

| Array Type    | Dimensions | `for` loop speed | Crate needed | Docs |
| ------------- | ---------- | ---------------- | ------------ | ---- |
| Box           | 1          | fast             | none         | [doc](https://doc.rust-lang.org/std/boxed/struct.Box.html)
| Vec           | 1          | fast             | none         |[doc](https://doc.rust-lang.org/std/vec/struct.Vec.html)
| mdarray       | N          | slower           | `mdarray`    | [doc](https://docs.rs/mdarray/latest/mdarray/)
| ndarray       | N          | slowest          | `ndarray`    | [doc](https://docs.rs/ndarray/latest/ndarray/), [handbook](https://towardsdatascience.com/the-ultimate-ndarray-handbook-mastering-the-art-of-scientific-computing-with-rust-ef5ab767212a/)
| ndarray w/ zip | N          | fastest          | `ndarray`    | [zip doc](https://docs.rs/ndarray/latest/ndarray/struct.Zip.html)

### Zip ndarray

In `zip_ndarray`, ndarray array types are used and filled using the `zip` function in place of nested for loops. This has shown to be the fastest method in serial. Documentation is listed in the table above.

### Zip ndarray rayon

In `zip_ndarray_rayon`, the code from `zip_ndarray` is parallelized with shared memory using the `rayon` crate. This is done by replacing [`.for_each`](https://docs.rs/ndarray/latest/ndarray/struct.Zip.html#method.for_each) with [`.par_for_each`](https://docs.rs/ndarray/latest/ndarray/struct.Zip.html#method.par_for_each) in the zip function calls. The number of threads can by specified at compile time with

`RAYON_NUM_THREADS={num}`

The `rayon` crate documentation can be found [here](https://docs.rs/rayon/latest/rayon/). More information on the `rayon` crate feature for `ndarrray` can be found [here](https://docs.rs/ndarray/latest/ndarray/parallel/index.html).

More general information on Rust threading using the standard library can be found [here](https://doc.rust-lang.org/std/thread/). Much of this is called under the hood in `rayon`.

### GPU cudarc

In `gpu-cudarc`, the crate `cudarc` is used to run code on an NVIDIA GPU. Rust code is called from the CPU to run CUDA kernels on the GPU. More info on the crate can be found [here](https://github.com/chelsea0x3b/cudarc). Example 07 is used as a base for this implementation.

### Grid size

Grid size is a compile time argument. Default is 256x256. For example, to run on a 2048x2048 grid : 

`M=2048 N=2048 cargo run --release`

The crate [`const_env`](https://crates.io/crates/const_env) is used to configure compile-time constants in Rust.

## GPU Support

GPU support in Rust is in active development. Some crates are more stable and user-friendly than others. Here are some current options:

| Crate     | Kernel language | Portability  | Development Stage | Docs |
| --------- | --------------- | ------------ | ----------------- | ---- |
| cudarc    | CUDA C          | NVIDIA only  | stable            | [doc](https://docs.rs/cudarc/latest/cudarc/), [github](https://github.com/chelsea0x3b/cudarc)
| rust-cuda | Rust            | NVIDIA only  | early             | [doc](https://rust-gpu.github.io/rust-cuda/), [github](https://github.com/rust-gpu/rust-cuda)
| opencl3   | OpenCL C        | GPU portable | stable            | [doc](https://docs.rs/opencl3/latest/opencl3/), [github](https://github.com/kenba/opencl3)
| CubeCL    | Rust            | GPU portable | early             | [doc](https://docs.rs/cubecl/0.10.0/cubecl/), [github](https://github.com/tracel-ai/cubecl)

In general, the github repos were more helpful in understanding how the crate works, in particular looking at the `examples`.

A note for `rust-cuda`: attempts were made to run multiple example codes on Derecho. None were successful in compiling the code. It seemed that the `rustc_codegen_nvvm` rustc backend was the culprit for the compile errors. This is the crate that generates PTX code from the Rust kernels to run on the GPU. This binding seems to be the biggest barrier in getting native Rust implmentations to a user-friendly state.

## Other things of note

### Ownership rules

1. Each value in Rust has a owner.
1. There can only be one owner at a time.
1. When the owner goes out of scope, the value will be dropped.

*Scope* is the range within a code where a variable is valid. It typically looks like a code block defined by braces { }. This includes functions. A simple example is
```
{
   let a = 5;  // a is owner
   let b = a;  // b is owner
}
               // out of scope,
 	           // value dropped
```

### Mutability

In Rust, variables (and references) are immutable by default. You must explicitly state if a variable should be mutable, i.e. its value will change later in the code. A simple example is
```
let x = 5;      // immutable, cannot be changed

let mut x = 5;  // mutable, can be changed
x = 7;          // value is changed, old value dropped
```

Example of calling variables in a function:
```
func(&x)      // x is referenced but not changed
func(&mut x)  // x is a mutable reference, can be changed
```

### Type Array

The Array type in Rust is a bit misleading. This is not what we commonly think of as an array in scientific computing. `array` is a fixed size, one-dimensional, list of elements that is stored on the stack. This means overflow can occur at modest sizes, >~256 on a standard laptop. Use `vec` or `box` to create one-dimensional "arrays" that can be stored on the heap.

## Performance related notes

Cargo does support a test suite if needed, however I don't envision using that for performance testing work. 

`swm_rust.rs` currently contains an example for timing. This should be comparable to how walltime is computed in the C examples.

We will want to look for CUDA/GPU support crates eventually. See [NERSC docs](https://docs.nersc.gov/development/languages/rust/) for additional information on how to setup Rust to run efficiently on an HPC system. 
