# AMReX Shallow Water Equations Mini-App

Based on the NCAR SWM mini-app found [here](https://github.com/NCAR/SWM). 

## Prerequisites
 - g++ (a version with support for c++20)
    - Other compilers should also work fine but I have been using gcc/g++ for testing so far.
 - make
 - [AMReX](https://github.com/AMReX-Codes/amrex)
 - [yt](https://yt-project.org/) (Only needed if you want to run the postprocessing script [plot_with_yt.py](plot_with_yt.py))

## Build

- Set AMREX_HOME
  - The make file requires this environment variable to be set. It should point to the directory where you pulled [AMReX](https://github.com/AMReX-Codes/amrex). I set mine in my shell startup script. You can also set it at the top of this file: [GNUmakefile](./GNUmakefile).

- make
    - If AMREX_HOME is set correctly then make should just work. This will produce and executable with a name similar to, main2d.gnu.ex . The exact name will vary depending on compile time options (compiler, MPI usage, openMP usage, Debugging options, etc.)

## Run
- Note that you need to supply [inputs](./inputs) file on command line when running the executable. For example:
    - ./main2d.gnu.ex inputs


## Postprocess
- python plot.py
    - Currently just generates images of plots the x velocity (u), y velocity (v), and pressure (p) for first and last time step. TODO - add check against the output from the other versions of the [mini-app](https://github.com/NCAR/SWM).


## Convenience Script
- [run_local.sh](./run_local.sh)
    - The environment variable SWM_AMREX_ROOT must be set to use this script. It should point to the top level of this directoy. This script will: 
        - Build. 
        - Create a subdirectory called SWM_AMREX_ROOT/run_dir if it does not already exist. 
        - Run and postprocess in that subdirectory. 
