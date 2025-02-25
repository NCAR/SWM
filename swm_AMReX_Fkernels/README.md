# AMReX Shallow Water Equations Mini-App

Based on the NCAR SWM mini-app found [here](https://github.com/NCAR/SWM). 

## Prerequisites
 - A compiler for C, C++, and Fortran 
    - Just about any compilers should work but I have been using gcc/g++/gfortran for testing so far.
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
- Solution Verification
    - One time setup for building AMReX provided tool to diff plotfiles (fcompare). Full instructions on how to build AMReX postprocessing tools [here](https://amrex-codes.github.io/amrex/docs_html/Post_Processing.html).
        - cd $AMREX_HOME/Tools/Plotfile/
        - make
            - This should produce an executable: 
                - $AMREX_HOME/Tools/Plotfile/fcompare.gnu.ex
            - Note the exact name will vary depending on the compiler you used. The above example is when gcc was used.
    - Now you can diff two plot files using the fcompare tool:
        -    $AMREX_HOME/Tools/Plotfile/fcompare.gnu.ex plotfile1 plotfile2

- python plot.py
    - Currently just generates images of plots the x velocity (u), y velocity (v), and pressure (p) for first and last time step. TODO - add check against the output from the other versions of the [mini-app](https://github.com/NCAR/SWM).

## Convenience Script
- [run_local.sh](./run_local.sh)
    - The environment variable SWM_AMREX_ROOT must be set to use this script. It should point to this directory. This script will: 
        - Build. 
        - Create a subdirectory called SWM_AMREX_ROOT/run_dir if it does not already exist. 
        - Run in that subdirectory.
        - Verify the solution matches a reference solution. 
