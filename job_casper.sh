
#!/bin/bash
 
module purge
module load intel/2021.2 ncarenv/1.3 ncarcompilers/0.5.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/glade/u/apps/dav/opt/gnu/9.1.0/lib64
 
source /glade/u/apps/opt/intel/2021.2/setvars.sh > /dev/null 2>&1
cd $PBS_O_WORKDIR
 
echo "##" $(whoami) is compiling DPCPP_Essentials Module1 -- oneAPI Intro sample - simple-vector-incr.cpp
d:pcpp --gcc-toolchain=/glade/u/apps/dav/opt/gnu/9.1.0 -std=c++17 simple-vector-incr.cpp -o a.out
 
if [ $? -eq 0 ]; then ./a.out ; fi
