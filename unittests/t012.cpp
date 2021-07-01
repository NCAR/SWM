// 1D range - Create i and j inside the flattened array
//   This updates all nodes except for the ghost ones
//   parallelization inside update()

#include <vector>
#include <chrono>
#include <fstream>
#include <iostream>
#include <cmath>
#include "dpc_common.hpp"
#include <limits>
#include <CL/sycl.hpp>
#if FPGA || FPGA_EMULATOR
#include <CL/sycl/INTEL/fpga_extensions.hpp>
#endif

using namespace std;
using namespace sycl;

#define M 4
#define N 5
#define DOMAIN_SIZE (M+2)*(N+2)

int main () {
    
    default_selector d_selector;
    queue q(d_selector);
    
    double u[3][DOMAIN_SIZE];
    double *u_[3] = {u[0], u[1], u[2]};
    
    double **_u_ = malloc_shared<double *>(3*DOMAIN_SIZE, q);
    _u_ = u_;
    
    for (int i=0;i<DOMAIN_SIZE; i++) _u_[0][i] = i;
    for (int i=0;i<DOMAIN_SIZE; i++) _u_[1][i] = 2*i;
    for (int i=0;i<DOMAIN_SIZE; i++) _u_[2][i] = 3*i;
    
    for (int i=0; i<DOMAIN_SIZE; i++) {
        std::cout << "u[0][" << i << "] = " << u[0][i] << std::endl;
        std::cout << "u[1][" << i << "] = " << u[1][i] << std::endl;
        std::cout << "u[2][" << i << "] = " << u[2][i] << std::endl;
    }
    
    return(0);
}

/*
Compile
    dpcpp -O2 -g -std=c++17 t012.cpp -o test012
Run
    ./test012 >test012-check.sh
Check
    diff t012.sh test012-check.sh
All
dpcpp -O2 -g -std=c++17 t012.cpp -o test012; ./test012 >test012-check.sh; diff t012.sh test012-check.sh
*/