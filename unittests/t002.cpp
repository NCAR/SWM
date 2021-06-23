// 2D range

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
    
    //int dim = 2;
    auto R = range<2>{M+2, N+2};
    
    double u[DOMAIN_SIZE];
    buffer<double, 2> u_buf(u, R);
    
    for (int i=0;i<DOMAIN_SIZE; i++) u[i] = 0.;
    
    default_selector d_selector;
    queue q(d_selector);
    
    
    
    q.submit([&](handler &h) {
        
        auto u = u_buf.get_access(h, write_only);
        
        h.parallel_for(R, [=](auto index) {
            
            u[index[0]][index[1]] =  index[0] + index[1];
            
        }); 
    });
    
    host_accessor u_read(u_buf, read_only);
    for (int i=0; i<M+2; i++)
        for (int j=0; j<N+2; j++)
            std::cout << "u[" << i << "][" << j << "] = " << u_read[i][j] << std::endl;
    
    return(0);
}

// Compile
// dpcpp -O2 -g -std=c++17 t002.cpp -o test002

// Here is the warning I get:
/*
In file included from t002.cpp:6:
In file included from /glob/development-tools/versions/oneapi/2021.2/inteloneapi/dev-utilities/2021.2.0/include/dpc_common.hpp:15:
In file included from /glob/development-tools/versions/oneapi/2021.2/inteloneapi/compiler/2021.2.0/linux/bin/../include/sycl/CL/sycl.hpp:11:
In file included from /glob/development-tools/versions/oneapi/2021.2/inteloneapi/compiler/2021.2.0/linux/bin/../include/sycl/CL/sycl/ONEAPI/atomic.hpp:11:
In file included from /glob/development-tools/versions/oneapi/2021.2/inteloneapi/compiler/2021.2.0/linux/bin/../include/sycl/CL/sycl/ONEAPI/atomic_accessor.hpp:14:
/glob/development-tools/versions/oneapi/2021.2/inteloneapi/compiler/2021.2.0/linux/bin/../include/sycl/CL/sycl/accessor.hpp:883:5: warning: loop not unrolled: the optimizer was unable to perform the requested transformation; the transformation might be disabled or specified as part of an unsupported transformation ordering [-Wpass-failed=transform-warning]
    for (int I = 0; I < AdjustedDim; ++I) {
    ^
1 warning generated.
*/

// NOTE: 
/*
A 2D range was an obstruction for good performance
*/