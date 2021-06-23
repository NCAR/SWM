// Create buffer and accessor
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
    
    // Option 1
    double u[DOMAIN_SIZE];
    buffer<double, 1> u_buf(u, range<1>{DOMAIN_SIZE});

    // Option 2
    //vector<double> u; u.resize(DOMAIN_SIZE);
    //buffer u_buf(u);
    
    for (int i=0;i<DOMAIN_SIZE; i++) u[i] = 0.;
    
    default_selector d_selector;
    queue q(d_selector);
    
    
    
    q.submit([&](handler &h) {
        
        //auto u = u_buf.get_access<access::mode::read_write>(h);
        auto u = u_buf.get_access(h, write_only);
        //accessor u(u_buf, h, write_only);
        
        h.parallel_for(range(DOMAIN_SIZE), [=](auto m_) {
            
            u[m_] =  m_;
            
        }); 
    });
    
    host_accessor u_read(u_buf, read_only);
    for (int i=0; i<DOMAIN_SIZE; i++)
        std::cout << "u[" << i << "] = " << u_read[i] << std::endl;
    
    return(0);
}

// Compile
// dpcpp -O2 -g -std=c++17 t001.cpp -o test001