// Test adding two vectors with malloc_shared called 4 times for each 2d array

#include <CL/sycl.hpp>
#include <array>
#include <iostream>
#if FPGA || FPGA_EMULATOR
#include <CL/sycl/INTEL/fpga_extensions.hpp>
#endif

using namespace sycl;

#define M 4
#define N 5
#define M_LEN (M + 2)
#define N_LEN (N + 2)
#define DOMAIN_SIZE M_LEN*N_LEN
#define DIM 1


void VecAdd(queue &q, range<DIM> R, const int a[DOMAIN_SIZE], const int b[DOMAIN_SIZE], int sum[DOMAIN_SIZE]) {

  auto e = q.parallel_for(R, [=](auto i) { 
      sum[i] = a[i] + b[i]; 
  });

  e.wait();
}

int main() {
    auto R = range<1>{DOMAIN_SIZE};
    default_selector d_selector;
    queue q(d_selector);
    std::cout << "Device: " << q.get_device().get_info<info::device::name>() << std::endl;
    
    int **u = malloc_shared<int *>(3*DOMAIN_SIZE, q);
    int **v = malloc_shared<int *>(3*DOMAIN_SIZE, q);
    int **p = malloc_shared<int *>(3*DOMAIN_SIZE, q);
    
    for (int i=0; i<3; i++) {
        u[i] = malloc_shared<int>(DOMAIN_SIZE, q);
        v[i] = malloc_shared<int>(DOMAIN_SIZE, q);
        p[i] = malloc_shared<int>(DOMAIN_SIZE, q);
    }
    
    auto e = q.parallel_for(R, [=](auto i) { 
        u[0][i] = i;
        v[0][i] = 2*i;
        
        u[1][i] = 3*i;
        v[1][i] = 6*i;
        
        u[2][i] = 5*i;
        v[2][i] = 10*i;
    });
    
    VecAdd(q, R, u[0], v[0], p[0]);
    VecAdd(q, R, u[1], v[1], p[1]);
    VecAdd(q, R, u[2], v[2], p[2]);
    
    for (int j=0; j<3; j++)
        for (int i=0; i<DOMAIN_SIZE; i++)
            std::cout << "p[" << j << "][" << i << "] = " << p[j][i] << std::endl;
    
    free(u, q);
    free(v, q);
    free(p, q);
    
    return 0;
}

// Compile: dpcpp -O2 -g -std=c++17 t015.cpp -o test015

