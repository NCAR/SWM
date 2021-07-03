// Test adding two vectors with malloc_shared

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


void VectorAdd(queue &q, range<DIM> R, const int *a, const int *b, int *sum) {

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
    
    int *u = malloc_shared<int>(DOMAIN_SIZE, q);
    int *v = malloc_shared<int>(DOMAIN_SIZE, q);
    int *p = malloc_shared<int>(DOMAIN_SIZE, q);

    auto e = q.parallel_for(R, [=](auto i) { 
        u[i] = i;
        v[i] = 2*i;
    });
    
    VectorAdd(q, R, u, v, p);
    
    for (int i=0; i<DOMAIN_SIZE; i++)
        if (p[i] != 3*i)
            std::cout << "p[" << i << "] = " << p[i] << std::endl;
    
    free(u, q);
    free(v, q);
    free(p, q);
    
    return 0;
}

// Compile: dpcpp -O2 -g -std=c++17 t013.cpp -o test013
