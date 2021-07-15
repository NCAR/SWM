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
constexpr size_t  DOMAIN_SIZE = M_LEN*N_LEN;
#define DIM 1

int main() {
    auto R = range<1>{DOMAIN_SIZE};
    default_selector d_selector;
    queue q(d_selector);

    int **u = malloc_device<int *>(DOMAIN_SIZE, q);
    int **v = malloc_device<int *>(DOMAIN_SIZE, q);
    int **p = malloc_device<int *>(DOMAIN_SIZE, q);
    for(int i=0;i<3;i++) {
            u[i] = malloc_device<int>(DOMAIN_SIZE, q);
            v[i] = malloc_device<int>(DOMAIN_SIZE, q);
            p[i] = malloc_device<int>(DOMAIN_SIZE, q);
    }
   free(u,q);
   free(v,q);
   free(p,q);
   return 0;
}

/* 
Compile: 
    dpcpp -O2 -g -std=c++17 t017.cpp -o test017
*/