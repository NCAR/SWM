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
    
    int u_[3][DOMAIN_SIZE]; int *_u_[3] = {u_[0], u_[1], u_[2]}; u = _u_;
    int v_[3][DOMAIN_SIZE]; int *_v_[3] = {v_[0], v_[1], v_[2]}; v = _v_;
    int p_[3][DOMAIN_SIZE]; int *_p_[3] = {p_[0], p_[1], p_[2]}; p = _p_;
    
    auto e = q.parallel_for(R, [=](auto i) { 
        u[0][i] = i;
        v[0][i] = 2*i;
    });
    
    VecAdd(q, R, u[0], v[0], p[0]);
    
    for (int i=0; i<DOMAIN_SIZE; i++)
        std::cout << "p[0][" << i << "] = " << p[0][i] << std::endl;
    
    free(u, q);
    free(v, q);
    free(p, q);
    
    return 0;
}

// Compile: dpcpp -O2 -g -std=c++17 t014.cpp -o test014
