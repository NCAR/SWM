// Test adding two vectors with malloc_device

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

void CopyVector(queue &q, int src[DOMAIN_SIZE], int dest[DOMAIN_SIZE]) {
    
    q.submit([&](handler &h) {
      h.memcpy(dest, &src[0], DOMAIN_SIZE * sizeof(int)); 
    });
}

void VecAdd(queue &q, range<DIM> R, int a_h[DOMAIN_SIZE], int b_h[DOMAIN_SIZE], int sum_h[DOMAIN_SIZE], int a_d[DOMAIN_SIZE], int b_d[DOMAIN_SIZE], int sum_d[DOMAIN_SIZE]) {
    
    CopyVector(q, a_h, a_d);
    CopyVector(q, b_h, b_d);
    q.wait();
    q.parallel_for(R, [=](auto i) { 
        sum_d[i] = a_d[i] + b_d[i]; 
    });
    q.wait();
    CopyVector(q, sum_d, sum_h);
    q.wait();
}

int main() {
    auto R = range<1>{DOMAIN_SIZE};
    host_selector h_selector;
    queue q_h(h_selector);
    default_selector d_selector;
    queue q_d(d_selector);
    std::cout << "Device: " << q_d.get_device().get_info<info::device::name>() << std::endl;
    
    // Allocating memory on the host
    int **u_h = malloc_host<int *>(3*DOMAIN_SIZE, q_h);
    int **v_h = malloc_host<int *>(3*DOMAIN_SIZE, q_h);
    int **p_h = malloc_host<int *>(3*DOMAIN_SIZE, q_h);
    for (int i=0; i<3; i++) {
        u_h[i] = malloc_host<int>(DOMAIN_SIZE, q_h);
        v_h[i] = malloc_host<int>(DOMAIN_SIZE, q_h);
        p_h[i] = malloc_host<int>(DOMAIN_SIZE, q_h);
    }
    
    // Allocating memory on the device
    // Allocating memory on the device
    int *u_d = malloc_device<int>(DOMAIN_SIZE, q_d);
    int *v_d = malloc_device<int>(DOMAIN_SIZE, q_d);
    int *p_d = malloc_device<int>(DOMAIN_SIZE, q_d);
    
    // Initialize arrays on the host
    auto e = q_h.parallel_for(R, [=](auto i) { 
        u_h[0][i] = i;
        v_h[0][i] = 2*i;
        
        u_h[1][i] = 3*i;
        v_h[1][i] = 6*i;
        
        u_h[2][i] = 5*i;
        v_h[2][i] = 10*i;
    });
    e.wait();
    
    VecAdd(q_d, R, u_h[0], v_h[0], p_h[0], u_d, v_d, p_d);
    VecAdd(q_d, R, u_h[1], v_h[1], p_h[1], u_d, v_d, p_d);
    VecAdd(q_d, R, u_h[2], v_h[2], p_h[2], u_d, v_d, p_d);
    
    for (int j=0; j<3; j++)
        for (int i=0; i<DOMAIN_SIZE; i++)
            std::cout << "p[" << j << "][" << i << "] = " << p_h[j][i] << std::endl;
    
    // Cleanup
    free(u_h, q_h);
    free(v_h, q_h);
    free(p_h, q_h);
    
    free(u_d, q_d);
    free(v_d, q_d);
    free(p_d, q_d);
    
    return 0;
}

/* 
Compile: 
    dpcpp -O2 -g -std=c++17 t016.cpp -o test016
Run
    ./test016 >test016-check.sh
Check
    diff t016.sh test016-check.sh
All
    rm test016; dpcpp -O2 -g -std=c++17 t016.cpp -o test016; ./test016 >test016-check.sh; diff t016.sh test016-check.sh
*/
