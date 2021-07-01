// Define queue in the main function

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

void halo_update(queue q, const int m, const int n, double u[DOMAIN_SIZE]);

int main () {
    
    int m = M;
    int n = N;
    
    default_selector d_selector;
    queue q(d_selector);
    
    double *u = malloc_shared<double>(DOMAIN_SIZE, q);
    
    //double u[DOMAIN_SIZE];
    for (int i=0;i<DOMAIN_SIZE; i++) u[i] = i;
    

    halo_update(q, m, n, u);
    
    for (int i=0; i<DOMAIN_SIZE; i++)
        std::cout << "u[" << i << "] = " << u[i] << std::endl;
    
    return(0);
}

void halo_update(queue q, const int m, const int n, double u[DOMAIN_SIZE]){

        auto R = range<1>{DOMAIN_SIZE};
        
        int hsw = 0;
        int hne = (m+1)*(n+2)+n+1;
        int hse = n+1;
        int hnw = (m+1)*(n+2);
    
        int ine = m*(n+2)+n;
        int ise = 1*(n+2)+n;
        int isw = 1*(n+2)+1;
        int inw = m*(n+2)+1;
        
        auto e = q.parallel_for(R, [=](auto index) {
            int j = index%(N+2);
            int i = (int) (index - j)/(N+2);
            if (i==0 || j==0 || i == M+1 || j== N+1) {}
            // for (int j=1; j<n+1; j++)
            else {
                
                int hnorth = (m+1)*(n+2) + j;
                int hsouth = j;
                int inorth = m*(n+2) + j;
                int isouth = (n+2) + j;
        
                u[hsouth] = u[inorth];
                u[hnorth] = u[isouth];
                
                int hwest = i*(n+2);
                int heast = i*(n+2) + n+1;
                int iwest = i*(n+2) + 1;
                int ieast = i*(n+2) + n;
        
                u[hwest] = u[ieast];
                u[heast] = u[iwest];
 
                if (i==m) {
                    u[hsw]   = u[ine];
                    u[hnw]   = u[ise];
                    u[hne]   = u[isw];
                    u[hse]   = u[inw];
                }
            }
        });
}

/*
Compile
    dpcpp -O2 -g -std=c++17 t011.cpp -o test011
Run
    ./test011 >test011-check.sh
Check
    diff t007.sh test011-check.sh
All
dpcpp -O2 -g -std=c++17 t011.cpp -o test011; ./test011 >test011-check.sh; diff t007.sh test011-check.sh
*/
