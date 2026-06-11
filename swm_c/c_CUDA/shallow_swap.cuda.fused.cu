/*
 * SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
 * SPDX-License-Identifier: Apache-2.0
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/* CUDA version of shallow water equations benchmark
 * Converted from OpenACC version
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>
#include <nvtx3/nvToolsExt.h>

#define MIN(x,y) ((x)>(y)?(y):(x))
#define MAX(x,y) ((x)>(y)?(x):(y))

#define TRUE 1
#define FALSE 0
#define M 256
#define N 256
#define M_LEN (M + 1)
#define N_LEN (N + 1)
#define SIZE ((M_LEN)*(N_LEN))
#define ITMAX 4000
#define L_OUT TRUE
#define TILE_I 4
#define TILE_J 64

extern "C" double wtime();

#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__, \
                    cudaGetErrorString(err)); \
            exit(EXIT_FAILURE); \
        } \
    } while(0)

__global__ void init_psi_p(double* __restrict__ psi, double* __restrict__ p,
                           double a, double di, double dj, double pcf) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < SIZE) {
        int i = idx / N_LEN;
        int j = idx % N_LEN;
        psi[idx] = a * sin((i + 0.5) * di) * sin((j + 0.5) * dj);
        p[idx] = pcf * (cos(2.0 * i * di) + cos(2.0 * j * dj)) + 50000.0;
    }
}

__global__ void init_velocity(double* __restrict__ u, double* __restrict__ v,
                              const double* __restrict__ psi, double dx, double dy) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < M * N) {
        int i = idx / N;
        int j = idx % N;
        int idx01 = i * N_LEN + j + 1;
        int idx10 = (i + 1) * N_LEN + j;
        int idx11 = (i + 1) * N_LEN + j + 1;
        u[idx10] = -(psi[idx11] - psi[idx10]) / dy;
        v[idx01] = (psi[idx11] - psi[idx01]) / dx;
    }
}

__global__ void periodic_uv_j(double* __restrict__ u, double* __restrict__ v) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j < N) {
        u[j] = u[M * N_LEN + j];
        v[M * N_LEN + j + 1] = v[j + 1];
        if (j == 0) {
            u[N] = u[M * N_LEN];
            v[M * N_LEN] = v[N];
        }
    }
}

__global__ void periodic_uv_i(double* __restrict__ u, double* __restrict__ v) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < M) {
        u[(i + 1) * N_LEN + N] = u[(i + 1) * N_LEN];
        v[i * N_LEN] = v[i * N_LEN + N];
    }
}

__global__ void copy_old(double* __restrict__ uold, double* __restrict__ vold,
                         double* __restrict__ pold,
                         const double* __restrict__ u, const double* __restrict__ v,
                         const double* __restrict__ p) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < SIZE) {
        uold[idx] = u[idx];
        vold[idx] = v[idx];
        pold[idx] = p[idx];
    }
}

__global__ void periodic_new_j(double* __restrict__ unew, double* __restrict__ vnew,
                               double* __restrict__ pnew) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j < N) {
        unew[j] = unew[M * N_LEN + j];
        vnew[M * N_LEN + j + 1] = vnew[j + 1];
        pnew[M * N_LEN + j] = pnew[j];
        if (j == 0) {
            unew[N] = unew[M * N_LEN];
            vnew[M * N_LEN] = vnew[N];
            pnew[M * N_LEN + N] = pnew[0];
        }
    }
}

__global__ void periodic_new_i(double* __restrict__ unew, double* __restrict__ vnew,
                               double* __restrict__ pnew) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < M) {
        unew[(i + 1) * N_LEN + N] = unew[(i + 1) * N_LEN];
        vnew[i * N_LEN] = vnew[i * N_LEN + N];
        pnew[i * N_LEN + N] = pnew[i * N_LEN];
    }
}

__device__ __forceinline__ int idx(int i, int j) { return i * N_LEN + j; }

// For cu: row 0 copies from row M, col N copies from col 0
// cu is computed at rows 1..M, cols 0..N-1, then boundaries are copied
__device__ double compute_cu_at(int row, int col, const double* __restrict__ p, const double* __restrict__ u) {
    // Handle periodic: row 0 -> use formula for row M, col N -> use formula for col 0
    int eff_row = (row == 0) ? M : row;
    int eff_col = (col == N) ? 0 : col;
    return 0.5 * (p[idx(eff_row, eff_col)] + p[idx(eff_row - 1, eff_col)]) * u[idx(eff_row, eff_col)];
}

// For cv: row M copies from row 0, col 0 copies from col N
// cv is computed at rows 0..M-1, cols 1..N, then boundaries are copied
__device__ double compute_cv_at(int row, int col, const double* __restrict__ p, const double* __restrict__ v) {
    int eff_row = (row == M) ? 0 : row;
    int eff_col = (col == 0) ? N : col;
    return 0.5 * (p[idx(eff_row, eff_col)] + p[idx(eff_row, eff_col - 1)]) * v[idx(eff_row, eff_col)];
}

// For z: row 0 copies from row M, col 0 copies from col N
// z is computed at rows 1..M, cols 1..N, then boundaries are copied
__device__ double compute_z_at(int row, int col, const double* __restrict__ p,
                                const double* __restrict__ u, const double* __restrict__ v,
                                double fsdx, double fsdy) {
    int eff_row = (row == 0) ? M : row;
    int eff_col = (col == 0) ? N : col;
    int r_m1 = eff_row - 1;
    int c_m1 = eff_col - 1;
    double numer = fsdx * (v[idx(eff_row, eff_col)] - v[idx(r_m1, eff_col)]) -
                   fsdy * (u[idx(eff_row, eff_col)] - u[idx(eff_row, c_m1)]);
    double denom = p[idx(r_m1, c_m1)] + p[idx(eff_row, c_m1)] + p[idx(eff_row, eff_col)] + p[idx(r_m1, eff_col)];
    return numer / denom;
}

// For h: row M copies from row 0, col N copies from col 0
// h is computed at rows 0..M-1, cols 0..N-1, then boundaries are copied
__device__ double compute_h_at(int row, int col, const double* __restrict__ p,
                                const double* __restrict__ u, const double* __restrict__ v) {
    int eff_row = (row == M) ? 0 : row;
    int eff_col = (col == N) ? 0 : col;
    int r_p1 = eff_row + 1;
    int c_p1 = eff_col + 1;
    double u_p1 = u[idx(r_p1, eff_col)];
    double u_00 = u[idx(eff_row, eff_col)];
    double v_p1 = v[idx(eff_row, c_p1)];
    double v_00 = v[idx(eff_row, eff_col)];
    return p[idx(eff_row, eff_col)] + 0.25 * (u_p1*u_p1 + u_00*u_00 + v_p1*v_p1 + v_00*v_00);
}

// Fused kernel: computes cucvzh + new + time_smooth in one pass
__global__ void compute_fused_smooth(double* __restrict__ unew, double* __restrict__ vnew,
                                     double* __restrict__ pnew,
                                     double* __restrict__ uold, double* __restrict__ vold,
                                     double* __restrict__ pold,
                                     const double* __restrict__ p, const double* __restrict__ u,
                                     const double* __restrict__ v,
                                     double fsdx, double fsdy, double tdts8, double tdtsdx, double tdtsdy,
                                     double alpha, int do_smooth) {
    int j = blockIdx.x * TILE_J + threadIdx.x;
    int i = blockIdx.y * TILE_I + threadIdx.y;

    if (i < M && j < N) {
        double cu_00 = compute_cu_at(i, j, p, u);
        double cu_01 = compute_cu_at(i, j + 1, p, u);
        double cu_10 = compute_cu_at(i + 1, j, p, u);
        double cu_11 = compute_cu_at(i + 1, j + 1, p, u);

        double cv_00 = compute_cv_at(i, j, p, v);
        double cv_01 = compute_cv_at(i, j + 1, p, v);
        double cv_10 = compute_cv_at(i + 1, j, p, v);
        double cv_11 = compute_cv_at(i + 1, j + 1, p, v);

        double z_10 = compute_z_at(i + 1, j, p, u, v, fsdx, fsdy);
        double z_11 = compute_z_at(i + 1, j + 1, p, u, v, fsdx, fsdy);
        double z_01 = compute_z_at(i, j + 1, p, u, v, fsdx, fsdy);

        double h_00 = compute_h_at(i, j, p, u, v);
        double h_10 = compute_h_at(i + 1, j, p, u, v);
        double h_01 = compute_h_at(i, j + 1, p, u, v);

        int idx00 = i * N_LEN + j;
        int idx01 = i * N_LEN + j + 1;
        int idx10 = (i + 1) * N_LEN + j;

        double unew_val = uold[idx10] + tdts8 * (z_11 + z_10) * (cv_11 + cv_01 + cv_00 + cv_10) -
                          tdtsdx * (h_10 - h_00);
        double vnew_val = vold[idx01] - tdts8 * (z_11 + z_01) * (cu_11 + cu_01 + cu_00 + cu_10) -
                          tdtsdy * (h_01 - h_00);
        double pnew_val = pold[idx00] - tdtsdx * (cu_10 - cu_00) - tdtsdy * (cv_01 - cv_00);

        // Write new values (needed for pointer swap)
        unew[idx10] = unew_val;
        vnew[idx01] = vnew_val;
        pnew[idx00] = pnew_val;

        // Fused time_smooth: update old values using the new values we just computed
        if (do_smooth) {
            uold[idx10] = u[idx10] + alpha * (unew_val - 2.0 * u[idx10] + uold[idx10]);
            vold[idx01] = v[idx01] + alpha * (vnew_val - 2.0 * v[idx01] + vold[idx01]);
            pold[idx00] = p[idx00] + alpha * (pnew_val - 2.0 * p[idx00] + pold[idx00]);
        }
    }
}

// Copy kernel for ncycle == 1 (when time_smooth is not called)
__global__ void copy_old_fused(double* __restrict__ uold, double* __restrict__ vold,
                               double* __restrict__ pold,
                               const double* __restrict__ u, const double* __restrict__ v,
                               const double* __restrict__ p) {
    int j = blockIdx.x * TILE_J + threadIdx.x;
    int i = blockIdx.y * TILE_I + threadIdx.y;

    if (i < M && j < N) {
        int idx00 = i * N_LEN + j;
        int idx01 = i * N_LEN + j + 1;
        int idx10 = (i + 1) * N_LEN + j;

        uold[idx10] = u[idx10];
        vold[idx01] = v[idx01];
        pold[idx00] = p[idx00];
    }
}

int main(int argc, char **argv) {
    double *u, *v, *p, *unew, *vnew, *pnew;
    double *d_u, *d_v, *d_p, *d_unew, *d_vnew, *d_pnew;
    double *d_uold, *d_vold, *d_pold;
    double *d_psi;

    u = (double *)malloc(sizeof(double) * SIZE);
    v = (double *)malloc(sizeof(double) * SIZE);
    p = (double *)malloc(sizeof(double) * SIZE);
    unew = (double *)calloc(SIZE, sizeof(double));
    vnew = (double *)calloc(SIZE, sizeof(double));
    pnew = (double *)calloc(SIZE, sizeof(double));

    CUDA_CHECK(cudaMalloc(&d_u, sizeof(double) * SIZE));
    CUDA_CHECK(cudaMalloc(&d_v, sizeof(double) * SIZE));
    CUDA_CHECK(cudaMalloc(&d_p, sizeof(double) * SIZE));
    CUDA_CHECK(cudaMalloc(&d_unew, sizeof(double) * SIZE));
    CUDA_CHECK(cudaMalloc(&d_vnew, sizeof(double) * SIZE));
    CUDA_CHECK(cudaMalloc(&d_pnew, sizeof(double) * SIZE));
    CUDA_CHECK(cudaMalloc(&d_uold, sizeof(double) * SIZE));
    CUDA_CHECK(cudaMalloc(&d_vold, sizeof(double) * SIZE));
    CUDA_CHECK(cudaMalloc(&d_pold, sizeof(double) * SIZE));
    CUDA_CHECK(cudaMalloc(&d_psi, sizeof(double) * SIZE));

    double dt, tdt, dx, dy, a, alpha, el, pi;
    double tpi, di, dj, pcf;
    double tdts8, tdtsdx, tdtsdy, fsdx, fsdy;

    int mnmin, ncycle;
    int i;

    double mfs100, mfs200, mfs300;
    double t100, t200, t300;
    double tstart, tcyc, time, ptime;
    double c1, c2;

    dt = 90.;
    tdt = dt;

    dx = 100000.;
    dy = 100000.;
    fsdx = 4. / dx;
    fsdy = 4. / dy;

    a = 1000000.;
    alpha = .001;

    el = N * dx;
    pi = 4. * atan(1.);
    tpi = pi + pi;
    di = tpi / M;
    dj = tpi / N;
    pcf = pi * pi * a * a / (el * el);

    int blockSize = 256;
    int numBlocksSize = (SIZE + blockSize - 1) / blockSize;
    int numBlocksMN = (M * N + blockSize - 1) / blockSize;
    int numBlocksN = (N + blockSize - 1) / blockSize;
    int numBlocksM = (M + blockSize - 1) / blockSize;

    // 2D grid and block for tiled kernels:
    dim3 blockTiled(TILE_J, TILE_I);
    dim3 gridTiled((N + TILE_J - 1) / TILE_J, (M + TILE_I - 1) / TILE_I);


    CUDA_CHECK(cudaMemset(d_unew, 0, sizeof(double) * SIZE));
    CUDA_CHECK(cudaMemset(d_vnew, 0, sizeof(double) * SIZE));
    CUDA_CHECK(cudaMemset(d_pnew, 0, sizeof(double) * SIZE));

    init_psi_p<<<numBlocksSize, blockSize>>>(d_psi, d_p, a, di, dj, pcf);
    CUDA_CHECK(cudaMemcpy(p, d_p, sizeof(double) * SIZE, cudaMemcpyDeviceToHost));

    init_velocity<<<numBlocksMN, blockSize>>>(d_u, d_v, d_psi, dx, dy);

    periodic_uv_j<<<numBlocksN, blockSize>>>(d_u, d_v);
    periodic_uv_i<<<numBlocksM, blockSize>>>(d_u, d_v);

    CUDA_CHECK(cudaMemcpy(u, d_u, sizeof(double) * SIZE, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(v, d_v, sizeof(double) * SIZE, cudaMemcpyDeviceToHost));

    copy_old<<<numBlocksSize, blockSize>>>(d_uold, d_vold, d_pold, d_u, d_v, d_p);
    CUDA_CHECK(cudaDeviceSynchronize());

    if (L_OUT) {
        printf(" number of points in the x direction %d\n", N);
        printf(" number of points in the y direction %d\n", M);
        printf(" grid spacing in the x direction     %f\n", dx);
        printf(" grid spacing in the y direction     %f\n", dy);
        printf(" time step                           %f\n", dt);
        printf(" time filter parameter               %f\n", alpha);

        mnmin = MIN(M, N);
        printf(" initial diagonal elements of p\n");
        for (i = 0; i < mnmin; i++) {
            printf("%f ", p[i * N_LEN + i]);
        }
        printf("\n initial diagonal elements of u\n");
        for (i = 0; i < mnmin; i++) {
            printf("%f ", u[i * N_LEN + i]);
        }
        printf("\n initial diagonal elements of v\n");
        for (i = 0; i < mnmin; i++) {
            printf("%f ", v[i * N_LEN + i]);
        }
        printf("\n");
    }

    tstart = wtime();
    time = 0.;
    t100 = 0.;
    t200 = 0.;
    t300 = 0.;

    for (ncycle = 1; ncycle <= ITMAX; ncycle++) {

        tdts8 = tdt / 8.;
        tdtsdx = tdt / dx;
        tdtsdy = tdt / dy;

        int do_smooth = (ncycle > 1) ? 1 : 0;

        nvtxRangePush("FullyFusedKernel");
        c1 = wtime();

        // Fully fused kernel: cucvzh + new + time_smooth all in one
        compute_fused_smooth<<<gridTiled, blockTiled>>>(d_unew, d_vnew, d_pnew,
                                                         d_uold, d_vold, d_pold,
                                                         d_p, d_u, d_v,
                                                         fsdx, fsdy, tdts8, tdtsdx, tdtsdy,
                                                         alpha, do_smooth);

        periodic_new_j<<<numBlocksN, blockSize>>>(d_unew, d_vnew, d_pnew);
        periodic_new_i<<<numBlocksM, blockSize>>>(d_unew, d_vnew, d_pnew);

        CUDA_CHECK(cudaDeviceSynchronize());
        c2 = wtime();
        t100 = t100 + (c2 - c1);
        nvtxRangePop();

        time = time + dt;

        if (ncycle == 1) {
            tdt = tdt + tdt;
            // For ncycle==1, copy u/v/p to uold/vold/pold at the positions compute_fused writes
            copy_old_fused<<<gridTiled, blockTiled>>>(d_uold, d_vold, d_pold, d_u, d_v, d_p);
            CUDA_CHECK(cudaDeviceSynchronize());
        }

        // Swap pointers
        double *tmp;
        tmp = d_u; d_u = d_unew; d_unew = tmp;
        tmp = d_v; d_v = d_vnew; d_vnew = tmp;
        tmp = d_p; d_p = d_pnew; d_pnew = tmp;
    }

    double tstop = wtime();
    double loop_time = tstop - tstart;

    CUDA_CHECK(cudaMemcpy(unew, d_unew, sizeof(double) * SIZE, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(vnew, d_vnew, sizeof(double) * SIZE, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(pnew, d_pnew, sizeof(double) * SIZE, cudaMemcpyDeviceToHost));

    if (L_OUT) {
        ptime = time / 3600.;
        printf(" cuda cycle number %d model time in hours %f\n", ITMAX, ptime);
        printf(" diagonal elements of p\n");
        for (i = 0; i < mnmin; i++) {
            printf("%f ", pnew[i * N_LEN + i]);
        }
        printf("\n diagonal elements of u\n");
        for (i = 0; i < mnmin; i++) {
            printf("%f ", unew[i * N_LEN + i]);
        }
        printf("\n diagonal elements of v\n");
        for (i = 0; i < mnmin; i++) {
            printf("%f ", vnew[i * N_LEN + i]);
        }
        printf("\n");

        mfs100 = 0.0;
        mfs200 = 0.0;
        mfs300 = 0.0;
        if (t100 > 0) { mfs100 = ITMAX * 24. * M * N / t100 / 1000000; }
        if (t200 > 0) { mfs200 = ITMAX * 26. * M * N / t200 / 1000000; }
        if (t300 > 0) { mfs300 = ITMAX * 15. * M * N / t300 / 1000000; }

        tcyc = loop_time / ITMAX;

        printf(" cycle number %d total computer time %f time per cycle %f\n", ITMAX, loop_time, tcyc);
        printf(" time and megaflops for loop 100 %.6f %.6f\n", t100, mfs100);
        printf(" time and megaflops for loop 200 %.6f %.6f\n", t200, mfs200);
        printf(" time and megaflops for loop 300 %.6f %.6f\n", t300, mfs300);
    }

    free(u);
    free(v);
    free(p);
    free(unew);
    free(vnew);
    free(pnew);

    cudaFree(d_u);
    cudaFree(d_v);
    cudaFree(d_p);
    cudaFree(d_unew);
    cudaFree(d_vnew);
    cudaFree(d_pnew);
    cudaFree(d_uold);
    cudaFree(d_vold);
    cudaFree(d_pold);
    cudaFree(d_psi);

    return 0;
}
