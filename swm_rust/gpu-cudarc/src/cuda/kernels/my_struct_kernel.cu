extern "C" __global__ void my_struct_kernel(double *data, const size_t TOT_LEN)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < TOT_LEN)
    {

        data[i] += 1.0;
    }
}

extern "C" __global__ void init_conds(
    double *u,
    double *v,
    double *p,
    double *psi,
    double dx,
    double dy,
    double a,
    int M_LEN,
    int N_LEN,
    int TOT_LEN,
    double di,
    double dj,
    double pcf
)
{
    // index
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < TOT_LEN)
    {
        int i = idx / N_LEN;
        int j = idx % N_LEN;

        double ii = (double)i;
        double jj = (double)j;

        // initialize stream function psi and pressure p
        psi[idx] = a * sin((ii + 0.5) * di) * sin((jj + 0.5) * dj);
        p[idx] = pcf * ( cos(2.0 * ii * di) + cos(2.0 * jj * dj) ) + 50000.0;

        // initialize velocities u and v
        int idx01 = (i*N_LEN) + j+1; //[i][j+1]
        int idx10 = ((i+1)*N_LEN) + j; //[i+1][j]
        int idx11 = ((i+1)*N_LEN) + j+1; //[i+1][j+1]
        u[idx10] = -(psi[idx11] - psi[idx10]) / dy;
        v[idx01] = (psi[idx11] - psi[idx01]) / dx;
    }

}

extern "C" __global__ void apply_uv_bcs(
    double *u,
    double *v,
    int M,
    int N,
    int N_LEN
)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Apply left/right periodic boundaries
    if (idx < N)
    {
        int j = idx;

        u[j] = u[M * N_LEN + j];
        v[M * N_LEN + (j + 1)] = v[j + 1];
    }

    // Apply top/bottom periodic boundaries
    if (idx < M)
    {
        int i = idx;

        u[(i + 1) * N_LEN + N] = u[(i + 1) * N_LEN];
        v[i * N_LEN] = v[i * N_LEN + N];
    }

    // Corner values
    if (idx == 0)
    {
        u[N] = u[M * N_LEN];
        v[M * N_LEN] = v[N];
    }
}

extern "C" __global__ void init_olds(
    double *uold,
    double *vold,
    double *pold,
    const double *u,
    const double *v,
    const double *p,
    int TOT_LEN
)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < TOT_LEN)
    {
        uold[idx] = u[idx];
        vold[idx] = v[idx];
        pold[idx] = p[idx];
    }
}