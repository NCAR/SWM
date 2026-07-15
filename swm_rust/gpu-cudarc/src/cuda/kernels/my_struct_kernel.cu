extern "C" __global__ void my_struct_kernel(double *data, const size_t n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < n)
    {
        int idx = 4 * i;

        data[idx + 0] += 1.0;
        data[idx + 1] += 1.0;
        data[idx + 2] += 1.0;
        data[idx + 3] += 1.0;
    }
}