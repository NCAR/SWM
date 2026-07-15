extern "C" __global__ void my_struct_kernel(double *data, const size_t n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < n)
    {

        data[i] += 1.0;
    }
}