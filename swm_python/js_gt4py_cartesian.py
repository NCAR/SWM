import gt4py.cartesian.gtscript as gtscript
import gt4py.next as gtx
from gt4py.next import Field
import numpy as np
import time

dtype = np.float64

M = 3000      # row dimension
N = 4000     # column dimension
nx = M + 1
ny = N + 1
nz = 1
coeff1 = dtype(0.5)
coeff2 = dtype(0.25)
dx = 100000.
dy = 100000.
fsdx = dtype(4. / dx)
fsdy = dtype(4. / dy)

allocator = gtx.gtfn_cpu       # gtx.gtfn_cpu, gtx.gtfn_gpu
backend = "gt:cpu_ifirst"    # "gt:gpu" # "gt:cpu_ifirst"

###############################
# Naive Python implementation #
###############################
 
start_time = time.perf_counter()

p = np.ones((nx, ny))
u = np.ones((nx, ny))
v = np.ones((nx, ny))
cu_ori = np.zeros((nx, ny))
cv_ori = np.zeros((nx, ny))
z_ori = np.ones((nx, ny))
h_ori = np.zeros((nx, ny))
 
for i in range(M):
    for j in range(N):
        cu_ori[i + 1, j] = coeff1 * (p[i + 1, j] + p[i, j]) * u[i + 1, j]
        cv_ori[i, j + 1] = coeff1 * (p[i, j + 1] + p[i, j]) * v[i, j + 1]
        z_ori[i + 1, j + 1] = ( fsdx * (v[i + 1, j + 1] - v[i, j + 1]) -
                                fsdy * (u[i + 1, j + 1] - u[i, j + 1]) 
                              ) / (p[i, j] + p[i + 1, j] + p[i + 1, j + 1] + p[i, j + 1])
        h_ori[i, j] = p[i, j] + coeff2 * (u[i + 1, j] * u[i + 1, j] +
                                         u[i, j] * u[i, j] +
                                         v[i, j + 1] * v[i, j + 1] +
                                         v[i, j] * v[i, j])

end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Naive python, elapsed time: {elapsed_time} seconds")

##################################
# GT4PY Cartesian implementation #
##################################

start_time = time.perf_counter()

I = gtx.Dimension("I")
J = gtx.Dimension("J")
K = gtx.Dimension("K", kind=gtx.DimensionKind.VERTICAL)
 
domain = gtx.domain({I: nx, J: ny, K: nz})
 
p_cart = gtx.ones(domain, dtype, allocator=allocator)
u_cart = gtx.as_field(domain, np.ones((nx,ny,nz),dtype=dtype), dtype, allocator=allocator)
v_cart = gtx.as_field(domain, np.ones((nx,ny,nz),dtype=dtype), dtype, allocator=allocator)
cu_cart = gtx.zeros(domain, dtype, allocator=allocator)
cv_cart = gtx.zeros(domain, dtype, allocator=allocator)
z_cart = gtx.ones(domain, dtype, allocator=allocator)
h_cart = gtx.zeros(domain, dtype, allocator=allocator)

@gtscript.stencil(backend=backend, rebuild=True)
def cart_calc_region(
    p: gtscript.Field[dtype],
    u: gtscript.Field[dtype],
    v: gtscript.Field[dtype],
    cu: gtscript.Field[dtype],
    cv: gtscript.Field[dtype],
    z: gtscript.Field[dtype],
    h: gtscript.Field[dtype],
    coeff1: dtype,
    fsdx: dtype,
    fsdy: dtype,
    coeff2: dtype
):
    with computation(PARALLEL), interval(...):
        with horizontal(region[1:,:-1]):    
            cu = coeff1 * (p + p[I-1]) * u
        with horizontal(region[:-1,1:]):
            cv = coeff1 * (p + p[J-1]) * v
        with horizontal(region[1:,1:]):
            z = (fsdx * (v - v[I-1]) - fsdy * (u - u[I-1])) / (p[I-1,J-1] + p[J-1] + p + p[I-1])
        with horizontal(region[:-1,:-1]):
            h = p + coeff2 * (u[I+1] * u[I+1] + u * u + v[J+1] * v[J+1] + v * v)

# compute cu, cv, z and h; domain here refers to the number of points in each dimension to be computed
cart_calc_region(p_cart, u_cart, v_cart, cu_cart, cv_cart, z_cart, h_cart, 
                 coeff1, fsdx, fsdy, coeff2, origin=(0,0,0), domain=(nx,ny,nz))

end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"GT4PY Cartesian, elapsed time: {elapsed_time} seconds")

####################################
# Check the GT4PY Cartesian result #
####################################

assert np.array_equal(cu_cart[:,:,0].asnumpy(), cu_ori)
assert np.array_equal(cv_cart[:,:,0].asnumpy(), cv_ori)
assert np.array_equal(z_cart[:,:,0].asnumpy(), z_ori)
assert np.array_equal(h_cart[:,:,0].asnumpy(), h_ori)
