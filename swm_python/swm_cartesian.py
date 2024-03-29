import numpy as np
import numpy as np
import gt4py.next as gtx
import gt4py.cartesian.gtscript as gtscript
from time import perf_counter
#import cupy
#import gt4py
import initial_conditions
import utils
import config

I = gtx.Dimension("I")
J = gtx.Dimension("J")
K = gtx.Dimension("K", kind = gtx.DimensionKind.VERTICAL)
    
dtype = np.float64

cartesian_backend = "numpy"
# cartesian_backend = "gt:gpu"

@gtscript.stencil(backend=cartesian_backend)
def calc_h(
    p: gtscript.Field[dtype],
    u: gtscript.Field[dtype],
    v: gtscript.Field[dtype],
    h: gtscript.Field[dtype]
):
    with computation(PARALLEL), interval(...):
        h = p + 0.25 * (u[1,0,0] * u[1,0,0] + u * u + v[0,1,0] * v[0,1,0] + v * v)

@gtscript.stencil(backend=cartesian_backend)
def update_bndry(
    old:   gtscript.Field[dtype],
    new:   gtscript.Field[dtype],
    north:  int,
    south:  int,
    west:   int,
    east:   int):
    with computation(PARALLEL), interval(...):
         with horizontal(region[:,0,0]):
              if(south==1):
                new = old[0,1,0]
         with horizontal(region[:,16,0]):
              if(north==1):
                new = old[0,-1,0]
         with horizontal(region[0,:,0]):
              if(west==1):
                new = old[1,0,0]
         with horizontal(region[16,:,0]):
              if(east==1):
                new = old[-1,0,0]


@gtscript.stencil(backend=cartesian_backend)
def calc_z(
    fsdx: float,
    fsdy: float,
    u: gtscript.Field[dtype],
    v: gtscript.Field[dtype],
    p: gtscript.Field[dtype],
    z: gtscript.Field[dtype]
):
    with computation(PARALLEL), interval(...):
        z = (fsdx * (v - v[-1,0,0]) - fsdy * (u - u[0,-1,0])) / (p[-1,-1,0] + p[0,-1,0] + p + p[-1,0,0])
        
@gtscript.stencil(backend=cartesian_backend)
def calc_cu(
    u: gtscript.Field[dtype],
    p: gtscript.Field[dtype],
    cu: gtscript.Field[dtype]
):
    with computation(PARALLEL), interval(...):
        cu = .5 * (p + p[-1,0,0]) * u
        
@gtscript.stencil(backend=cartesian_backend)
def calc_cv(
    v: gtscript.Field[dtype],
    p: gtscript.Field[dtype],
    cv: gtscript.Field[dtype]
):
    with computation(PARALLEL), interval(...):
        cv = .5 * (p[0,-1,0] + p) * v

#recheck this section
#pnew[i,j,0] = pold[i,j,0] - tdtsdx * (cu[i+1,j,0] - cu[i,j,0]) - tdtsdy * (cv[i,j+1,0] - cv[i,j,0])
@gtscript.stencil(backend=cartesian_backend)
def calc_pnew(
    tdtsdx: float,
    tdtsdy: float,
    pold: gtscript.Field[dtype],
    cu: gtscript.Field[dtype],
    cv: gtscript.Field[dtype],
    pnew: gtscript.Field[dtype]
):
    with computation(PARALLEL), interval(...):
        pnew = pold - tdtsdx * (cu[1,0,0] - cu) - tdtsdy * (cv[0,1,0] - cv)

#unew[i+1,j,0] = uold[i+1,j,0] + tdts8 * (z[i+1,j+1,0] + z[i+1,j,0]) * (cv[i+1,j+1,0] + cv[i+1,j,0] + cv[i,j+1,0] + cv[i,j,0]) - tdtsdx * (h[i+1,j,0] - h[i,j,0])
@gtscript.stencil(backend=cartesian_backend)
def calc_unew(
    tdts8: float,
    tdtsdx: float,
    uold: gtscript.Field[dtype],
    cu: gtscript.Field[dtype],
    cv: gtscript.Field[dtype],
    z: gtscript.Field[dtype],
    h: gtscript.Field[dtype],
    unew: gtscript.Field[dtype]
):
    with computation(PARALLEL), interval(...):
        unew = uold + tdts8 * (z + z[0,1,0]) * (cv[0,1,0] + cv + cv[-1,1,0] + cv[-1,0,0]) - tdtsdx * (h - h[-1,0,0])

    #vnew[i,j+1,0] = vold[i,j+1,0] - tdts8 * (z[i+1,j+1,0] + z[i,j+1,0]) * (cu[i+1,j+1,0] + cu[i+1,j,0] + cu[i,j+1,0] + cu[i,j,0]) - tdtsdy * (h[i,j+1,0] - h[i,j,0])    
@gtscript.stencil(backend=cartesian_backend)
def calc_vnew(
    tdts8: float,
    tdtsdy: float,
    vold: gtscript.Field[dtype],
    cu: gtscript.Field[dtype],
    cv: gtscript.Field[dtype],
    z: gtscript.Field[dtype],
    h: gtscript.Field[dtype],
    vnew: gtscript.Field[dtype]
):
    with computation(PARALLEL), interval(...):
        vnew = vold - tdts8 * (z[1,0,0] + z) * (cu[1,0,0] + cu[1,-1,0] + cu + cu[0,-1,0]) - tdtsdy * (h - h[0,-1,0])
        
@gtscript.stencil(backend=cartesian_backend)
def calc_pold(
    p: gtscript.Field[dtype],
    alpha: float,
    pnew: gtscript.Field[dtype],
    pold: gtscript.Field[dtype]
):
    with computation(PARALLEL), interval(...):
        pold = p + alpha * (pnew - 2 * p + pold)

@gtscript.stencil(backend=cartesian_backend)
def calc_uold(
    u: gtscript.Field[dtype],
    alpha: float,
    unew: gtscript.Field[dtype],
    uold: gtscript.Field[dtype]
):
    with computation(PARALLEL), interval(...):
        uold = u + alpha * (unew - 2 * u + uold)

@gtscript.stencil(backend=cartesian_backend)
def calc_vold(
    v: gtscript.Field[dtype],
    alpha: float,
    vnew: gtscript.Field[dtype],
    vold: gtscript.Field[dtype]
):
    with computation(PARALLEL), interval(...):
        vold = v + alpha * (vnew - 2 * v + vold)
        
@gtscript.stencil(backend=cartesian_backend)
def copy_var(
    inp: gtscript.Field[dtype],
    out: gtscript.Field[dtype]
):
    with computation(PARALLEL), interval(...):
        out = inp

def main():
    dt0  = 0.
    dt05 = 0.
    dt1  = 0.
    dt15 = 0.
    dt2  = 0.
    dt25 = 0.
    dt3  = 0.
    dt35 = 0.

    M_LEN = config.M_LEN
    N_LEN = config.N_LEN
    M = config.M
    N = config.N
    ITMAX = config.ITMAX


    unew = np.zeros((M_LEN, N_LEN, 1))
    vnew = np.zeros((M_LEN, N_LEN, 1))
    pnew = np.zeros((M_LEN, N_LEN, 1))
    uold = np.zeros((M_LEN, N_LEN, 1))
    vold = np.zeros((M_LEN, N_LEN, 1))
    pold = np.zeros((M_LEN, N_LEN, 1))

    cu = np.zeros((M_LEN, N_LEN, 1))
    cu_new = np.zeros((M_LEN, N_LEN, 1))

    cv = np.zeros((M_LEN, N_LEN, 1))
    z = np.zeros((M_LEN, N_LEN, 1))
    h = np.zeros((M_LEN, N_LEN, 1))

    u, v, p = initial_conditions.initialize(M, N, config.dx, config.dy, config.a)
    u = u[:,:,np.newaxis]
    v = v[:,:,np.newaxis]
    p = p[:,:,np.newaxis]


    # Save initial conditions
    uold = np.copy(u)
    vold = np.copy(v)
    pold = np.copy(p)

    # Print initial conditions
    if config.L_OUT:
        print(" Number of points in the x direction: ", M)
        print(" Number of points in the y direction: ", N)
        print(" grid spacing in the x direction: ", config.dx)
        print(" grid spacing in the y direction: ", config.dy)
        print(" time step: ", config.dt)
        print(" time filter coefficient: ", config.alpha)
        print(" Initial p:\n", p[:,:,0].diagonal()[:-1])
        print(" Initial u:\n", u[:,:,0].diagonal()[:-1])
        print(" Initial v:\n", v[:,:,0].diagonal()[:-1])

    nx = M
    ny = N
    nz = 1
    gt4py_type = "cartesian"
    #gt4py_type = "next"
    allocator = gtx.itir_python
    # allocator = gtx.gtfn_gpu

    domain = gtx.domain({I:nx+1, J:ny+1, K:nz})

    h_gt = gtx.as_field(domain,h,allocator=allocator)
    z_gt = gtx.as_field(domain,z,allocator=allocator)

    cu_gt = gtx.as_field(domain,cu,allocator=allocator)
    cu_new_gt = gtx.as_field(domain,cu_new,allocator=allocator)

    cv_gt = gtx.as_field(domain,cv,allocator=allocator)
    pnew_gt = gtx.as_field(domain,pnew,allocator=allocator)
    unew_gt = gtx.as_field(domain,unew,allocator=allocator)
    vnew_gt = gtx.as_field(domain,vnew,allocator=allocator)
    uold_gt = gtx.as_field(domain,uold,allocator=allocator)
    vold_gt = gtx.as_field(domain,vold,allocator=allocator)
    pold_gt = gtx.as_field(domain,pold,allocator=allocator)

    u_gt = gtx.as_field(domain,u,allocator=allocator)
    p_gt = gtx.as_field(domain,p,allocator=allocator)
    v_gt = gtx.as_field(domain,v,allocator=allocator)


    t0_start = perf_counter()
    time = 0.0
    tdt = config.dt
    # Main time loop
    for ncycle in range(ITMAX):
        t05_start = perf_counter()
        # p_gt = gtx.as_field(domain,p,allocator=allocator)
        # u_gt = gtx.as_field(domain,u,allocator=allocator)
        # v_gt = gtx.as_field(domain,v,allocator=allocator)
        t05_stop = perf_counter()
        dt05 = dt05 + (t05_stop - t05_start)

        if((ncycle%100==0) & (config.VIS==False)):
            print(f"cycle number{ncycle} and gt4py type {gt4py_type}")
        
        if config.VAL_DEEP and ncycle <= 3:
            u = u_gt.asnumpy()
            v = v_gt.asnumpy()
            p = p_gt.asnumpy()
            utils.validate_uvp(u, v, p, M, N, ncycle, 'init')
            utils.validate_uvp(u_gt.asnumpy(), v_gt.asnumpy(), p_gt.asnumpy(), M, N, ncycle, 'init')
        
        t1_start = perf_counter()
        # Calculate cu, cv, z, and h
        calc_h(p=p_gt, u=u_gt, v=v_gt, h=h_gt, origin=(0,0,0), domain=(nx,ny,nz)) 
        calc_z(fsdx=config.fsdx, fsdy=config.fsdy, u=u_gt, v=v_gt, p=p_gt, z=z_gt, origin=(1,1,0), domain=(nx,ny,nz)) # domain(nx+1,ny+1,nz) gives error why?
        calc_cu(u=u_gt, p=p_gt, cu=cu_gt, origin=(1,0,0), domain=(nx,ny,nz)) # (nx,ny+1,nz)-->works domain(nx+1,ny+1,nz) gives error why? try removing ny+1
        calc_cv(v=v_gt, p=p_gt, cv=cv_gt, origin=(0,1,0), domain=(nx,ny,nz)) #(nx+1,ny,nz)--> works domain(nx+1,ny+1,nz) gives error why?
        t1_stop = perf_counter()
        dt1 = dt1 + (t1_stop - t1_start)

        t15_start = perf_counter()
        # # Periodic Boundary conditions
        #try region


        cu = cu_gt.asnumpy()

        # update_bndry(old=cu_gt,new=cu_new_gt,north=0,south=0,west=1,east=0,domain=(nx,ny,nz))
        # np.testing.assert_allclose(cu_gt, cu_new_gt)
        cu_new= cu_new_gt.asnumpy()
        # print('numpy: old, new:[0,0] ',cu[0,0],cu_new[0,0])
        # print('numpy: old, new:[0,16] ',cu[0,16],cu_new[0,16])
        # print('numpy: old, new:[16,0] ',cu[16,0],cu_new[16,0])
        # print('numpy: old, new:[16,16] ',cu[16,16],cu_new[16,16])
        # exit()
        #            i[0:M]       
        #        0 -----------------
        #        |X X X X X X X X  |
        #        |---------------  |
        # j[0:N] |              |  |
        #        |              |  |
        #        |              |  |
        #        |              |  |
        #        |------------------
        cu[0, :,0] = cu[M, :,0]


         #           i[0:M]       
         #        0 -----------------
         #        |                 |
         #        |--------------- X|
         # j[0:N] |              | X|
         #        |              | X|
         #        |              | X|
         #        |              | X|
         #        |------------------
        cu[1:, N,0] = cu[1:, 0,0]

         #            i[0:M]       
         #        0 -----------------
         #        |                X|
         #        |---------------  |
         # j[0:N] |              |  |
         #        |              |  |
         #        |              |  |
         #        |              |  |
         #        |------------------
        cu[0, N,0] = cu[M, 0,0]

        h = h_gt.asnumpy()
        h[M, :,0] = h[0, :,0]
        h[:, N,0] = h[:, 0,0]
        h[M, N,0] = h[0, 0,0]

        cv = cv_gt.asnumpy()
        cv[M, 1:,0] = cv[0, 1:,0]
        cv[:, 0,0] = cv[:, N,0]
        cv[M, 0,0] = cv[0, N,0]

        z = z_gt.asnumpy()
        z[0, 1:,0] = z[M, 1:,0]
        z[1:, 0,0] = z[1:, N,0]
        z[0, 0,0] = z[M, N,0]

        t15_stop = perf_counter()
        dt15 = dt15 + (t15_stop - t15_start)
            
        if config.VAL_DEEP and ncycle <=1:
            utils.validate_cucvzh(cu, cv, z, h, M, N, ncycle, 't100')
            
        # Calclulate new values of u,v, and p
        tdts8 = tdt / 8.
        tdtsdx = tdt / config.dx
        tdtsdy = tdt / config.dy
        #print(tdts8, tdtsdx, tdtsdy)

        t2_start = perf_counter()
        calc_unew(tdts8=tdts8, tdtsdx=tdtsdx, uold=uold_gt, cu=cu_gt, cv=cv_gt, z=z_gt, h=h_gt, unew=unew_gt, origin=(1,0,0), domain=(nx,ny,nz))
        calc_vnew(tdts8=tdts8, tdtsdy=tdtsdy, vold=vold_gt, cu=cu_gt, cv=cv_gt, z=z_gt, h=h_gt, vnew=vnew_gt, origin=(0,1,0), domain=(nx,ny,nz))
        calc_pnew(tdtsdx=tdtsdx, tdtsdy=tdtsdy, pold=pold_gt, cu=cu_gt, cv=cv_gt, pnew=pnew_gt, origin=(0,0,0), domain=(nx,ny,nz))
        t2_stop = perf_counter()
        dt2 = dt2 + (t2_stop - t2_start)
            
        # for i in range(M):
        #     for j in range(N):
        #         unew[i+1,j,0] = uold[i+1,j,0] + tdts8 * (z[i+1,j+1,0] + z[i+1,j,0]) * (cv[i+1,j+1,0] + cv[i+1,j,0] + cv[i,j+1,0] + cv[i,j,0]) - tdtsdx * (h[i+1,j,0] - h[i,j,0])
        #         vnew[i,j+1,0] = vold[i,j+1,0] - tdts8 * (z[i+1,j+1,0] + z[i,j+1,0]) * (cu[i+1,j+1,0] + cu[i+1,j,0] + cu[i,j+1,0] + cu[i,j,0]) - tdtsdy * (h[i,j+1,0] - h[i,j,0])
        #         pnew[i,j,0] = pold[i,j,0] - tdtsdx * (cu[i+1,j,0] - cu[i,j,0]) - tdtsdy * (cv[i,j+1,0] - cv[i,j,0])
                            
        t25_start = perf_counter()
        # Periodic Boundary conditions

        unew = unew_gt.asnumpy()
        unew[0, :,0] = unew[M, :,0]
        unew[1:, N,0] = unew[1:, 0,0]
        unew[0, N,0] = unew[M, 0,0]

        pnew = pnew_gt.asnumpy()
        pnew[M, :,0] = pnew[0, :,0]
        pnew[:, N,0] = pnew[:, 0,0]
        pnew[M, N,0] = pnew[0, 0,0]

        vnew = vnew_gt.asnumpy()
        vnew[M, 1:,0] = vnew[0, 1:,0]
        vnew[:, 0,0] = vnew[:, N,0]
        vnew[M, 0,0] = vnew[0, N,0]

        t25_stop = perf_counter()
        dt25 = dt25 + (t25_stop - t25_start)
        
        if config.VAL_DEEP and ncycle <= 1:
            utils.validate_uvp(unew, vnew, pnew, M, N, ncycle, 't200')
        
        time = time + config.dt

        if(ncycle > 0):
            t3_start = perf_counter()
            calc_pold(p=p_gt, alpha=config.alpha, pnew=pnew_gt, pold=pold_gt, origin=(0,0,0), domain=(nx+1,ny+1,nz))
            calc_uold(u=u_gt, alpha=config.alpha, unew=unew_gt, uold=uold_gt, origin=(0,0,0), domain=(nx+1,ny+1,nz))
            calc_vold(v=v_gt, alpha=config.alpha, vnew=vnew_gt, vold=vold_gt, origin=(0,0,0), domain=(nx+1,ny+1,nz))

            # I don't think we need these .asnumpy() calls 
            #pold = pold_gt.asnumpy()
            #uold = uold_gt.asnumpy()
            #vold = vold_gt.asnumpy()
            
            copy_var(unew_gt, u_gt, origin=(0,0,0), domain=(nx+1,ny+1,nz))
            copy_var(vnew_gt, v_gt, origin=(0,0,0), domain=(nx+1,ny+1,nz))
            copy_var(pnew_gt, p_gt, origin=(0,0,0), domain=(nx+1,ny+1,nz))
            t3_stop = perf_counter()
            dt3 = dt3 + (t3_stop - t3_start)

            t35_start = perf_counter()
            # u = u_gt.asnumpy()
            # v = v_gt.asnumpy()
            # p = p_gt.asnumpy()
            t35_stop = perf_counter()
            dt35 = dt35 + (t35_stop - t35_start)


        else:
            tdt = tdt+tdt

            copy_var(u_gt,uold_gt, origin=(0,0,0), domain=(nx+1,ny+1,nz))
            copy_var(v_gt,vold_gt, origin=(0,0,0), domain=(nx+1,ny+1,nz))
            copy_var(p_gt,pold_gt, origin=(0,0,0), domain=(nx+1,ny+1,nz))
            # uold = np.copy(u[...])
            # vold = np.copy(v[...])
            # pold = np.copy(p[...])

            copy_var(unew_gt,u_gt, origin=(0,0,0), domain=(nx+1,ny+1,nz))
            copy_var(vnew_gt,v_gt, origin=(0,0,0), domain=(nx+1,ny+1,nz))
            copy_var(pnew_gt,p_gt, origin=(0,0,0), domain=(nx+1,ny+1,nz))
            # u = np.copy(unew[...])
            # v = np.copy(vnew[...])
            # p = np.copy(pnew[...])

        if((config.VIS == True) & (ncycle%config.VIS_DT==0)):
            u = u_gt.asnumpy() 
            v = v_gt.asnumpy() 
            p = p_gt.asnumpy() 
            utils.live_plot3(u, v, p, "ncycle: " + str(ncycle))
            
    t0_stop = perf_counter()
    dt0 = dt0 + (t0_stop - t0_start)
    # Print initial conditions
    if config.L_OUT:
            print("cycle number ", ITMAX)
            print(" diagonal elements of p:\n", pnew[:,:,0].diagonal()[:-1])
            print(" diagonal elements of u:\n", unew[:,:,0].diagonal()[:-1])
            print(" diagonal elements of v:\n", vnew[:,:,0].diagonal()[:-1])
    print("total: ",dt0)
    print("t050: ",dt05)
    print("t100: ",dt1)
    print("t150: ",dt15)
    print("t200: ",dt2)
    print("t250: ",dt25)
    print("t300: ",dt3)
    print("t350: ",dt35)

    if config.VAL:
        utils.final_validation(u, v, p, ITMAX=ITMAX, M=M, N=N)

if __name__ == "__main__":
    main()
