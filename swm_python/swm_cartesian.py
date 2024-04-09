import numpy as np
import numpy as np
import gt4py.next as gtx
import gt4py.cartesian.gtscript as gtscript
import gt4py.storage as gts
from time import perf_counter
# import cupy
# import gt4py
import initial_conditions
import utils
import config

I = gtx.Dimension("I")
J = gtx.Dimension("J")
K = gtx.Dimension("K", kind = gtx.DimensionKind.VERTICAL)

dtype = np.float64

cartesian_backend = config.backend
allocator = gtx.gtfn_cpu
if cartesian_backend in ("gt:gpu", "cuda", "dace:gpu"):
    allocator = gtx.gtfn_gpu

print(f"Using {cartesian_backend} backend with {allocator.__name__} allocator.")

@gtscript.stencil(backend=cartesian_backend)
def calc_cucvzh(u: gtscript.Field[dtype],
                v: gtscript.Field[dtype],
                p: gtscript.Field[dtype],
                cu: gtscript.Field[dtype],
                cv: gtscript.Field[dtype],
                z: gtscript.Field[dtype],
                h: gtscript.Field[dtype],
                fsdx:float,
                fsdy: float):

    with computation(PARALLEL), interval(...):
        cu = .5 * (p[1,0,0] + p) * u
        cv = .5 * (p[0,1,0] + p) * v
        z = (config.fsdx * (v[1,0,0] - v) - config.fsdy * (u[0,1,0] - u)) / (p[1,1,0] + p[0,1,0] + p + p[1,0,0])
        h = p + 0.25 * (u[-1,0,0] * u[-1,0,0] + u * u + v[0,-1,0] * v[0,-1,0] + v * v)

@gtscript.stencil(backend=cartesian_backend)
def calc_uvp(
    tdts8: float,
    tdtsdx: float,
    tdtsdy: float,
    uold: gtscript.Field[dtype],
    vold: gtscript.Field[dtype],
    pold: gtscript.Field[dtype],
    cu: gtscript.Field[dtype],
    cv: gtscript.Field[dtype],
    z: gtscript.Field[dtype],
    h: gtscript.Field[dtype],
    unew: gtscript.Field[dtype],
    vnew: gtscript.Field[dtype],
    pnew: gtscript.Field[dtype]):

    with computation(PARALLEL), interval(...):
        unew = uold + tdts8 * (z + z[0,-1,0]) * (cv[1,0,0] + cv + cv[0,-1,0] + cv[1,-1,0]) - tdtsdx * (h[1,0,0] - h)
        vnew = vold - tdts8 * (z + z[-1,0,0]) * (cu[-1,0,0] + cu[-1,1,0] + cu + cu[0,1,0]) - tdtsdy * (h[0,1,0] - h)
        pnew = pold - tdtsdx * (cu- cu[-1,0,0]) - tdtsdy * (cv-cv[0,-1,0])

@gtscript.stencil(backend=cartesian_backend)
def calc_uvp_old(
    alpha: float,
    v: gtscript.Field[dtype],
    vnew: gtscript.Field[dtype],
    vold: gtscript.Field[dtype],
    u: gtscript.Field[dtype],
    unew: gtscript.Field[dtype],
    uold: gtscript.Field[dtype],
    p: gtscript.Field[dtype],
    pnew: gtscript.Field[dtype],
    pold: gtscript.Field[dtype]):

    with computation(PARALLEL), interval(...):
        uold = u + alpha * (unew - 2 * u + uold)
        vold = v + alpha * (vnew - 2 * v + vold)
        pold = p + alpha * (pnew - 2 * p + pold)


@gtscript.stencil(backend=cartesian_backend)
def copy_3var(inp0: gtscript.Field[dtype],
              inp1: gtscript.Field[dtype],
              inp2: gtscript.Field[dtype],
              out0: gtscript.Field[dtype],
              out1: gtscript.Field[dtype],
              out2: gtscript.Field[dtype]):

    with computation(PARALLEL), interval(...):
        out0 = inp0
        out1 = inp1
        out2 = inp2


def main():
    dt0 = 0.
    dt1 = 0.
    dt15 = 0.
    dt2 = 0.
    dt25 = 0.
    dt3 = 0.

    M_LEN = config.M_LEN
    N_LEN = config.N_LEN
    M = config.M
    N = config.N
    
    _u, _v, _p = initial_conditions.initialize(M, N, config.dx, config.dy, config.a)
    _u = _u[:,:,np.newaxis]
    _v = _v[:,:,np.newaxis]
    _p = _p[:,:,np.newaxis]

    domain = gtx.domain({I:M+1, J:N+1, K:1})

    # if cartesian_backend in ("gpu","gt:cpu_ifirst","cuda", "dace:gpu"):
    if cartesian_backend in ("cuda", "dace:gpu"):
        h_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        z_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        cu_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        cv_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        pnew_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        unew_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        vnew_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        uold_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        vold_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        pold_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        u_gt = gtx.as_field(domain,_u,allocator=allocator)
        v_gt = gtx.as_field(domain,_v,allocator=allocator)
        p_gt = gtx.as_field(domain,_p,allocator=allocator)
    else:
        shape = (M+1,N+1,1)
        h_gt  = gts.empty(dtype=dtype,backend=cartesian_backend,shape=shape)
        z_gt  = gts.empty(dtype=dtype,backend=cartesian_backend,shape=shape)
        cu_gt = gts.empty(dtype=dtype,backend=cartesian_backend,shape=shape)
        cv_gt = gts.empty(dtype=dtype,backend=cartesian_backend,shape=shape)
        pnew_gt = gts.empty(dtype=dtype,backend=cartesian_backend,shape=shape)
        unew_gt = gts.empty(dtype=dtype,backend=cartesian_backend,shape=shape)
        vnew_gt = gts.empty(dtype=dtype,backend=cartesian_backend,shape=shape)
        uold_gt = gts.empty(dtype=dtype,backend=cartesian_backend,shape=shape)
        vold_gt = gts.empty(dtype=dtype,backend=cartesian_backend,shape=shape)
        pold_gt = gts.empty(dtype=dtype,backend=cartesian_backend,shape=shape)
        u_gt    = gts.from_array(_u,dtype=dtype,backend=cartesian_backend)
        v_gt    = gts.from_array(_v,dtype=dtype,backend=cartesian_backend)
        p_gt    = gts.from_array(_p,dtype=dtype,backend=cartesian_backend)
        #
        #  Create different fields.  These fields are:
        #      (1) Much slower than the regular fields 
        #      (2) they support the asnumpy() method needed for validation
        #
        utmp_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        vtmp_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        ptmp_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        cutmp_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        cvtmp_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        htmp_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)
        ztmp_gt = gtx.empty(domain,dtype=dtype,allocator=allocator)

    # Save initial conditions
    uold_gt[...] = u_gt[...]
    vold_gt[...] = v_gt[...]
    pold_gt[...] = p_gt[...]

    # Print initial conditions
    if config.L_OUT:
        print(" Number of points in the x direction: ", M)
        print(" Number of points in the y direction: ", N)
        print(" grid spacing in the x direction: ", config.dx)
        print(" grid spacing in the y direction: ", config.dy)
        print(" time step: ", config.dt)
        print(" time filter coefficient: ", config.alpha)
        #
        # asnumpy() method is not working with from_array allocator
        #
        ptmp_gt[...] = p_gt[...]
        utmp_gt[...] = u_gt[...]
        vtmp_gt[...] = v_gt[...]
        print(" Initial p:\n", ptmp_gt.asnumpy()[:,:,0].diagonal()[:-1])
        print(" Initial u:\n", utmp_gt.asnumpy()[:,:,0].diagonal()[:-1])
        print(" Initial v:\n", vtmp_gt.asnumpy()[:,:,0].diagonal()[:-1])

    t0_start = perf_counter()
    time = 0.0
    tdt = config.dt

    u_origin=(1,0,0)
    v_origin=(0,1,0)
    p_origin=(0,0,0)
    z_origin=(1,1,0)
    # Main time loop
    for ncycle in range(config.ITMAX):

        if((ncycle%100==0) & (config.VIS==False)):
            print(f"cycle number{ncycle} and gt4py type cartesian")

        if config.VAL_DEEP and ncycle <= 3:
            utmp_gt[...]=u_gt[...]
            vtmp_gt[...]=v_gt[...]
            ptmp_gt[...]=p_gt[...]
            utils.validate_uvp(utmp_gt.asnumpy(), vtmp_gt.asnumpy(), ptmp_gt.asnumpy(), M, N, ncycle, 'init')

        t1_start = perf_counter()

        calc_cucvzh(
            u=u_gt,
            v=v_gt,
            p=p_gt,
            cu=cu_gt,
            cv=cv_gt,
            z=z_gt,
            h=h_gt,
            fsdx=config.fsdx,
            fsdy=config.fsdy,
            origin={"u":u_origin, "v":v_origin, "p":p_origin, "z":z_origin, "h":p_origin, "cu":u_origin, "cv":v_origin},
            domain=(M, N, 1),
        )

        t1_stop = perf_counter()
        t15_start = perf_counter()
        dt1 = dt1 + (t1_stop - t1_start)

        t15_start = perf_counter()
        # # Periodic Boundary conditions
        # try region
        cu_gt[0, :,0] = cu_gt[M, :,0]
        # update_boundary(cu_gt, cu_gt, M, N)

        h_gt[M, :,0] = h_gt[0, :,0]
        cv_gt[M, 1:,0] = cv_gt[0, 1:,0]
        z_gt[0, 1:,0] = z_gt[M, 1:,0]

        cv_gt[:, 0,0] = cv_gt[:, N,0]
        h_gt[:, N,0] = h_gt[:, 0,0]
        cu_gt[1:, N,0] = cu_gt[1:, 0,0]
        z_gt[1:, 0,0] = z_gt[1:, N,0]

        cu_gt[0, N,0] = cu_gt[M, 0,0]
        cv_gt[M, 0,0] = cv_gt[0, N,0]
        z_gt[0, 0,0] = z_gt[M, N,0]
        h_gt[M, N,0] = h_gt[0, 0,0]

        t15_stop = perf_counter()
        dt15 = dt15 + (t15_stop - t15_start)

        if config.VAL_DEEP and ncycle <=1:
            cutmp_gt[...] = cu_gt[...]
            cvtmp_gt[...] = cv_gt[...]
            ztmp_gt[...]  = z_gt[...]
            htmp_gt[...]  = h_gt[...]
            utils.validate_cucvzh(cutmp_gt.asnumpy(), cvtmp_gt.asnumpy(), ztmp_gt.asnumpy(), htmp_gt.asnumpy(), M, N, ncycle, 't100')

        # Calclulate new values of u,v, and p
        tdts8 = tdt / 8.
        tdtsdx = tdt / config.dx
        tdtsdy = tdt / config.dy
        # print(tdts8, tdtsdx, tdtsdy)

        t2_start = perf_counter()
        
        calc_uvp(
            tdts8=tdts8,
            tdtsdx=tdtsdx,
            tdtsdy=tdtsdy,
            uold=uold_gt,
            vold=vold_gt,
            pold=pold_gt,
            cu=cu_gt,
            cv=cv_gt,
            z=z_gt,
            h=h_gt,
            unew=unew_gt,
            vnew=vnew_gt,
            pnew=pnew_gt,
            origin={
                "uold": u_origin,
                "vold": v_origin,
                "pold": p_origin,
                "cu": u_origin,
                "cv": v_origin,
                "z": z_origin,
                "h": p_origin,
                "unew": u_origin,
                "vnew": v_origin,
                "pnew": p_origin,
            },
            domain=(M, N, 1),
        )

        t2_stop = perf_counter()
        t25_start = perf_counter()
        dt2 = dt2 + (t2_stop - t2_start)

        # Periodic Boundary conditions
        unew_gt[0, :,0] = unew_gt[M, :,0]
        pnew_gt[M, :,0] = pnew_gt[0, :,0]
        vnew_gt[M, 1:,0] = vnew_gt[0, 1:,0]
        unew_gt[1:, N,0] = unew_gt[1:, 0,0]
        vnew_gt[:, 0,0] = vnew_gt[:, N,0]
        pnew_gt[:, N,0] = pnew_gt[:, 0,0]

        unew_gt[0, N,0] = unew_gt[M, 0,0]
        vnew_gt[M, 0,0] = vnew_gt[0, N,0]
        pnew_gt[M, N,0] = pnew_gt[0, 0,0]

        t25_stop = perf_counter()
        dt25 = dt25 + (t25_stop - t25_start)

        
        if config.VAL_DEEP and ncycle <= 1:
            utmp_gt[...]=unew_gt[...]
            vtmp_gt[...]=vnew_gt[...]
            ptmp_gt[...]=pnew_gt[...]
            utils.validate_uvp(utmp_gt.asnumpy(), vtmp_gt.asnumpy(), ptmp_gt.asnumpy(), M, N, ncycle, 't200')

        time = time + config.dt

        if(ncycle > 0):
            t3_start = perf_counter()
            calc_uvp_old(alpha=config.alpha, v=v_gt, vnew=vnew_gt, vold=vold_gt, u=u_gt, unew=unew_gt, uold=uold_gt, p=p_gt, pnew=pnew_gt, pold=pold_gt, domain=(M+1, N+1, 1))
            copy_3var(unew_gt, vnew_gt, pnew_gt, u_gt, v_gt, p_gt, origin=(0,0,0), domain=(M+1,N+1,1))

            t3_stop = perf_counter()
            dt3 = dt3 + (t3_stop - t3_start)

        else:
            tdt = tdt+tdt

            uold_gt[...] = u_gt[...]
            vold_gt[...] = v_gt[...]
            pold_gt[...] = p_gt[...]
            u_gt[...] = unew_gt[...]
            v_gt[...] = vnew_gt[...]
            p_gt[...] = pnew_gt[...]

        if((config.VIS == True) & (ncycle%config.VIS_DT==0)):
            utmp_gt[...]=u_gt[...]
            vtmp_gt[...]=v_gt[...]
            ptmp_gt[...]=p_gt[...]
            utils.live_plot3(utmp_gt.asnumpy(), vtmp_gt.asnumpy(), ptmp_gt.asnumpy(), "ncycle: " + str(ncycle))

    t0_stop = perf_counter()
    dt0 = dt0 + (t0_stop - t0_start)
    # Print initial conditions
    if config.L_OUT:
        utmp_gt[...]=u_gt[...]
        vtmp_gt[...]=v_gt[...]
        ptmp_gt[...]=p_gt[...]
        print("cycle number ", config.ITMAX)
        print(" diagonal elements of p:\n", ptmp_gt.asnumpy()[:,:,0].diagonal()[:-1])
        print(" diagonal elements of u:\n", utmp_gt.asnumpy()[:,:,0].diagonal()[:-1])
        print(" diagonal elements of v:\n", vtmp_gt.asnumpy()[:,:,0].diagonal()[:-1])
    print("total: ",dt0)
    print("t100: ",dt1)
    print("t150: ",dt15)
    print("t200: ",dt2)
    print("t250: ",dt25)
    print("t300: ",dt3)

    if config.VAL:
        utmp_gt[...]=u_gt[...]
        vtmp_gt[...]=v_gt[...]
        ptmp_gt[...]=p_gt[...]
        utils.final_validation(utmp_gt.asnumpy(), vtmp_gt.asnumpy(), ptmp_gt.asnumpy(), ITMAX=config.ITMAX, M=M, N=N)

if __name__ == "__main__":
    main()
