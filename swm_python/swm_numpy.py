import numpy as np
import argparse
from time import perf_counter
import utils
import initial_conditions
import config

def main():

    dt0=0.
    dt1=0.
    dt2=0.
    dt3=0.

    M_LEN = config.M_LEN
    N_LEN = config.N_LEN
    M = config.M
    N = config.N
    ITMAX = config.ITMAX

    # Model Variables
    unew = np.zeros((M_LEN, N_LEN))
    vnew = np.zeros((M_LEN, N_LEN))
    pnew = np.zeros((M_LEN, N_LEN))
    cu = np.zeros((M_LEN, N_LEN))
    cv = np.zeros((M_LEN, N_LEN))
    z = np.zeros((M_LEN, N_LEN))
    h = np.zeros((M_LEN, N_LEN))
    
    u, v, p = initial_conditions.initialize(M, N, config.dx, config.dy, config.a)
    
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
        print(" Initial p:\n", p.diagonal()[:-1])
        print(" Initial u:\n", u.diagonal()[:-1])
        print(" Initial v:\n", v.diagonal()[:-1])
    
    t0_start = perf_counter()
    time = 0.0 

    tdt = config.dt
    # Main time loop
    for ncycle in range(ITMAX):
        if config.VAL_DEEP and ncycle <= 3:
            utils.validate_uvp(u, v, p, M, N, ncycle, 'init')

        if ncycle % 100 == 0:
            print("cycle number ", ncycle)
                
        t1_start = perf_counter()
        
        cu[1:,:-1] = .5 * (p[1:,:-1] + p[:-1,:-1]) * u[1:,:-1]
        cv[:-1,1:] = .5 * (p[:-1,1:] + p[:-1,:-1]) * v[:-1,1:]
        z[1:,1:] = (config.fsdx * (v[1:,1:] - v[:-1,1:]) -
                    config.fsdy * (u[1:,1:] - u[1:,:-1] )
                    ) / (p[:-1,:-1] + p[1:,:-1] + p[1:,1:] + p[:-1,1:])
        h[:-1,:-1] = p[:-1,:-1] + 0.25 * (u[1:,:-1] * u[1:, :-1] + u[:-1,:-1] * u[:-1,:-1] +
                                          v[:-1,1:] * v[:-1,1:] + v[:-1,:-1] * v[:-1,:-1])
        t1_stop = perf_counter()
        dt1 = dt1 + (t1_stop - t1_start)
        
        # # Periodic Boundary conditions
        cu[0, :] = cu[M, :]
        h[M, :] = h[0, :]
        cv[M, 1:] = cv[0, 1:]
        z[0, 1:] = z[M, 1:]
        
        cv[:, 0] = cv[:, N]
        h[:, N] = h[:, 0]
        cu[1:, N] = cu[1:, 0]
        z[1:, 0] = z[1:, N]
            
        cu[0, N] = cu[M, 0]
        cv[M, 0] = cv[0, N]
        z[0, 0] = z[M, N]
        h[M, N] = h[0, 0]
        
        if config.VAL_DEEP and ncycle <=1:
            utils.validate_cucvzh(cu, cv, z, h, M, N, ncycle, 't100')

        # Calclulate new values of u,v, and p
        tdts8 = tdt / np.float64("8.")
        tdtsdx = tdt / config.dx
        tdtsdy = tdt / config.dy
        
        t2_start = perf_counter()
        
        unew[1:,:-1] = uold[1:,:-1] + tdts8 * (z[1:,1:] + z[1:,:-1]) * (cv[1:,1:] + cv[1:,:-1] + cv[:-1,1:] + cv[:-1,:-1]) - tdtsdx * (h[1:,:-1] - h[:-1,:-1])
        vnew[:-1,1:] = vold[:-1,1:] - tdts8 * (z[1:,1:] + z[:-1,1:]) * (cu[1:,1:] + cu[1:,:-1] + cu[:-1,1:] + cu[:-1,:-1]) - tdtsdy * (h[:-1,1:] - h[:-1,:-1])
        pnew[:-1,:-1] = pold[:-1,:-1] - tdtsdx * (cu[1:,:-1] - cu[:-1,:-1]) - tdtsdy * (cv[:-1,1:] - cv[:-1,:-1])
        t2_stop = perf_counter()
        dt2 = dt2 + (t2_stop - t2_start)
                
        # Periodic Boundary conditions
        unew[0, :] = unew[M, :]
        pnew[M, :] = pnew[0, :]
        vnew[M, 1:] = vnew[0, 1:]
        unew[1:, N] = unew[1:, 0]
        vnew[:, 0] = vnew[:, N]
        pnew[:, N] = pnew[:, 0]
        
        unew[0, N] = unew[M, 0]
        vnew[M, 0] = vnew[0, N]
        pnew[M, N] = pnew[0, 0]

        if config.VAL_DEEP and ncycle <= 1:
            utils.validate_uvp(unew, vnew, pnew, M, N, ncycle, 't200')

        time = time + config.dt

        if(ncycle > 0):
            t3_start = perf_counter()
            
            uold[...] = u + config.alpha * (unew - 2. * u + uold)
            vold[...] = v + config.alpha * (vnew - 2. * v + vold)
            pold[...] = p + config.alpha * (pnew - 2. * p + pold)

            u[...] = unew
            v[...] = vnew
            p[...] = pnew
            t3_stop = perf_counter()
            dt3 = dt3 + (t3_stop - t3_start)
        else:
            tdt = tdt + tdt
            uold = np.copy(u)
            vold = np.copy(v)
            pold = np.copy(p)
            u = np.copy(unew)
            v = np.copy(vnew)
            p = np.copy(pnew)
        
        # Print initial conditions
    t0_stop = perf_counter()
    dt0 = dt0 + (t0_stop - t0_start)
    if config.L_OUT:
        print("cycle number ", ITMAX)
        print(" diagonal elements of p:\n", pnew.diagonal()[:-1])
        print(" diagonal elements of u:\n", unew.diagonal()[:-1])
        print(" diagonal elements of v:\n", vnew.diagonal()[:-1])
    print("total: ",dt0)
    print("t100: ",dt1)
    print("t200: ",dt2)
    print("t300: ",dt3)


    if config.VAL:
        utils.final_validation(u, v, p, ITMAX=ITMAX, M=M, N=N)
        

if __name__ == "__main__":
    main()
