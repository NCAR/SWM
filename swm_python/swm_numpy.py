import numpy as np
import argparse
from time import perf_counter
from IPython.display import clear_output
from matplotlib import pyplot as plt

def read_arrays(u_file, v_file, p_file):
    u_v = np.loadtxt(u_file,)
    v_v = np.loadtxt(v_file)
    p_v = np.loadtxt(p_file)
    return u_v, v_v, p_v

def live_plot_val(fu, fv, fp, title=''):
    mxu = fu.max()
    mxv = fv.max()
    mxp = fp.max()
    clear_output(wait=True)
    fig, (ax1, ax2, ax3) = plt.subplots(figsize=(13, 3), ncols=3)

    pos1 = ax1.imshow(fp, cmap='Blues', vmin=-mxp, vmax=mxp,interpolation='none')
    ax1.set_title('p')
    plt.colorbar(pos1,ax=ax1)
    pos2 = ax2.imshow(fu, cmap='Reds', vmin=-mxu, vmax=mxu,interpolation='none')
    ax2.set_title('u')
    plt.colorbar(pos2,ax=ax2)
    pos3 = ax3.imshow(fv, cmap='Greens',vmin=-mxv, vmax=mxv,interpolation='none')
    ax3.set_title('v')
    plt.colorbar(pos3, ax=ax3)

    fig.suptitle(title)
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Shallow Water Model")
    parser.add_argument('--M', type=int, default=64, help='Number of points in the x direction')
    parser.add_argument('--N', type=int, default=64, help='Number of points in the y direction')
    parser.add_argument('--L_OUT', type=bool, default=True, help='a boolean for L_OUT')

    
    args = parser.parse_args()
    print(args.M, args.N)
    
    # Initialize model parameters
    M = args.M
    N = args.N
    #M = 64
    #N = 64
    M_LEN = M + 1
    N_LEN = N + 1
    L_OUT = args.L_OUT
    VAL=True
    ITMAX = 4000
    dt = 90.
    tdt = dt
    dx = 100000.
    dy = 100000.
    fsdx = 4. / (dx)
    fsdy = 4. / (dy)
    a = 1000000.
    alpha = 0.001
    el = N * dx
    pi = 4. * np.arctan(1.)
    tpi = 2. * pi
    d_i = tpi / M
    d_j = tpi / N
    dtheta = tpi / M
    pcf = (pi * pi * a * a) / (el * el)
    SIZE = M_LEN * N_LEN
    dt0=0.
    dt1=0.
    dt2=0.
    dt3=0.

    # Model Variables
    u = np.zeros((M_LEN, N_LEN))
    v = np.zeros((M_LEN, N_LEN))
    p = np.zeros((M_LEN, N_LEN))
    unew = np.zeros((M_LEN, N_LEN))
    vnew = np.zeros((M_LEN, N_LEN))
    pnew = np.zeros((M_LEN, N_LEN))
    uold = np.zeros((M_LEN, N_LEN))
    vold = np.zeros((M_LEN, N_LEN))
    pold = np.zeros((M_LEN, N_LEN))
    cu = np.zeros((M_LEN, N_LEN))
    cv = np.zeros((M_LEN, N_LEN))
    z = np.zeros((M_LEN, N_LEN))
    h = np.zeros((M_LEN, N_LEN))
    psi = np.zeros((M_LEN, N_LEN))
    
    # Initial values of the stream function and p
    # for i in range(M + 1):
    #    for j in range(N + 1):
    #        psi[i, j] = a * np.sin((i + .5) * d_i) * np.sin((j + .5) * d_j)
    #        p[i, j] = pcf * (np.cos(2. * (i) * d_i) + np.cos(2. * (j) * d_j)) + 50000.
    psi[...] = a * np.sin((np.arange(0,M_LEN)[:, np.newaxis]+0.5) * d_i) * np.sin((np.arange(0,N_LEN)+ .5) * d_j)
    p[...]   = pcf * (np.cos(2. * np.arange(0,M_LEN)[:, np.newaxis] * d_i) + np.cos(2. * np.arange(0,N_LEN) * d_j)) + 50000.

            
    # Calculate initial u and v
    #for i in range(M):
    #    for j in range(N):
    #        u[i+1,j] = -(psi[i+1,j+1] - psi[i+1,j]) / dy
    #        v[i,j+1] = (psi[i+1,j+1] - psi[i,j+1]) / dx
    u[1:,:-1] = -(psi[1:,1:] - psi[1:,:-1]) / dy
    v[:-1,1:] = (psi[1:,1:] - psi[:-1,1:]) / dx
            
    # # Periodic Boundary conditions
    u[0, :] = u[M, :]
    v[M, 1:] = v[0, 1:]
    u[1:, N] = u[1:, 0]
    v[:, 0] = v[:, N]

    u[0, N] = u[M, 0]
    v[M, 0] = v[0, N]
    
    # Save initial conditions
    uold = np.copy(u)
    vold = np.copy(v)
    pold = np.copy(p)
    
    # Print initial conditions
    if L_OUT:
        print(" Number of points in the x direction: ", M)
        print(" Number of points in the y direction: ", N)
        print(" grid spacing in the x direction: ", dx)
        print(" grid spacing in the y direction: ", dy)
        print(" time step: ", dt)
        print(" time filter coefficient: ", alpha)
        print(" Initial p:\n", p.diagonal()[:-1])
        print(" Initial u:\n", u.diagonal()[:-1])
        print(" Initial v:\n", v.diagonal()[:-1])
    
    t0_start = perf_counter()
    time = 0.0 
    # Main time loop
    for ncycle in range(ITMAX):
        if ncycle % 100 == 0:
            print("cycle number ", ncycle)
        # Calculate cu, cv, z, and h
        # for i in range(1, M):
        #     for j in range(N):
        #         cu[i, j] = .5 * (p[i, j] + p[i - 1, j]) * u[i, j]
        
        # for i in range(M):
        #     for j in range(1, N):
        #         cv[i, j] = .5 * (p[i, j] + p[i, j - 1]) * v[i, j]
                
        # for i in range(1, M):
        #     for j in range(1, N):
        #         z[i, j] = (fsdx * (v[i, j] - v[i - 1, j]) - 
        #                    fsdy * (u[i, j] - u[i, j - 1])
        #                    ) / (p[i - 1, j - 1] + p[i, j - 1] + p[i, j] + p[i - 1, j])
                
        # for i in range(M):
        #     for j in range(N):
        #         h[i, j] = p[i, j] + 0.25 * (u[i + 1, j] * u[i + 1, j] + u[i, j] * u[i, j] + 
        #                                     v[i, j + 1] * v[i, j + 1] + v[i, j] * v[i, j])
                
        t1_start = perf_counter()
        # loop:  i+1 --> 1:
        #        i   --> :-1
        #        j+1 --> 1:
        #        j   --> :-1
        #for i in range(M):
        #    for j in range(N):
        #        cu[i + 1, j] = .5 * (p[i + 1, j] + p[i, j]) * u[i + 1, j]
        #        cv[i, j + 1] = .5 * (p[i, j + 1] + p[i, j]) * v[i, j + 1]
        #        z[i + 1, j + 1] = (fsdx * (v[i + 1, j + 1] - v[i, j + 1]) -
        #                        fsdy * (u[i + 1, j + 1] - u[i+1, j] )
        #                        ) / (p[i, j] + p[i + 1, j] + p[i + 1, j + 1] + p[i, j + 1])
        #        h[i, j] = p[i, j] + 0.25 * (u[i + 1, j] * u[i + 1, j] + u[i, j] * u[i, j] +
        #                                v[i, j + 1] * v[i, j + 1] + v[i, j] * v[i, j])
        cu[1:,:-1] = .5 * (p[1:,:-1] + p[:-1,:-1]) * u[1:,:-1]
        cv[:-1,1:] = .5 * (p[:-1,1:] + p[:-1,:-1]) * v[:-1,1:]
        z[1:,1:] = (fsdx * (v[1:,1:] - v[:-1,1:]) -
                    fsdy * (u[1:,1:] - u[1:,:-1] )
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
        z[1:, N] = z[1:, 0]
            
        cu[0, N] = cu[M, 0]
        cv[M, 0] = cv[0, N]
        z[0, 0] = z[M, N]
        h[M, N] = h[0, 0]

        # Calclulate new values of u,v, and p
        tdts8 = tdt / 8.
        tdtsdx = tdt / dx
        tdtsdy = tdt / dy
        
        t2_start = perf_counter()
        # loop:  i+1 --> 1:
        #        i   --> :-1
        #        j+1 --> 1:
        #        j   --> :-1
        #for i in range(M):
        #    for j in range(N):
        #        unew[i+1,j] = uold[i+1,j] + tdts8 * (z[i+1,j+1] + z[i+1,j]) * (cv[i+1,j+1] + cv[i+1,j] + cv[i,j+1] + cv[i,j]) - tdtsdx * (h[i+1,j] - h[i,j])
        #        vnew[i,j+1] = vold[i,j+1] - tdts8 * (z[i+1,j+1] + z[i,j+1]) * (cu[i+1,j+1] + cu[i+1,j] + cu[i,j+1] + cu[i,j]) - tdtsdy * (h[i,j+1] - h[i,j])
        #        pnew[i,j] = pold[i,j] - tdtsdx * (cu[i+1,j] - cu[i,j]) - tdtsdy * (cv[i,j+1] - cv[i,j])
        unew[1:,:-1] = uold[1:,:-1] + tdts8 * (z[1:,1:] + z[1:,:-1]) * (cv[1:,1:] + cv[1:,:-1] + cv[:-1,1:] + cv[:-1,:-1]) - tdtsdx * (h[1:,:-1] - h[:-1,:-1])
        vnew[:-1,1:] = vold[:-1,1:] - tdts8 * (z[1:,1:] + z[:-1,1:]) * (cu[1:,1:] + cu[1:,:-1] + cu[:-1,1:] + cu[:-1,:-1]) - tdtsdy * (h[:-1,1:] - h[:-1,:-1])
        pnew[:-1,:-1] = pold[:-1,:-1] - tdtsdx * (cu[1:,:-1] - cu[:-1,:-1]) - tdtsdy * (cv[:-1,1:] - cv[:-1,:-1])
        t2_stop = perf_counter()
        dt2 = dt2 + (t2_stop - t2_start)
        # for i in range(M):
        #     for j in range(N):
        #         unew[i + 1, j] = (uold[i + 1, j] + tdts8 * (z[i + 1, j + 1] + z[i + 1, j]) *
        #                         (cv[i + 1, j + 1] + cv[i + 1, j] + cv[i, j + 1] + cv[i, j]) -
        #                         tdtsdx * (h[i + 1, j] - h[i, j])
        #                         )
        # # for i in range(1, M):
        # #     for j in range(N):
        # #         unew[i, j] = (uold[i, j] + tdts8 * (z[i, j + 1] + z[i, j]) * 
        # #                       (cv[i, j + 1] + cv[i, j] + cv[i - 1, j + 1] + cv[i - 1, j]) - 
        # #                       tdtsdx * (h[i, j] - h[i - 1, j])
        # #                       )
        # for i in range(M):
        #     for j in range(N):        
        #         vnew[i, j + 1] = (vold[i, j + 1] - tdts8 * (z[i + 1, j + 1] + z[i, j + 1]) *
        #                         (cu[i + 1, j + 1] + cu[i + 1, j] + cu[i, j + 1] + cu[i, j]) -
        #                         tdtsdy * (h[i, j + 1] - h[i, j])
        #                         )
        # # for i in range(M):
        # #     for j in range(1, N):        
        # #         vnew[i, j] = (vold[i, j] - tdts8 * (z[i + 1, j] + z[i, j]) * 
        # #                       (cu[i + 1, j] + cu[i + 1, j - 1] + cu[i, j] + cu[i, j - 1]) - 
        # #                       tdtsdy * (h[i, j] - h[i, j - 1])
        # #                       )
        # for i in range(M):
        #     for j in range(N):
        #         pnew[i, j] = (pold[i, j] - tdtsdx * (cu[i + 1, j] - cu[i, j]) -
        #                         tdtsdy * (cv[i, j + 1] - cv[i, j])
        #                         )
                
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
        
        time = time + dt

        if(ncycle > 0):
            t3_start = perf_counter()
            #for i in range(M_LEN):
            #    for j in range(N_LEN):
            #        uold[i,j] = u[i,j] + alpha * (unew[i,j] - 2 * u[i,j] + uold[i,j])
            #        vold[i,j] = v[i,j] + alpha * (vnew[i,j] - 2 * v[i,j] + vold[i,j])
            #        pold[i,j] = p[i,j] + alpha * (pnew[i,j] - 2 * p[i,j] + pold[i,j])
            uold[...] = u + alpha * (unew - 2 * u + uold)
            vold[...] = v + alpha * (vnew - 2 * v + vold)
            pold[...] = p + alpha * (pnew - 2 * p + pold)
                    
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
    if L_OUT:
        print("cycle number ", ITMAX)
        print(" diagonal elements of p:\n", pnew.diagonal()[:-1])
        print(" diagonal elements of u:\n", unew.diagonal()[:-1])
        print(" diagonal elements of v:\n", vnew.diagonal()[:-1])
    print("total: ",dt0)
    print("t100: ",dt1)
    print("t200: ",dt2)
    print("t300: ",dt3)


    if VAL:

        u_val_f = 'ref/u.64.64.IT4000.txt'
        v_val_f = 'ref/v.64.64.IT4000.txt'
        p_val_f = 'ref/p.64.64.IT4000.txt'
        uval = np.zeros((M_LEN, N_LEN))
        vval = np.zeros((M_LEN, N_LEN))
        pval = np.zeros((M_LEN, N_LEN))

        uref, vref, pref = read_arrays(u_val_f, v_val_f, p_val_f)
        uval = uref-unew
        vval = vref-vnew
        pval = pref-pnew
        
        uLinfN= np.linalg.norm(uval, np.inf)
        vLinfN= np.linalg.norm(vval, np.inf)
        pLinfN= np.linalg.norm(pval, np.inf)

        

        live_plot_val(uval, vval, pval, "Val")
        print("uLinfN: ", uLinfN)
        print("vLinfN: ", vLinfN)
        print("pLinfN: ", pLinfN)
        print("udiff max: ",uval.max())
        print("vdiff max: ",vval.max())
        print("pdiff max: ",pval.max())

if __name__ == "__main__":
    main()
