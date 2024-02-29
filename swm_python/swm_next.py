import numpy as np
import argparse
import matplotlib.pyplot as plt
import numpy as np
import gt4py.next as gtx
import gt4py.cartesian.gtscript as gtscript
from time import perf_counter
from gt4py.next import Field
import time as t
#import cupy
#import gt4py

# Initialize model parameters
M = 64 # args.M
N = 64 # args.N
M_LEN = M + 1
N_LEN = N + 1
L_OUT = True # args.L_OUT
VIS = False
VIS_DT=10
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
pcf = (pi * pi * a * a) / (el * el)
SIZE = M_LEN * N_LEN
dt0=0.
dt1=0.
dt2=0.
dt3=0.


# Model Variables
u = np.zeros((M_LEN, N_LEN, 1))
v = np.zeros((M_LEN, N_LEN, 1))
p = np.zeros((M_LEN, N_LEN, 1))
unew = np.zeros((M_LEN, N_LEN, 1))
vnew = np.zeros((M_LEN, N_LEN, 1))
pnew = np.zeros((M_LEN, N_LEN, 1))
uold = np.zeros((M_LEN, N_LEN, 1))
vold = np.zeros((M_LEN, N_LEN, 1))
pold = np.zeros((M_LEN, N_LEN, 1))
uvis = np.zeros((M_LEN, N_LEN, 1))
vvis = np.zeros((M_LEN, N_LEN, 1))
pvis = np.zeros((M_LEN, N_LEN, 1))
cu = np.zeros((M_LEN, N_LEN, 1))
cv = np.zeros((M_LEN, N_LEN, 1))
z = np.zeros((M_LEN, N_LEN, 1))
h = np.zeros((M_LEN, N_LEN, 1))
psi = np.zeros((M_LEN, N_LEN, 1))

from IPython.display import clear_output
from matplotlib import pyplot as plt
    

def live_plot3(fu, fv, fp, title=''):
    clear_output(wait=True)
    fig, (ax1, ax2, ax3) = plt.subplots(figsize=(13, 3), ncols=3)

    pos1 = ax1.imshow(fp, cmap='Blues', vmin=49999, vmax=50001,interpolation='none')
    ax1.set_title('p')
    pos2 = ax2.imshow(fu, cmap='Reds', vmin=-1, vmax=1,interpolation='none')
    ax2.set_title('u')
    pos3 = ax3.imshow(fv, cmap='Greens',vmin=-1, vmax=1,interpolation='none')
    ax3.set_title('v')

    fig.suptitle(title)
    #plt.xlabel('x')
    #plt.ylabel('y')
    plt.show()

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

# Initial values of the stream function and p
psi[:,:,0] = a * np.sin((np.arange(0,M_LEN)[:,np.newaxis]+0.5) * d_i) * np.sin((np.arange(0,N_LEN)+ .5) * d_j)
p[:,:,0]   = pcf * (np.cos(2. * np.arange(0,M_LEN)[:, np.newaxis] * d_i) + np.cos(2. * np.arange(0,N_LEN) * d_j)) + 50000.
            
# Calculate initial u and v
u[1:,:-1,0] = -(psi[1:,1:,0] - psi[1:,:-1,0]) / dy
v[:-1,1:,0] =  (psi[1:,1:,0] - psi[:-1,1:,0]) / dx

            

if VIS==True:
    live_plot3(u,v,p, "init")
    print(p.max())
    print(p.min())
    print(u.max())
    print(u.min())
    print(v.max())
    print(v.min())

# Periodic Boundary conditions

u[0, :,0] = u[M, :,0]
v[M, 1:,0] = v[0, 1:,0]
u[1:, N,0] = u[1:, 0,0]
v[:, 0,0] = v[:, N,0]

u[0, N,0] = u[M, 0,0]
v[M, 0,0] = v[0, N,0]


if VIS==True:
    live_plot3(u,v,p, "Periodic Bounday Conditions")
    
# Save initial conditions
uold = np.copy(u[...])
vold = np.copy(v[...])
pold = np.copy(p[...])

# Print initial conditions
if L_OUT:
    print(" Number of points in the x direction: ", M)
    print(" Number of points in the y direction: ", N)
    print(" grid spacing in the x direction: ", dx)
    print(" grid spacing in the y direction: ", dy)
    print(" time step: ", dt)
    print(" time filter coefficient: ", alpha)
    print(" Initial p:\n", p[:,:,0].diagonal()[:-1])
    print(" Initial u:\n", u[:,:,0].diagonal()[:-1])
    print(" Initial v:\n", v[:,:,0].diagonal()[:-1])

nx = M
ny = N
nz = 1
dtype = np.float64
#gt4py_type = "cartesian"
gt4py_type = "next"
allocator = gtx.itir_python

I = gtx.Dimension("I")
J = gtx.Dimension("J")
K = gtx.Dimension("K", kind = gtx.DimensionKind.VERTICAL)

domain = gtx.domain({I:nx+1, J:ny+1, K:nz})

h_gt = gtx.as_field(domain,h,allocator=allocator)
p_gt = gtx.as_field(domain,p,allocator=allocator)
u_gt = gtx.as_field(domain,u,allocator=allocator)
v_gt = gtx.as_field(domain,v,allocator=allocator)
z_gt = gtx.as_field(domain,z,allocator=allocator)
cu_gt = gtx.as_field(domain,cu,allocator=allocator)
cv_gt = gtx.as_field(domain,cv,allocator=allocator)
pnew_gt = gtx.as_field(domain,pnew,allocator=allocator)
unew_gt = gtx.as_field(domain,unew,allocator=allocator)
vnew_gt = gtx.as_field(domain,vnew,allocator=allocator)
uold_gt = gtx.as_field(domain,uold,allocator=allocator)
vold_gt = gtx.as_field(domain,vold,allocator=allocator)
pold_gt = gtx.as_field(domain,pold,allocator=allocator)

cartesian_backend = "numpy"
next_backend = gtx.itir_python

if gt4py_type == "cartesian":
    # i --> 0,M
    #j --> 0,N
    # at nx+1 its the boundary region, mask one region and 
    @gtscript.stencil(backend=cartesian_backend)
    def calc_h(
        p: gtscript.Field[dtype],
        u: gtscript.Field[dtype],
        v: gtscript.Field[dtype],
        h: gtscript.Field[dtype]
    ):
        with computation(PARALLEL), interval(...):
            h = p + 0.25 * u[1,0,0] * u[1,0,0] + u * u + v[0,1,0] * v[0,1,0] + v * v
    
    #nx = M
    #ny = N
    #nz = 1
    # i --> 1,M+1  (1,1,0) (nx,ny,nz)
    #j --> 1,N+1
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
            cu = .5 * (p + p) * u
            
    @gtscript.stencil(backend=cartesian_backend)
    def calc_cv(
        v: gtscript.Field[dtype],
        p: gtscript.Field[dtype],
        cv: gtscript.Field[dtype]
    ):
        with computation(PARALLEL), interval(...):
            cv = .5 * (p + p) * v
    
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
    
    t0_start = perf_counter()
    time = 0.0
    # Main time loop
    for ncycle in range(ITMAX):
        if((ncycle%100==0) & (VIS==False)):
            print(f"cycle number{ncycle} and gt4py type {gt4py_type}")
        t1_start = perf_counter()
        # Calculate cu, cv, z, and h
        calc_h(p=p_gt, u=u_gt, v=v_gt, h=h_gt, origin=(0,0,0), domain=(nx,ny,nz)) 
        h = h_gt.asnumpy()

        calc_z(fsdx=fsdx, fsdy=fsdy, u=u_gt, v=v_gt, p=p_gt, z=z_gt, origin=(1,1,0), domain=(nx,ny,nz)) # domain(nx+1,ny+1,nz) gives error why?
        z = z_gt.asnumpy()

        calc_cu(u=u_gt, p=p_gt, cu=cu_gt, origin=(1,0,0), domain=(nx,ny,nz)) # (nx,ny+1,nz)-->works domain(nx+1,ny+1,nz) gives error why? try removing ny+1
        cu = cu_gt.asnumpy()

        calc_cv(v=v_gt, p=p_gt, cv=cv_gt, origin=(0,1,0), domain=(nx,ny,nz)) #(nx+1,ny,nz)--> works domain(nx+1,ny+1,nz) gives error why?
        cv = cv_gt.asnumpy()
        t1_stop = perf_counter()
        dt1 = dt1 + (t1_stop - t1_start)
    
        # # Periodic Boundary conditions
        #try region
        cu[0, :,0] = cu[M, :,0]
        h[M, :,0] = h[0, :,0]
        cv[M, 1:,0] = cv[0, 1:,0]
        z[0, 1:,0] = z[M, 1:,0]
        
        cv[:, 0,0] = cv[:, N,0]
        h[:, N,0] = h[:, 0,0]
        cu[1:, N,0] = cu[1:, 0,0]
        z[1:, N,0] = z[1:, 0,0]
            
        cu[0, N,0] = cu[M, 0,0]
        cv[M, 0,0] = cv[0, N,0]
        z[0, 0,0] = z[M, N,0]
        h[M, N,0] = h[0, 0,0]
            
        # Calclulate new values of u,v, and p
        tdts8 = tdt / 8.
        tdtsdx = tdt / dx
        tdtsdy = tdt / dy
        #print(tdts8, tdtsdx, tdtsdy)
    
        t2_start = perf_counter()
        calc_unew(tdts8=tdts8, tdtsdx=tdtsdx, uold=uold_gt, cu=cu_gt, cv=cv_gt, z=z_gt, h=h_gt, unew=unew_gt, origin=(1,0,0), domain=(nx,ny,nz))
        unew = unew_gt.asnumpy()
        
        calc_vnew(tdts8=tdts8, tdtsdy=tdtsdy, vold=vold_gt, cu=cu_gt, cv=cv_gt, z=z_gt, h=h_gt, vnew=vnew_gt, origin=(0,1,0), domain=(nx,ny,nz))
        vnew = vnew_gt.asnumpy()
        
        calc_pnew(tdtsdx=tdtsdx, tdtsdy=tdtsdy, pold=pold_gt, cu=cu_gt, cv=cv_gt, pnew=pnew_gt, origin=(0,0,0), domain=(nx,ny,nz))
        pnew = pnew_gt.asnumpy()
        t2_stop = perf_counter()
        dt2 = dt2 + (t2_stop - t2_start)
            
        # for i in range(M):
        #     for j in range(N):
        #         unew[i+1,j,0] = uold[i+1,j,0] + tdts8 * (z[i+1,j+1,0] + z[i+1,j,0]) * (cv[i+1,j+1,0] + cv[i+1,j,0] + cv[i,j+1,0] + cv[i,j,0]) - tdtsdx * (h[i+1,j,0] - h[i,j,0])
        #         vnew[i,j+1,0] = vold[i,j+1,0] - tdts8 * (z[i+1,j+1,0] + z[i,j+1,0]) * (cu[i+1,j+1,0] + cu[i+1,j,0] + cu[i,j+1,0] + cu[i,j,0]) - tdtsdy * (h[i,j+1,0] - h[i,j,0])
        #         pnew[i,j,0] = pold[i,j,0] - tdtsdx * (cu[i+1,j,0] - cu[i,j,0]) - tdtsdy * (cv[i,j+1,0] - cv[i,j,0])
                           
        # Periodic Boundary conditions
        unew[0, :,0] = unew[M, :,0]
        pnew[M, :,0] = pnew[0, :,0]
        vnew[M, 1:,0] = vnew[0, 1:,0]
        unew[1:, N,0] = unew[1:, 0,0]
        vnew[:, 0,0] = vnew[:, N,0]
        pnew[:, N,0] = pnew[:, 0,0]
        
        unew[0, N,0] = unew[M, 0,0]
        vnew[M, 0,0] = vnew[0, N,0]
        pnew[M, N,0] = pnew[0, 0,0]
        
        time = time + dt
    
        if(ncycle > 0):
            t3_start = perf_counter()
            calc_pold(p=p_gt, alpha=alpha, pnew=pnew_gt, pold=pold_gt, origin=(0,0,0), domain=(nx,ny,nz))
            pold = pold_gt.asnumpy()
            calc_uold(u=u_gt, alpha=alpha, unew=unew_gt, uold=uold_gt, origin=(0,0,0), domain=(nx,ny,nz))
            uold = uold_gt.asnumpy()
            calc_vold(v=v_gt, alpha=alpha, vnew=vnew_gt, vold=vold_gt, origin=(0,0,0), domain=(nx,ny,nz))
            vold = vold_gt.asnumpy()
            
            copy_var(unew_gt, u_gt, origin=(0,0,0), domain=(nx,ny,nz))
            copy_var(vnew_gt, v_gt, origin=(0,0,0), domain=(nx,ny,nz))
            copy_var(pnew_gt, p_gt, origin=(0,0,0), domain=(nx,ny,nz))

            t3_stop = perf_counter()
            dt3 = dt3 + (t3_stop - t3_start)
    
        else:
            tdt = tdt+tdt
    
            uold = np.copy(u[...])
            vold = np.copy(v[...])
            pold = np.copy(p[...])
            u = np.copy(unew[...])
            v = np.copy(vnew[...])
            p = np.copy(pnew[...])
    
        if((VIS == True) & (ncycle%VIS_DT==0)):
            live_plot3(u, v, p, "ncycle: " + str(ncycle))
            
    t0_stop = perf_counter()
    dt0 = dt0 + (t0_stop - t0_start)
    # Print initial conditions
    if L_OUT:
           print("cycle number ", ITMAX)
           print(" diagonal elements of p:\n", pnew[:,:,0].diagonal()[:-1])
           print(" diagonal elements of u:\n", unew[:,:,0].diagonal()[:-1])
           print(" diagonal elements of v:\n", vnew[:,:,0].diagonal()[:-1])
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
        uval = uref-unew[:,:,0]
        vval = vref-vnew[:,:,0]
        pval = pref-pnew[:,:,0]
        
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


# gt4py NEXT part!!!!!!!!!!!!!!!!!!!!


if gt4py_type == "next":
    Ioff = gtx.FieldOffset("I", source=I, target=(I,))
    Joff = gtx.FieldOffset("J", source=J, target=(J,))

    @gtx.field_operator
    def calc_h(p: Field[[I,J,K],dtype], u: Field[[I,J,K],dtype], v: Field[[I,J,K],dtype]) -> Field[[I,J,K],dtype]:
        return (p + 0.25 * u(Ioff[1]) * u(Ioff[1]) + u * u + v(Joff[1]) * v(Joff[1]) + v * v)
    #u[1,0,0] --> u(Ioff[1])
    #u[1,1,0] --> u(Ioff[1])(Joff[-1])
    
    @gtx.program(backend=next_backend)
    def calc_h_program(
        p: Field[[I,J,K],dtype],
        u: Field[[I,J,K],dtype],
        v: Field[[I,J,K],dtype],
        h: Field[[I,J,K],dtype]
    ):
        calc_h(p,u,v, out=h[:-1,:-1,:])
      
    @gtx.field_operator
    def calc_z(p: Field[[I,J,K],dtype], u: Field[[I,J,K],dtype], v: Field[[I,J,K],dtype],fsdx: float,fsdy: float) -> Field[[I,J,K],dtype]:
        return ((fsdx * (v - v(Ioff[-1])) - fsdy * (u - u(Joff[-1]))) / (p(Ioff[-1])(Joff[-1]) + p(Joff[-1]) + p + p(Ioff[-1])))
    
    @gtx.program(backend=next_backend)
    def calc_z_program(
        p: Field[[I,J,K],dtype],
        u: Field[[I,J,K],dtype],
        v: Field[[I,J,K],dtype],
        z: Field[[I,J,K],dtype],
        fsdx: float,
        fsdy: float
    ):
        calc_z(p,u,v,fsdx,fsdy, out=z[1:-1,1:-1,:])
        
    @gtx.field_operator
    def calc_cu(p: Field[[I,J,K],dtype], u: Field[[I,J,K],dtype]) -> Field[[I,J,K],dtype]:
        return (.5 * (p + p) * u)
    @gtx.program(backend=next_backend)
    def calc_cu_program(
        p: Field[[I,J,K],dtype],
        u: Field[[I,J,K],dtype],
        cu: Field[[I,J,K],dtype],
    ):
        calc_cu(p,u, out=cu[1:-1,:-1,:])
    
    @gtx.field_operator
    def calc_cv(p: Field[[I,J,K],dtype], v: Field[[I,J,K],dtype]) -> Field[[I,J,K],dtype]:
        return (.5 * (p + p) * v)
    @gtx.program(backend=next_backend)
    def calc_cv_program(
        p: Field[[I,J,K],dtype],
        v: Field[[I,J,K],dtype],
        cv: Field[[I,J,K],dtype],
    ):
        calc_cv(p,v, out=cv[:-1,1:-1,:])
            
    @gtx.field_operator
    def calc_pnew(
        pold: Field[[I,J,K],dtype], 
        cu: Field[[I,J,K],dtype], 
        cv: Field[[I,J,K],dtype], 
        tdtsdx: float, 
        tdtsdy: float) -> Field[[I,J,K],dtype]:
        return (pold - tdtsdx * (cu(Ioff[1]) - cu) - tdtsdy * (cv(Joff[1]) - cv))
    @gtx.program(backend=next_backend)
    def calc_pnew_program(
        pold: Field[[I,J,K],dtype],
        cu: Field[[I,J,K],dtype],
        cv: Field[[I,J,K],dtype],
        pnew: Field[[I,J,K],dtype],
        tdtsdx: float,
        tdtsdy: float
    ):
        calc_pnew(pold, cu, cv, tdtsdx, tdtsdy, out=pnew[:-1,:-1,:])
        
    @gtx.field_operator
    def calc_unew(
        uold: Field[[I,J,K],dtype], 
        cu: Field[[I,J,K],dtype], 
        cv: Field[[I,J,K],dtype], 
        z: Field[[I,J,K],dtype], 
        h: Field[[I,J,K],dtype], 
        tdts8: float, 
        tdtsdx: float) -> Field[[I,J,K],dtype]:
        return (uold + tdts8 * (z + z(Joff[1])) * (cv(Joff[1]) + cv + cv(Ioff[-1])(Joff[1]) + cv(Ioff[-1])) - tdtsdx * (h - h(Ioff[-1])))
    
    @gtx.program(backend=next_backend)
    def calc_unew_program(
        uold: Field[[I,J,K],dtype],
        cu: Field[[I,J,K],dtype],
        cv: Field[[I,J,K],dtype],
        z: Field[[I,J,K],dtype],
        h: Field[[I,J,K],dtype],
        unew: Field[[I,J,K],dtype],
        tdts8: float,
        tdtsdx: float
    ):
        calc_unew(uold, cu, cv, z, h, tdts8, tdtsdx, out=unew[1:-1,:-1,:])
        
    @gtx.field_operator
    def calc_vnew(
        vold: Field[[I,J,K],dtype], 
        cu: Field[[I,J,K],dtype], 
        cv: Field[[I,J,K],dtype], 
        z: Field[[I,J,K],dtype], 
        h: Field[[I,J,K],dtype], 
        tdts8: float, 
        tdtsdy: float) -> Field[[I,J,K],dtype]:
        return (vold - tdts8 * (z(Ioff[1]) + z) * (cu(Ioff[1]) + cu(Ioff[1])(Joff[-1]) + cu + cu(Joff[-1])) - tdtsdy * (h - h(Joff[-1])))
    @gtx.program(backend=next_backend)
    def calc_vnew_program(
        vold: Field[[I,J,K],dtype],
        cu: Field[[I,J,K],dtype],
        cv: Field[[I,J,K],dtype],
        z: Field[[I,J,K],dtype],
        h: Field[[I,J,K],dtype],
        vnew: Field[[I,J,K],dtype],
        tdts8: float,
        tdtsdy: float
    ):
        calc_vnew(vold, cu, cv, z, h, tdts8, tdtsdy, out=vnew[:-1,1:-1,:])
        
    @gtx.field_operator
    def calc_pold(
        p: Field[[I,J,K],dtype], 
        alpha: float, 
        pnew: Field[[I,J,K],dtype],
        pold: Field[[I,J,K],dtype]) -> Field[[I,J,K],dtype]:
        return (p + alpha * (pnew - 2. * p + pold))
    @gtx.program(backend=next_backend)
    def calc_pold_program(
        p: Field[[I,J,K],dtype],
        alpha: float,
        pnew: Field[[I,J,K],dtype],
        pold: Field[[I,J,K],dtype]
    ):
        calc_pold(p, alpha, pnew, pold, out=pold[:-1,:-1,:])
        
    @gtx.field_operator
    def calc_uold(
        u: Field[[I,J,K],dtype], 
        alpha: float, 
        unew: Field[[I,J,K],dtype],
        uold: Field[[I,J,K],dtype]) -> Field[[I,J,K],dtype]:
        return (u + alpha * (unew - 2. * u + uold))
    @gtx.program(backend=next_backend)
    def calc_uold_program(
        u: Field[[I,J,K],dtype],
        alpha: float,
        unew: Field[[I,J,K],dtype],
        uold: Field[[I,J,K],dtype]
    ):
        calc_uold(u, alpha, unew, uold, out=uold[:-1,:-1,:])
        
    @gtx.field_operator
    def calc_vold(
        v: Field[[I,J,K],dtype], 
        alpha: float, 
        vnew: Field[[I,J,K],dtype],
        vold: Field[[I,J,K],dtype]) -> Field[[I,J,K],dtype]:
        return (v + alpha * (vnew - 2. * v + vold))
    @gtx.program(backend=next_backend)
    def calc_vold_program(
        v: Field[[I,J,K],dtype],
        alpha: float,
        vnew: Field[[I,J,K],dtype],
        vold: Field[[I,J,K],dtype]
    ):
        calc_vold(v, alpha, vnew, vold, out=vold[:-1,:-1,:])
    
    time = 0.0 
    # Main time loop
    for ncycle in range(ITMAX):
        if((ncycle%100==0) & (VIS==False)):
            print(f"cycle number {ncycle} and gt4py type {gt4py_type}")
        # Calculate cu, cv, z, and h
        #for i in range(M):
        #    for j in range(N):
        #        h[i, j,0] = p[i, j,0] + 0.25 * (u[i + 1, j,0] * u[i + 1, j,0] + u[i, j,0] * u[i, j,0] +
        #                                v[i, j + 1,0] * v[i, j + 1,0] + v[i, j,0] * v[i, j,0])

        start_time = t.time()
        calc_h_program(p=p_gt, u=u_gt, v=v_gt, h=h_gt, offset_provider={"Ioff":I, "Joff":J})
        h = h_gt.asnumpy()
        calc_z_program(p=p_gt, u=u_gt, v=v_gt, z=z_gt, fsdx=fsdx, fsdy=fsdy, offset_provider={"Ioff":I, "Joff":J})
        z = z_gt.asnumpy()
        calc_cu_program(p=p_gt, u=u_gt, cu=cu_gt, offset_provider={"Ioff":I, "Joff":J})
        cu = cu_gt.asnumpy()
        calc_cv_program(p=p_gt, v=v_gt, cv=cv_gt, offset_provider={"Ioff":I, "Joff":J})
        cv = cv_gt.asnumpy()
    
        #cu[1:,:-1] = .5 * (p[1:,:-1] + p[:-1,:-1]) * u[1:,:-1]
        #cv[:-1,1:] = .5 * (p[:-1,1:] + p[:-1,:-1]) * v[:-1,1:]
        #z[1:,1:] = (fsdx * (v[1:,1:] - v[:-1,1:]) -
        #            fsdy * (u[1:,1:] - u[1:,:-1] )
        #            ) / (p[:-1,:-1] + p[1:,:-1] + p[1:,1:] + p[:-1,1:])
        # for i in range(M):
        #     for j in range(N):
        #         cu[i + 1, j,0] = .5 * (p[i + 1, j,0] + p[i, j,0]) * u[i + 1, j,0]
        #         cv[i, j + 1,0] = .5 * (p[i, j + 1,0] + p[i, j,0]) * v[i, j + 1,0]
        #         z[i + 1, j + 1,0] = (fsdx * (v[i + 1, j + 1,0] - v[i, j + 1,0]) -
        #                         fsdy * (u[i + 1, j + 1,0] - u[i+1, j,0] )
        #                         ) / (p[i, j,0] + p[i + 1, j,0] + p[i + 1, j + 1,0] + p[i, j + 1,0])
    
            # # Periodic Boundary conditions
        cu[0, :,0] = cu[M, :,0]
        h[M, :,0] = h[0, :,0]
        cv[M, 1:,0] = cv[0, 1:,0]
        z[0, 1:,0] = z[M, 1:,0]
        
        cv[:, 0,0] = cv[:, N,0]
        h[:, N,0] = h[:, 0,0]
        cu[1:, N,0] = cu[1:, 0,0]
        z[1:, N,0] = z[1:, 0,0]
            
        cu[0, N,0] = cu[M, 0,0]
        cv[M, 0,0] = cv[0, N,0]
        z[0, 0,0] = z[M, N,0]
        h[M, N,0] = h[0, 0,0]
            
        # Calclulate new values of u,v, and p
        tdts8 = tdt / 8.
        tdtsdx = tdt / dx
        tdtsdy = tdt / dy
        #print(tdts8, tdtsdx, tdtsdy)
    
        #unew[1:,:-1] = uold[1:,:-1] + tdts8 * (z[1:,1:] + z[1:,:-1]) * (cv[1:,1:] + cv[1:,:-1] + cv[:-1,1:] + cv[:-1,:-1]) - tdtsdx * (h[1:,:-1] - h[:-1,:-1])
        #vnew[:-1,1:] = vold[:-1,1:] - tdts8 * (z[1:,1:] + z[:-1,1:]) * (cu[1:,1:] + cu[1:,:-1] + cu[:-1,1:] + cu[:-1,:-1]) - tdtsdy * (h[:-1,1:] - h[:-1,:-1])
        #pnew[:-1,:-1] = pold[:-1,:-1] - tdtsdx * (cu[1:,:-1] - cu[:-1,:-1]) - tdtsdy * (cv[:-1,1:] - cv[:-1,:-1])
        calc_unew_program(uold=uold_gt, cu=cu_gt, cv=cv_gt, z=z_gt, h=h_gt, tdts8=tdts8, tdtsdx=tdtsdx, unew=unew_gt, offset_provider={"Ioff":I, "Joff":J})
        unew = unew_gt.asnumpy()
        calc_vnew_program(vold=vold_gt, cu=cu_gt, cv=cv_gt, z=z_gt, h=h_gt, tdts8=tdts8, tdtsdy=tdtsdy, vnew=vnew_gt, offset_provider={"Ioff":I, "Joff":J})
        vnew = vnew_gt.asnumpy()
        calc_pnew_program(pold=p_gt, cu=cu_gt, cv=cv_gt, pnew=pnew_gt, tdtsdx=tdtsdx, tdtsdy=tdtsdy, offset_provider={"Ioff":I, "Joff":J})
        pnew = pnew_gt.asnumpy()

        # for i in range(M):
        #     for j in range(N):
        #         unew[i+1,j,0] = uold[i+1,j,0] + tdts8 * (z[i+1,j+1,0] + z[i+1,j,0]) * (cv[i+1,j+1,0] + cv[i+1,j,0] + cv[i,j+1,0] + cv[i,j,0]) - tdtsdx * (h[i+1,j,0] - h[i,j,0])
        #         vnew[i,j+1,0] = vold[i,j+1,0] - tdts8 * (z[i+1,j+1,0] + z[i,j+1,0]) * (cu[i+1,j+1,0] + cu[i+1,j,0] + cu[i,j+1,0] + cu[i,j,0]) - tdtsdy * (h[i,j+1,0] - h[i,j,0])
        #         pnew[i,j,0] = pold[i,j,0] - tdtsdx * (cu[i+1,j,0] - cu[i,j,0]) - tdtsdy * (cv[i,j+1,0] - cv[i,j,0])
                    
        
        # Periodic Boundary conditions
        unew[0, :,0] = unew[M, :,0]
        pnew[M, :,0] = pnew[0, :,0]
        vnew[M, 1:,0] = vnew[0, 1:,0]
        unew[1:, N,0] = unew[1:, 0,0]
        vnew[:, 0,0] = vnew[:, N,0]
        pnew[:, N,0] = pnew[:, 0,0]
        
        unew[0, N,0] = unew[M, 0,0]
        vnew[M, 0,0] = vnew[0, N,0]
        pnew[M, N,0] = pnew[0, 0,0]
        
        time = time + dt
    
        if(ncycle > 0):
            #uold[...] = u + alpha * (unew - 2 * u + uold)
            #vold[...] = v + alpha * (vnew - 2 * v + vold)
            #pold[...] = p + alpha * (pnew - 2 * p + pold)
            calc_uold_program(u=u_gt, alpha=alpha, unew=unew_gt, uold=uold_gt, offset_provider={"Ioff":I, "Joff":J})
            uold = uold_gt.asnumpy()
            calc_vold_program(v=v_gt, alpha=alpha, vnew=vnew_gt, vold=vold_gt, offset_provider={"Ioff":I, "Joff":J})
            vold = vold_gt.asnumpy()
            calc_pold_program(p=p_gt, alpha=alpha, pnew=pnew_gt, pold=pold_gt, offset_provider={"Ioff":I, "Joff":J})
            pold = pold_gt.asnumpy()

            u[...] = unew
            v[...] = vnew
            p[...] = pnew
    
        else:
            tdt = tdt+tdt
    
            uold = np.copy(u[...])
            vold = np.copy(v[...])
            pold = np.copy(p[...])
            u = np.copy(unew[...])
            v = np.copy(vnew[...])
            p = np.copy(pnew[...])

        end_time = t.time()
        elapsed_time = end_time - start_time
        print(f"Time for calc total: {elapsed_time}")
        if((VIS == True) & (ncycle%VIS_DT==0)):
            live_plot3(u, v, p, "ncycle: " + str(ncycle))
    # Print initial conditions
    if L_OUT:
           print("cycle number ", ITMAX)
           print(" diagonal elements of p:\n", pnew[:,:,0].diagonal()[:-1])
           print(" diagonal elements of u:\n", unew[:,:,0].diagonal()[:-1])
           print(" diagonal elements of v:\n", vnew[:,:,0].diagonal()[:-1])
           
    if VAL:
        print("val")
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
        
        uLinfN= np.linalg.norm(uval.flatten(), np.inf)
        vLinfN= np.linalg.norm(vval.flatten(), np.inf)
        pLinfN= np.linalg.norm(pval.flatten(), np.inf)

        

        #live_plot_val(uval, vval, pval, "Val")
        print("uLinfN: ", uLinfN)
        print("vLinfN: ", vLinfN)
        print("pLinfN: ", pLinfN)
        print("udiff max: ",uval.max())
        print("vdiff max: ",vval.max())
        print("pdiff max: ",pval.max())
