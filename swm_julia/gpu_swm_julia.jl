using LinearAlgebra
using BenchmarkTools
using Test
using CUDA

M = 64
N = 64
M_LEN = M + 1
N_LEN = N + 1
ITMAX = 4000
dt = 90.
dt = dt
dx = 100000.
dy = 100000.
fsdx = 4. / (dx)
fsdy = 4. / (dy)
a = 1000000.
alpha = 0.001
L_OUT = true

function initialize(u, v, p, M, N, dx, dy, a)
    pi = 4. * atan(1.)
    tpi = 2. * pi
    d_i = tpi / M
    d_j = tpi / N
    el = N * dx
    pcf = (pi * pi * a * a) / (el * el)

    M_LEN = M + 1
    N_LEN = N + 1

    psi = zeros(Float64, M_LEN, N_LEN)
    global p, u, v

    for i in 0:(M_LEN-1)
        for j in 0:(N_LEN-1)
            idx = i * N_LEN + j + 1
            psi[idx] = a * sin((i + 0.5) * d_i) * sin((j + 0.5) * d_j)
            p[idx] = pcf * (cos(2.0 * i * d_i) + cos(2.0 * j * d_j)) + 50000.0
        end
    end

    u[2:end, 1:end-1] = -(psi[2:end, 2:end] .- psi[2:end, 1:end-1]) / dy
    v[1:end-1, 2:end] = (psi[2:end, 2:end] .- psi[1:end-1, 2:end]) / dx

    u[1, :] = u[M, :]
    v[M, 2:end] = v[1, 2:end]
    u[2:end, N] = u[2:end, 1]
    v[:, 1] = v[:, N]

    u[1, N] = u[M, 1]
    v[M, 1] = v[1, N]
    return u, v, p
end

global unew = zeros(Float64, M_LEN, N_LEN)
global vnew = zeros(Float64, M_LEN, N_LEN)
global pnew = zeros(Float64, M_LEN, N_LEN)
global cu = zeros(Float64, M_LEN, N_LEN)
global cv = zeros(Float64, M_LEN, N_LEN)
global z = zeros(Float64, M_LEN, N_LEN)
global h = zeros(Float64, M_LEN, N_LEN)
global u = zeros(Float64, M_LEN, N_LEN)
global v = zeros(Float64, M_LEN, N_LEN)
global p = zeros(Float64, M_LEN, N_LEN)

@time begin
    u, v, p = initialize(u, v, p, M, N, dx, dy, a)
end

global uold = copy(u)
global vold = copy(v)
global pold = copy(p)

if L_OUT
    println(" Number of points in the x direction: ", M)
    println(" Number of points in the y direction: ", N)
    println(" grid spacing in the x direction: ", dx)
    println(" grid spacing in the y direction: ", dy)
    println(" time step: ", dt)
    println(" time filter coefficient: ", alpha)
    println(" Initial p:\n", diag(p)[1:end-1])
    println(" Initial u:\n", diag(u)[1:end-1])
    println(" Initial v:\n", diag(v)[1:end-1])
end

#Allocate and copy to device

global u_d = CUDA.CuArray(u)
global v_d = CUDA.CuArray(v)
global p_d = CUDA.CuArray(p)
global cu_d = CUDA.CuArray(cu)
global cv_d = CUDA.CuArray(cv)
global z_d = CUDA.CuArray(z)
global h_d = CUDA.CuArray(h)
global unew_d = CUDA.CuArray(unew)
global vnew_d = CUDA.CuArray(vnew)
global pnew_d = CUDA.CuArray(pnew)
global uold_d = CUDA.CuArray(uold)
global vold_d = CUDA.CuArray(vold)
global pold_d = CUDA.CuArray(pold)
global fsdx_d = CUDA.CuArray([fsdx])
global fsdy_d = CUDA.CuArray([fsdy])
global alpha_d = CUDA.CuArray([alpha])


CUDA.allowscalar(true)

function calc_cu_cv_z_h!(cu_d, cv_d, z_d, h_d, p_d, u_d, v_d, fsdx_d, fsdy_d)
    CUDA.@sync for i in 1:size(cu_d, 1)-1
        for j in 1:size(cu_d, 2)-1
            cu_d[i+1, j] = ((p_d[i+1, j] + p_d[i, j]) * u_d[i+1, j]) / 2.0
            cv_d[i, j+1] = ((p_d[i, j+1] + p_d[i, j]) * v_d[i, j+1]) / 2.0
            z_d[i+1, j+1] = (fsdx_d[1] * (v_d[i+1, j+1] - v_d[i, j+1]) - 
                            fsdy_d[1] * (u_d[i+1, j+1] - u_d[i+1, j])) /
                            (p_d[i, j] + p_d[i+1, j] + p_d[i+1, j+1] + p_d[i, j+1])
            h_d[i, j] = p_d[i, j] + 0.25 * (u_d[i+1, j] * u_d[i+1, j] +
                            u_d[i, j] * u_d[i, j] + v_d[i, j+1] * v_d[i, j+1] +
                            v_d[i, j] * v_d[i, j])
        end
    end
    return nothing
end

function calc_new_val(unew_s, vnew_d, pnew_d, uold_d, vold_d, pold_d, cu_d, cv_d, z_d, h_d, tdts8_d, tdtsdx_d, tdtsdy_d)
    CUDA.@sync for i in 1:size(unew_d, 1)-1
        for j in 1:size(unew_d, 2)-1
            unew_d[i+1, j] = uold_d[i+1, j] + tdts8_d[1] * 
                            (z_d[i+1, j+1] + z_d[i+1, j]) *
                            (cv_d[i+1, j+1] + cv_d[i+1, j] + 
                            cv_d[i, j+1] + cv_d[i, j]) -
                            tdtsdx_d[1] * (h_d[i+1, j] - h_d[i, j])
            vnew_d[i, j+1] = vold_d[i, j+1] - tdts8_d[1] *
                            (z_d[i+1, j+1] + z_d[i, j+1]) *
                            (cu_d[i+1, j+1] + cu_d[i+1, j] + 
                            cu_d[i, j+1] + cu_d[i, j]) -
                            tdtsdy_d[1] * (h_d[i, j+1] - h_d[i, j])
            pnew_d[i, j] = pold_d[i, j] - tdtsdx_d[1] *
                            (cu_d[i+1, j] - cu_d[i, j]) -
                            tdtsdy_d[1] * (cv_d[i, j+1] - cv_d[i, j])
        end
    end
end

function calc_old_val(uold_d, vold_d, pold_d, unew_d, vnew_d, pnew_d, alpha_d)
    CUDA.@sync for i in 1:size(uold_d, 1)-1
        for j in 1:size(uold_d, 2)-1
            uold_d[i, j] = uold_d[i, j] + alpha_d[1] * (unew_d[i, j] - 2.0 * uold_d[i, j])
            vold_d[i, j] = vold_d[i, j] + alpha_d[1] * (vnew_d[i, j] - 2.0 * vold_d[i, j])
            pold_d[i, j] = pold_d[i, j] + alpha_d[1] * (pnew_d[i, j] - 2.0 * pold_d[i, j])
        end
    end
    return nothing
end

global tdt = dt

for ncycle in 1:ITMAX
    
    global p, u, v, cu, cv, z, h, uold, vold, pold, unew, vnew, pnew, tdt, tdtsdx, tdtsdy, time, tdts8
    global cu_d, p_d, u_d, v_d, cv_d, z_d, h_d, unew_d, vnew_d, pnew_d, uold_d, vold_d, pold_d
    @time begin
        calc_cu_cv_z_h!(cu_d, cv_d, z_d, h_d, p_d, u_d, v_d, fsdx_d, fsdy_d)
    end

    cu_d[1, :] = cu_d[M, :]
    h_d[M,:] = h_d[1, :]
    cv_d[M, 2:end] = cv_d[1, 2:end]
    z_d[1, 2:end] = z_d[M, 2:end]

    cv_d[:, 1] = cv_d[:, N]
    h_d[:, N] = h_d[:, 1]
    cu_d[2:end, N] = cu_d[2:end, 1]
    z_d[2:end, 1] = z_d[2:end, N]

    cu_d[1, N] = cu_d[M, 1]
    cv_d[M, 1] = cv_d[1, N]
    z_d[1, 1] = z_d[M, N]
    h_d[M, N] = h_d[1, 1]

    tdts8 = tdt / 8.0
    tdtsdx = tdt / dx
    tdtsdy = tdt / dy
    
    tdts8_d = CUDA.CuArray([tdts8])
    tdtsdx_d = CUDA.CuArray([tdtsdx])
    tdtsdy_d = CUDA.CuArray([tdtsdy])

    @time begin
        calc_new_val(unew_d, vnew_d, pnew_d, uold_d, vold_d, pold_d, cu_d, cv_d, z_d, h_d, tdts8_d, tdtsdx_d, tdtsdy_d)
    end
    
    unew_d[1, :] = unew_d[M, :]
    pnew_d[M, :] = pnew_d[1, :]
    vnew_d[M, 2:end] = vnew_d[1, 2:end]

    unew_d[2:end, N] = unew_d[2:end, 1]
    vnew_d[:, 1] = vnew_d[:, N]
    pnew_d[:, N] = pnew_d[:, 1]

    unew_d[1, N] = unew_d[M, 1]
    vnew_d[M, 1] = vnew_d[1, N]
    pnew_d[M, N] = pnew_d[1, 1]

    if ncycle>1
        @time begin
            calc_old_val(uold_d, vold_d, pold_d, unew_d, vnew_d, pnew_d, alpha_d)
        end
        u_d = unew_d
        v_d = vnew_d
        p_d = pnew_d
    else
        tdt = tdt + tdt
        uold_d = copy(u_d)
        vold_d = copy(v_d)
        pold_d = copy(p)
        u_d = copy(unew_d)
        v_d = copy(vnew_d)
        p_d = copy(pnew_d)
    end

end

if L_OUT
    u = Array(u_d)
    v = Array(v_d)
    p = Array(p_d)
    println(" Cycle number: ", ITMAX)
    println(" u:\n", diag(u)[1:end-1])
    println(" v:\n", diag(v)[1:end-1])
    println(" p:\n", diag(p)[1:end-1])
end
