using LinearAlgebra

M = 32
N = 32
M_LEN = M + 1
N_LEN = N + 1
ITMAX = 1000
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

global tdt = dt

for ncycle in 1:ITMAX
    @time begin
    global p, u, v, cu, cv, z, h, uold, vold, pold, unew, vnew, pnew, tdt, tdtsdx, tdtsdy, time, tdts8
    # println(size(h[1:end-1, 1:end-1]))
    # println(size(p[1:end-1, 1:end-1]))
    # println(size(u[2:end, 1:end-1]))
    # println(size(u[2:end, 1:end-1]))

    # println(size(u[1:end-1, 1:end-1]))
    # println(size(u[1:end-1, 1:end-1]))
    # println(size(v[1:end-1, 2:end]))
    # println(size(v[1:end-1, 2:end]))
    # println(size(v[1:end-1, 1:end-1]))
    # println(size(v[1,end-1, 1:end-1]))
    #println((fsdx .* (v[2:end, 2:end] .- v[1:end-1, 2:end])))
    #println((fsdy .* (u[2:end, 2:end] .- u[2:end, 1:end-1])))
    #println((p[1:end-1, 1:end-1] .+ p[2:end,1:end-1] .+ p[2:end,2:end] .+ p[1:end-1,2:end]))
    cu[2:end, 1:end-1] = ((p[2:end, 1:end-1] .+ p[1:end-1, 1:end-1]) .* u[2:end,1:end-1]) ./ 2.0
    cv[1:end-1, 2:end] = ((p[1:end-1, 2:end] .+ p[1:end-1, 1:end-1]) .* v[1:end-1,2:end]) ./ 2.0
    z[2:end, 2:end] = (fsdx .* (v[2:end, 2:end] .- v[1:end-1, 2:end]) .- 
                      fsdy .* (u[2:end, 2:end] .- u[2:end, 1:end-1])) ./
                      (p[1:end-1, 1:end-1] .+ p[2:end,1:end-1] .+ 
                      p[2:end,2:end] .+ p[1:end-1,2:end])
    
    h[1:end-1, 1:end-1] = p[1:end-1, 1:end-1] .+ 0.25 .* 
                            (u[2:end, 1:end-1] .* u[2:end, 1:end-1] .+
                            u[1:end-1, 1:end-1] .* u[1:end-1, 1:end-1] .+ 
                            v[1:end-1, 2:end] .* v[1:end-1, 2:end] .+ 
                            v[1:end-1, 1:end-1] .* v[1:end-1, 1:end-1])
    
    end

    cu[1, :] = cu[M, :]
    h[M,:] = h[1, :]
    cv[M, 2:end] = cv[1, 2:end]
    z[1, 2:end] = z[M, 2:end]

    cv[:, 1] = cv[:, N]
    h[:, N] = h[:, 1]
    cu[2:end, N] = cu[2:end, 1]
    z[2:end, 1] = z[2:end, N]

    cu[1, N] = cu[M, 1]
    cv[M, 1] = cv[1, N]
    z[1, 1] = z[M, N]
    h[M, N] = h[1, 1]

    tdts8 = tdt / 8.0
    tdtsdx = tdt / dx
    tdtsdy = tdt / dy

    @time begin
        unew[2:end, 1:end-1] = uold[2:end, 1:end-1] .+ tdts8 .* 
                                (z[2:end, 2:end] .+ z[2:end, 1:end-1]) .*
                                (cv[2:end,2:end] .+ cv[2:end,1:end-1] .+ 
                                cv[1:end-1,2:end] .+ cv[1:end-1,1:end-1]) .-
                                tdtsdx .* (h[2:end, 1:end-1] .- h[1:end-1, 1:end-1])
        vnew[1:end-1, 2:end] = vold[1:end-1, 2:end] .- tdts8 .*
                                (z[2:end, 2:end] .+ z[1:end-1, 2:end]) .*
                                (cu[2:end,2:end] .+ cu[2:end,1:end-1] .+ 
                                cu[1:end-1,2:end] .+ cu[1:end-1,1:end-1]) .-
                                tdtsdy .* (h[1:end-1, 2:end] .- h[1:end-1, 1:end-1])
        pnew[1:end-1, 1:end-1] = pold[1:end-1, 1:end-1] .- tdtsdx .*
                                (cu[2:end, 1:end-1] .- cu[1:end-1, 1:end-1]) .-
                                tdtsdy .* (cv[1:end-1, 2:end] .- cv[1:end-1, 1:end-1])
    end

    unew[1, :] = unew[M, :]
    pnew[M, :] = pnew[1, :]
    vnew[M, 2:end] = vnew[1, 2:end]

    unew[2:end, N] = unew[2:end, 1]
    vnew[:, 1] = vnew[:, N]
    pnew[:, N] = pnew[:, 1]

    unew[1, N] = unew[M, 1]
    vnew[M, 1] = vnew[1, N]
    pnew[M, N] = pnew[1, 1]

    if ncycle>1
        uold[:, :] = u[:, :] .+ alpha .* (unew[:, :] .- 2.0 .* u[:, :] .+ uold[:, :])
        vold[:, :] = v[:, :] .+ alpha .* (vnew[:, :] .- 2.0 .* v[:, :] .+ vold[:, :])
        pold[:, :] = p[:, :] .+ alpha .* (pnew[:, :] .- 2.0 .* p[:, :] .+ pold[:, :])

        u[:, :] = unew[:, :]
        v[:, :] = vnew[:, :]
        p[:, :] = pnew[:, :]
    else
        tdt = tdt + tdt
        uold = copy(u)
        vold = copy(v)
        pold = copy(p)
        u = copy(unew)
        v = copy(vnew)
        p = copy(pnew)
    end

end

if L_OUT
    println(" Cycle number: ", ITMAX)
    println(" u:\n", diag(u)[1:end-1])
    println(" v:\n", diag(v)[1:end-1])
    println(" p:\n", diag(p)[1:end-1])
end
