using Profile
using LinearAlgebra
using LoopVectorization

const M::Int = 256
const N::Int = 256
const M_LEN::Int = M + 1
const N_LEN::Int = N + 1
const ITMAX::Int = 40

const dt::Float64 = 90.::Float64
const dx::Float64 = 100000.::Float64
const dy::Float64 = 100000.::Float64
const fsdx::Float64 = 4.::Float64 / (dx)
const fsdy::Float64 = 4.::Float64 / (dy)
const a::Float64 = 1000000.::Float64
const alpha::Float64 = 0.001::Float64
L_OUT = false


function initialize(u, v, p, M, N, dx, dy, a)
    pi = 4.::Float64 * atan(1.)
    tpi = 2. * pi
    d_i = tpi / M
    d_j = tpi / N
    el = N * dx
    pcf = (pi * pi * a * a) / (el * el)

    M_LEN = M + 1
    N_LEN = N + 1

    psi = zeros(Float64, M_LEN, N_LEN)

    for i in 0:(M_LEN-1)
        for j in 0:(N_LEN-1)
            idx = i * N_LEN + j + 1
            psi[idx] = a * sin((i + 0.5) * d_i) * sin((j + 0.5) * d_j)
            p[idx] = pcf * (cos(2.0 * i * d_i) + cos(2.0 * j * d_j)) + 50000.0::Float64
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

uold = u
vold = v
pold = p

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

function turbo_dot(x::AbstractVector{Float64}, A::AbstractMatrix{Float64}, y::AbstractVector{Float64})
    s = 0.0;
    @turbo for j in eachindex(y), i in eachindex(x)
         s += x[i]*A[i,j]*y[j];
    end
    return s;
end

function calc_cu!(p::Array{Float64},u::Array{Float64},cu::Array{Float64})
   
   # cu[2:end, 1:end-1] = ((p[2:end, 1:end-1] .+ p[1:end-1,1:end-1]) .* u[2:end,1:end-1]) ./ 2.0::Float64
   @turbo for j = 1:N,i = 2:M+1
      @inbounds cu[i,j] = ((p[i,j] + p[i-1,j]) * u[i,j]) / 2.0
   end
   return
end
function calc_cv!(p::Array{Float64},v::Array{Float64},cv::Array{Float64})
   
   # cv[1:end-1, 2:end] = ((p[1:end-1, 2:end] .+ p[1:end-1,1:end-1]) .* v[1:end-1,2:end]) ./ 2.0::Float64
   @turbo for j = 2:N+1,i = 1:M
      @inbounds cv[i,j] = ((p[i,j] + p[i,j-1]) * v[i,j]) / 2.0
   end
   return
end
function calc_z!(v::Array{Float64},u::Array{Float64},p::Array{Float64},z::Array{Float64},fsdx::Float64,fsdy::Float64)
     # z[2:end, 2:end] = (fsdx .* (v[2:end, 2:end] .- v[1:end-1, 2:end]) .- 
     #                  fsdy .* (u[2:end, 2:end] .- u[2:end, 1:end-1])) ./
     #                  (p[1:end-1,1:end-1] .+ p[2:end,1:end-1] .+ 
     #                 p[2:end,2:end] .+ p[1:end-1,2:end])
     @turbo for j=2:N,i=2:M
       @inbounds z[i,j] = (fsdx * (v[i,j] - v[i-1,j]) - fsdy * (u[i,j]-u[i,j-1]))/
			(p[i-1,j-1] + p[i,j-1] + p[i,j] + p[i-1,j])
     end
     return
end
function calc_h!(p::Array{Float64},u::Array{Float64},v::Array{Float64},h::Array{Float64})
     # h[1:end-1,1:end-1] = p[1:end-1,1:end-1] .+ 0.25::Float64 .* 
     #                        (u[2:end, 1:end-1] .* u[2:end, 1:end-1] .+
     #                        u[1:end-1,1:end-1] .* u[1:end-1,1:end-1] .+ 
     #                        v[1:end-1, 2:end] .* v[1:end-1, 2:end] .+ 
     #                       v[1:end-1,1:end-1] .* v[1:end-1,1:end-1])
     @turbo for j=1:N,i=1:M
       @inbounds h[i,j] = p[i,j] + 0.25* (u[i+1,j]^2 + u[i,j]^2 + v[i,j+2]^2 + v[i,j]^2)
     end
end
function calc_unew!(uold::Array{Float64},z::Array{Float64},cv::Array{Float64},h::Array{Float64},unew::Array{Float64},tdts8,tdtsdx::Float64)
        #unew[2:end,1:end-1]= uold[2:end,1:end-1] .+ tdts8 .* 
        #                        (z[2:end, 2:end] .+ z[2:end, 1:end-1]) .*
        #                        (cv[2:end,2:end] .+ cv[2:end,1:end-1] .+ 
        #                        cv[1:end-1,2:end] .+ cv[1:end-1,1:end-1]) .-
        #                        tdtsdx .* (h[2:end, 1:end-1] .- h[1:end-1,1:end-1])
        @turbo for j=2:N+1,i=1:M
          @inbounds unew[i,j] = uold[i,j] + tdts8 * (z[i+1,j] + z[i+1,j-1]) * 
				(cv[i+1,j] + cv[i+1,j-1] + cv[i,j] + cv[i,j-1]) -
				tdtsdx * (h[i+1,j-1]-h[i,j-1])				
	end
end
function calc_vnew!(vold::Array{Float64},z::Array{Float64},cu::Array{Float64},h::Array{Float64},vnew::Array{Float64},tdts8::Float64,tdtsdy::Float64)
        # vnew[1:end-1,2:end] = vold[1:end-1,2:end] .- tdts8 .*
        #                         (z[2:end, 2:end] .+ z[1:end-1, 2:end]) .*
        #                        (cu[2:end,2:end] .+ cu[2:end,1:end-1] .+ 
        #                        cu[1:end-1,2:end] .+ cu[1:end-1,1:end-1]) .-
        #                        tdtsdy .* (h[1:end-1, 2:end] .- h[1:end-1,1:end-1])
	@turbo for j=2:N+1,i=1:M
	   @inbounds vnew[i,j] = vold[i,j] - tdts8 * (z[i+1,j] + z[i,j]) *
                                 (cu[i+1,j] + cu[i+1,j-1] + cu[i,j] + cu[i,j-1]) -
				 tdtsdy * (h[i,j] - h[i,j-1])
	end
end
function calc_pnew!(pold::Array{Float64},cu::Array{Float64},cv::Array{Float64},pnew::Array{Float64},tdtsdx::Float64)
        # pnew[1:end-1,1:end-1] = pold[1:end-1,1:end-1] .- tdtsdx .*
        #                         (cu[2:end, 1:end-1] .- cu[1:end-1,1:end-1]) .-
        #                         tdtsdy .* (cv[1:end-1, 2:end] .- cv[1:end-1,1:end-1])
	@turbo for j=1:N,i=1:M
 	   @inbounds pnew[i,j] = pold[i,j] - tdtsdx * (cu[i+2,j] - cu[i,j]) -
				 tdtsdy * ( cv[i,j+1] - cv[i,j])
        end
end
function bndry1!(h::Array{Float64},z::Array{Float64},cu::Array{Float64},cv::Array{Float64})
    @inbounds cu[1,1]  = cu[M, 1]
    @inbounds h[M,1]   = h[1, 1]
    @turbo for i=2:N+1
       @inbounds cu[1,i]  = cu[M, i]
       @inbounds h[M,i]   = h[1, i]
       @inbounds cv[M, i] = cv[1, i]
       @inbounds z[1, i]  = z[M, i]
    end
end

#@views cenU(A) = A[2:end,1:end-1]
#@views cenV(A) = A[1:end-1,2:end]
#@views cenP(A) = A[1:end-1,1:end-1]
#@views cenZ(A) = A[2:end,2:end]

#@views im1U(A) = A[1:end-1,1:end-1]
#@views jm1V(A) = A[1:end-1,1:end-1]
#@views allJ(n,A) = A[n,1:end]
#@views allI(n,A) = A[1:end,n]


global tdt = dt
@time begin
for ncycle in 1:ITMAX
    global alpha, fsdx, fsdy
    global p, u, v, cu, cv, z, h, uold, vold, pold, unew, vnew, pnew, tdt, tdtsdx, tdtsdy, time, tdts
    
     #
     # cu[2:end, 1:end-1] = ((p[2:end, 1:end-1] .+ p[1:end-1,1:end-1]) .* u[2:end,1:end-1]) ./ 2.0::Float64
     calc_cu!(p,u,cu)
     #
     # cv[1:end-1, 2:end] = ((p[1:end-1, 2:end] .+ p[1:end-1,1:end-1]) .* v[1:end-1,2:end]) ./ 2.0::Float64
     calc_cv!(p,v,cv)
     #
     # end



     # z[2:end, 2:end] = (fsdx .* (v[2:end, 2:end] .- v[1:end-1, 2:end]) .- 
     #                 fsdy .* (u[2:end, 2:end] .- u[2:end, 1:end-1])) ./
     #                 (p[1:end-1,1:end-1] .+ p[2:end,1:end-1] .+ 
     #                 p[2:end,2:end] .+ p[1:end-1,2:end])
     calc_z!(v,u,p,z,fsdx,fsdy)

    
     # h[1:end-1,1:end-1] = p[1:end-1,1:end-1] .+ 0.25::Float64 .* 
     #                        (u[2:end, 1:end-1] .* u[2:end, 1:end-1] .+
     #                        u[1:end-1,1:end-1] .* u[1:end-1,1:end-1] .+ 
     #                        v[1:end-1, 2:end] .* v[1:end-1, 2:end] .+ 
     #                       v[1:end-1,1:end-1] .* v[1:end-1,1:end-1])
     calc_h!(p,u,v,h)
    

    # bndry1!(h,z,cu,cv)
    cu[1,1:N+1]  = cu[M, 1:N+1]
     h[M,1:N+1]  =  h[1, 1:N+1]
    cv[M, 2:N+1] = cv[1, 2:N+1]
     z[M, 2:N+1] =  z[1, 2:N+1]

    cv[1:M+1, 1] = cv[1:M+1, N]
     h[1:M+1, N] =  h[1:M+1, 1]
    cu[2:M+1, N] = cu[2:M+1, 1]
     z[2:M+1, 1] =  z[2:M+1, N]

    cu[1, N]     = cu[M, 1]
    cv[M, 1]     = cv[1, N]
     z[1, 1]     =  z[M, N]
     h[M, N]     =  h[1, 1]

    tdts8 = tdt / 8.0::Float64
    tdtsdx = tdt / dx
    tdtsdy = tdt / dy

    #unew[2:end,1:end-1]= uold[2:end,1:end-1] .+ tdts8 .* 
    #                        (z[2:end, 2:end] .+ z[2:end, 1:end-1]) .*
    #                        (cv[2:end,2:end] .+ cv[2:end,1:end-1] .+ 
    #                        cv[1:end-1,2:end] .+ cv[1:end-1,1:end-1]) .-
    #                        tdtsdx .* (h[2:end, 1:end-1] .- h[1:end-1,1:end-1])
    calc_unew!(uold,z,cv,h,unew,tdts8,tdtsdx)

    # vnew[1:end-1,2:end] = vold[1:end-1,2:end] .- tdts8 .*
    #                         (z[2:end, 2:end] .+ z[1:end-1, 2:end]) .*
    #                        (cu[2:end,2:end] .+ cu[2:end,1:end-1] .+ 
    #                        cu[1:end-1,2:end] .+ cu[1:end-1,1:end-1]) .-
    #                        tdtsdy .* (h[1:end-1, 2:end] .- h[1:end-1,1:end-1])
    calc_vnew!(vold,z,cu,h,vnew,tdts8,tdtsdy)

    # pnew[1:end-1,1:end-1] = pold[1:end-1,1:end-1] .- tdtsdx .*
    #                         (cu[2:end, 1:end-1] .- cu[1:end-1,1:end-1]) .-
    #                         tdtsdy .* (cv[1:end-1, 2:end] .- cv[1:end-1,1:end-1])
    calc_pnew!(pold,cu,cv,pnew,tdtsdx)


    # @time begin
    unew[1, 1:end] = unew[M, 1:end]
    pnew[M, 1:end] = pnew[1, 1:end]
    vnew[M, 2:end] = vnew[1, 2:end]

    unew[2:end, N] = unew[2:end, 1]
    vnew[1:end, 1] = vnew[1:end, N]
    pnew[1:end, N] = pnew[1:end, 1]

    unew[1, N]     = unew[M, 1]
    vnew[M, 1]     = vnew[1, N]
    pnew[M, N]     = pnew[1, 1]
    # end

    if ncycle>1
        uold = u .+ alpha .* (unew .- 2.0::Float64 .* u .+ uold)
        vold = v .+ alpha .* (vnew .- 2.0::Float64 .* v .+ vold)
        pold = p .+ alpha .* (pnew .- 2.0::Float64 .* p .+ pold)

        u = unew
        v = vnew
        p = pnew
    else
        tdt = tdt + tdt
        uold = u
        vold = v
        pold = p
        u    = unew
        v    = vnew
        p    = pnew
    end

end
end

if L_OUT
    println(" Cycle number: ", ITMAX)
    println(" u:\n", diag(u)[1:end-1])
    println(" v:\n", diag(v)[1:end-1])
    println(" p:\n", diag(p)[1:end-1])
end
