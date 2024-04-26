using Profile
using Printf
using LinearAlgebra
using LoopVectorization

const M::Int = 16
const N::Int = 16
const M_LEN::Int = M + 1
const N_LEN::Int = N + 1
const ITMAX::Int = 4

const dt::Float64 = 90.::Float64
const dx::Float64 = 100000.::Float64
const dy::Float64 = 100000.::Float64
const fsdx::Float64 = 4.::Float64 / (dx)
const fsdy::Float64 = 4.::Float64 / (dy)
const a::Float64 = 1000000.::Float64
const alpha::Float64 = 0.001::Float64
L_OUT = false

function meminfo_julia()
  # @printf "GC total:  %9.3f MiB\n" Base.gc_total_bytes(Base.gc_num())/2^20
  # Total bytes (above) usually underreports, thus I suggest using live bytes (below)
  @printf "GC live:   %9.3f MiB\n" Base.gc_live_bytes()/2^20
  @printf "JIT:       %9.3f MiB\n" Base.jit_total_bytes()/2^20
  @printf "Max. RSS:  %9.3f MiB\n" Sys.maxrss()/2^20
end


function validate_uvp(u::Array{Float64}, v::Array{Float64}, p::Array{Float64}, M::Int, N::Int, step::Int, suffix::String)
    u_ref, v_ref, p_ref = read_uvp(step, suffix, M, N)
    # a = norm((u.-u_ref),Inf)
    println("U_ref - U: ",norm((u.-u_ref),Inf))
    println("U_ref: ",u_ref[1,1]);
    println("U[1,1]: ",u[1,1]);
    println("U[1,N]: ",u[1,N]);
    println("U[M,1]: ",u[M,1]);
    println("U[M,N]: ",u[M,N]);

    #println("V_ref - V: ",norm((v.-v_ref),Inf))
    #println("P_ref - P: ",norm((p.-p_ref),Inf))
    # np.testing.assert_allclose(u, u_ref)
    # np.testing.assert_allclose(v, v_ref)
    # np.testing.assert_allclose(p, p_ref)
    # print("step {step} {suffix} values are correct.\n")
end

function read_uvp(step::Int, suffix::String, M::Int, N::Int)
   
    base = "../ref/" * string(M) * "x" * string(N) * "/"
    step = "step" * string(step) * "." * suffix * ".bin"
    u_file = base * "u." * step;
    v_file = base * "v." * step;
    p_file = base * "p." * step;

    # u_file = f'../ref/{M}x{N}/u.step{step}.{suffix}.bin'
    # v_file = f'../ref/{M}x{N}/v.step{step}.{suffix}.bin'
    # p_file = f'../ref/{M}x{N}/p.step{step}.{suffix}.bin'

    u = zeros(Float64,M+1,M+1)
    v = zeros(Float64,M+1,M+1)
    p = zeros(Float64,M+1,M+1)

    read!(u_file,u)
    read!(v_file,v)
    read!(p_file,p)
    return u, v, p
end


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
function calc_unew!(uold::Array{Float64},z::Array{Float64},cv::Array{Float64},h::Array{Float64},unew::Array{Float64},tdts8::Float64,tdtsdx::Float64)
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
function calc_uvpnew!(uold::Array{Float64},
		     vold::Array{Float64},
		     pold::Array{Float64},
		       cu::Array{Float64},
		       cv::Array{Float64},
			z::Array{Float64},
			h::Array{Float64},
	     	     unew::Array{Float64},
	     	     vnew::Array{Float64},
	     	     pnew::Array{Float64},
		    tdts8::Float64,tdtsdy::Float64)
   @turbo for j=2:N+1,i=1:M
      @inbounds unew[i,j] = uold[i,j] + tdts8 * (z[i+1,j] + z[i+1,j-1]) * 
			(cv[i+1,j] + cv[i+1,j-1] + cv[i,j] + cv[i,j-1]) -
			tdtsdx * (h[i+1,j-1]-h[i,j-1])				
      @inbounds vnew[i,j] = vold[i,j] - tdts8 * (z[i+1,j] + z[i,j]) *
                           (cu[i+1,j] + cu[i+1,j-1] + cu[i,j] + cu[i,j-1]) -
			 tdtsdy * (h[i,j] - h[i,j-1])
   end
   @turbo for j=1:N,i=1:M
      @inbounds pnew[i,j] = pold[i,j] - tdtsdx * (cu[i+2,j] - cu[i,j]) -
			 tdtsdy * ( cv[i,j+1] - cv[i,j])
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

@views d_yi(A) = A[2:end, 2:end] .- A[2:end, 1:end-1]
@views d_xi(A) = A[2:end, 2:end] .- A[1:end-1, 2:end]
@views dp(A) = p[1:end-1,1:end-1] .+ A[2:end,1:end-1] .+ A[2:end,2:end] .+ A[1:end-1,2:end]

global tdt = dt
@time begin
for ncycle in 1:ITMAX
   global alpha, fsdx, fsdy,tdts8
   global p, u, v, cu, cv, z, h, uold, vold, pold, unew, vnew, pnew, tdt, tdtsdx, tdtsdy, time, tdts
    

   # read_uvp(step::Int, suffix::String, M::Int, N::Int):
   if ncycle == 1
       # foo,bar,me = read_uvp(ncycle,"init",M,N)
       validate_uvp(u,v,p,M,N,ncycle,"init")
   end
   # if VAL_DEEP and ncycle <= 3
   #    validate_uvp(u,v,p,M,N,ncycle,'init')
   calc_cu!(p,u,cu)
   calc_cv!(p,v,cv)



   #  z = ((fsdx .* d_xi(v)) .- (fsdy .* d_yi(u)))./dp(p)
   # z[2:end, 2:end] = (fsdx .* (v[2:end, 2:end] .- v[1:end-1, 2:end]) .- 
   #                  fsdy .* (u[2:end, 2:end] .- u[2:end, 1:end-1])) ./
   #                  (p[1:end-1,1:end-1] .+ p[2:end,1:end-1] .+ 
   #                 p[2:end,2:end] .+ p[1:end-1,2:end])
   calc_z!(v,u,p,z,fsdx,fsdy)

   calc_h!(p,u,v,h)
    
   cu[1,1:N+1]  = cu[M, 1:N+1]
    h[M,1:N+1]  =  h[1, 1:N+1]
   cv[M, 2:N+1] = cv[1, 2:N+1]
    # z[M, 2:N+1] =  z[1, 2:N+1]

   cv[1:M+1, 1] = cv[1:M+1, N]
    h[1:M+1, N] =  h[1:M+1, 1]
   cu[2:M+1, N] = cu[2:M+1, 1]
    # z[2:M+1, 1] =  z[2:M+1, N]

   cu[1, N]     = cu[M, 1]
   cv[M, 1]     = cv[1, N]
    # z[1, 1]     =  z[M, N]
    h[M, N]     =  h[1, 1]

   tdts8 = tdt / 8.0::Float64
   tdtsdx = tdt / dx
   tdtsdy = tdt / dy

   # calc_unew!(uold,z,cv,h,unew,tdts8,tdtsdx)
   # calc_vnew!(vold,z,cu,h,vnew,tdts8,tdtsdy)
   # calc_pnew!(pold,cu,cv,pnew,tdtsdx)
   calc_uvpnew!(uold,vold,pold,cu,cv,z,h,unew,vnew,pnew,tdts8,tdtsdy)


   unew[1, 1:N+1] = unew[M, 1:N+1]
   pnew[M, 1:N+1] = pnew[1, 1:N+1]
   vnew[M, 2:N+1] = vnew[1, 2:N+1]

   unew[2:M+1, N] = unew[2:M+1, 1]
   vnew[1:M+1, 1] = vnew[1:M+1, N]
   pnew[1:M+1, N] = pnew[1:M+1, 1]

   unew[1, N]     = unew[M, 1]
   vnew[M, 1]     = vnew[1, N]
   pnew[M, N]     = pnew[1, 1]

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

   if mod(ncycle,10) == 0
      println(" NCYCLE is ",ncycle)
      meminfo_julia()
   end
end
end

if L_OUT
    println(" Cycle number: ", ITMAX)
    println(" u:\n", diag(u)[1:end-1])
    println(" v:\n", diag(v)[1:end-1])
    println(" p:\n", diag(p)[1:end-1])
end
