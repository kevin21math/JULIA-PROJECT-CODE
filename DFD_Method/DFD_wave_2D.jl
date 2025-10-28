using SparseArrays, BSplineKit, LinearAlgebra, QuadGK, Plots,FastGaussQuadrature
include("C:\\Users\\hp\\Downloads\\VS CODE PRG\\julia project\\DFD_Method\\galerkin_cross.jl")
include("C:\\Users\\hp\\Downloads\\VS CODE PRG\\julia project\\DFD_Method\\shifted_knots.jl")

function to_hat(U::Matrix{Float64}, LX, LY)  
    return LX\U/LY'
end

function from_hat(U::Matrix{Float64}, LX, LY)  
    return LX'\U/LY
end

function project_point_source(xbasis, ybasis,LX1, LY1, x_s, y_s, amplitude)
    nbx = length(xbasis); nby = length(ybasis)
    bx = [evaluate(xbasis, j, x_s) for j in 1:nbx]
    by = [evaluate(ybasis, j, y_s) for j in 1:nby]
    S = bx * by'                                 
    F = amplitude* S
    F_hat=to_hat(F,LX1,LY1)              
    return F_hat
end

custom_colors = cgrad([ RGB(0, 0, 1),RGB(0, 1, 0),RGB(1, 1, 0),RGB(1, 0, 0)], rev = false)

vp=3400
vs=1963
dt=1.18*10^(-7)
rho=2500
λ = rho * (vp^2 - 2*vs^2)             
μ = rho * vs^2
T=1.4*10^(-4)
n=Int(round(T/dt))

a, b = 0.0, 1.0
order = 4
nx_dom=10
ny_dom=10
nx_local = 26
ny_local = 26
Nx = nx_dom * (nx_local - 1) + 1
Ny = ny_dom * (ny_local - 1) + 1
x  = LinRange(a,b,Nx)
y  = LinRange(a,b,Ny)

subdomain_width = (b - a) /  nx_dom
subdomain_height = (b - a) / ny_dom

subdomains = []
for j in 0:(ny_dom - 1)  
    for i in 0:(nx_dom - 1) 
        
        x_min = a + i * subdomain_width
        x_max = x_min + subdomain_width
        y_min = a + j * subdomain_height
        y_max = y_min + subdomain_height
        
        domain = (
            id = j * nx_dom + i + 1,
            ix = i + 1, iy = j + 1, 
            x_min = x_min,x_max = x_max,
            y_min = y_min,y_max = y_max
            )

        push!(subdomains, domain)
    end
end
ndomains=length(subdomains)

function get_neighbors(ix, iy, Nx, Ny)
    left   = (ix > 1)  ? (iy - 1) * Nx + (ix - 1) : nothing
    right  = (ix < Nx) ? (iy - 1) * Nx + (ix + 1) : nothing
    bottom = (iy > 1)  ? (iy - 2) * Nx + ix       : nothing
    top    = (iy < Ny) ? iy * Nx + ix             : nothing
    return (left, right, bottom, top)
end

D1 = BSplineKit.Derivative(1)
D0 = BSplineKit.Derivative(0)
alphaxp = fill(0.5, nx_dom)
alphayp = fill(0.5, ny_dom)
alphaxm = fill(0.5, nx_dom)
alphaym = fill(0.5, ny_dom)
alphaxp[end] = 0.0 
alphayp[end] = 0.0
alphaxm[1] = 1.0 
alphaym[1] = 1.0

c1, d1 = -0.2,0
v1 = c1 .+ (d1 - c1) .* rand(ndomains)
c2, d2 = 0, 0.2 
v2 = c2 .+ (d2 - c2) .* rand(ndomains)

x1basis        = Vector{BSplineBasis}(undef, ndomains)
x2basis        = Vector{BSplineBasis}(undef, ndomains)
y1basis        = Vector{BSplineBasis}(undef, ndomains)
y2basis        = Vector{BSplineBasis}(undef, ndomains)
LX1            = Vector{Matrix{Float64}}(undef, ndomains)
LX2            = Vector{Matrix{Float64}}(undef, ndomains)
LY1            = Vector{Matrix{Float64}}(undef, ndomains)
LY2            = Vector{Matrix{Float64}}(undef, ndomains)
Dx0            = Vector{Matrix{Float64}}(undef, ndomains) 
Dy0            = Vector{Matrix{Float64}}(undef, ndomains)
Kx             = Vector{Matrix{Float64}}(undef, ndomains)
Ky             = Vector{Matrix{Float64}}(undef, ndomains)
qx1plus        = Vector{Vector{Float64}}(undef, ndomains)
qx1minus       = Vector{Vector{Float64}}(undef, ndomains)
qx2plus        = Vector{Vector{Float64}}(undef, ndomains)
qx2minus       = Vector{Vector{Float64}}(undef, ndomains)
qy1plus        = Vector{Vector{Float64}}(undef, ndomains)
qy1minus       = Vector{Vector{Float64}}(undef, ndomains)
qy2plus        = Vector{Vector{Float64}}(undef, ndomains)
qy2minus       = Vector{Vector{Float64}}(undef, ndomains)
bx1plus        = Vector{Vector{Float64}}(undef, ndomains)
bx2plus        = Vector{Vector{Float64}}(undef, ndomains)
bx1minus       = Vector{Vector{Float64}}(undef, ndomains)
bx2minus       = Vector{Vector{Float64}}(undef, ndomains)
by1plus        = Vector{Vector{Float64}}(undef, ndomains)
by2plus        = Vector{Vector{Float64}}(undef, ndomains) 
by1minus       = Vector{Vector{Float64}}(undef, ndomains)
by2minus       = Vector{Vector{Float64}}(undef, ndomains)
Dx1bottom      = Vector{Matrix{Float64}}(undef, ndomains)
Dx2bottom      = Vector{Matrix{Float64}}(undef, ndomains)
Dx1top         = Vector{Matrix{Float64}}(undef, ndomains)
Dx2top         = Vector{Matrix{Float64}}(undef, ndomains)
Dy1left        = Vector{Matrix{Float64}}(undef, ndomains)
Dy2left        = Vector{Matrix{Float64}}(undef, ndomains)
Dy1right       = Vector{Matrix{Float64}}(undef, ndomains)
Dy2right       = Vector{Matrix{Float64}}(undef, ndomains)
Dx1            = Vector{Matrix{Float64}}(undef, ndomains)
Dy1            = Vector{Matrix{Float64}}(undef, ndomains)
Ux             = Vector{Matrix{Float64}}(undef, ndomains)
Uy             = Vector{Matrix{Float64}}(undef, ndomains)
Uxhat          = Vector{Matrix{Float64}}(undef, ndomains)
Uyhat          = Vector{Matrix{Float64}}(undef, ndomains)
Uxhat1         = Vector{Matrix{Float64}}(undef, ndomains)
Uyhat1         = Vector{Matrix{Float64}}(undef, ndomains)
Uxhat_new      = Vector{Matrix{Float64}}(undef, ndomains)
Uyhat_new      = Vector{Matrix{Float64}}(undef, ndomains)
Uxhat_old      = Vector{Matrix{Float64}}(undef, ndomains)
Uyhat_old      = Vector{Matrix{Float64}}(undef, ndomains)
Sxx            = Vector{Matrix{Float64}}(undef, ndomains)
Sxy            = Vector{Matrix{Float64}}(undef, ndomains)
Syx            = Vector{Matrix{Float64}}(undef, ndomains)
Syy            = Vector{Matrix{Float64}}(undef, ndomains)
F_hat_x        = Vector{Matrix{Float64}}(undef, ndomains)
F_hat_y        = Vector{Matrix{Float64}}(undef, ndomains)
ux_out_minus   = Vector{Matrix{Float64}}(undef, ndomains)
ux_out_plus    = Vector{Matrix{Float64}}(undef, ndomains)
uy_out_minus   = Vector{Matrix{Float64}}(undef, ndomains)
uy_out_plus    = Vector{Matrix{Float64}}(undef, ndomains)
ux_out_bottom  = Vector{Matrix{Float64}}(undef, ndomains)
ux_out_top     = Vector{Matrix{Float64}}(undef, ndomains)
uy_out_bottom  = Vector{Matrix{Float64}}(undef, ndomains)
uy_out_top     = Vector{Matrix{Float64}}(undef, ndomains)
sxx_out_minus  = Vector{Matrix{Float64}}(undef, ndomains)
sxx_out_plus   = Vector{Matrix{Float64}}(undef, ndomains)
syx_out_minus  = Vector{Matrix{Float64}}(undef, ndomains)
syx_out_plus   = Vector{Matrix{Float64}}(undef, ndomains)
sxy_out_bottom = Vector{Matrix{Float64}}(undef, ndomains)
sxy_out_top    = Vector{Matrix{Float64}}(undef, ndomains)
syy_out_bottom = Vector{Matrix{Float64}}(undef, ndomains)
syy_out_top    = Vector{Matrix{Float64}}(undef, ndomains)
ux_in_left     = Vector{Matrix{Float64}}(undef, ndomains)
ux_in_right    = Vector{Matrix{Float64}}(undef, ndomains)
uy_in_left     = Vector{Matrix{Float64}}(undef, ndomains)
uy_in_right    = Vector{Matrix{Float64}}(undef, ndomains)
ux_in_bottom   = Vector{Matrix{Float64}}(undef, ndomains)
ux_in_top      = Vector{Matrix{Float64}}(undef, ndomains)
uy_in_bottom   = Vector{Matrix{Float64}}(undef, ndomains)
uy_in_top      = Vector{Matrix{Float64}}(undef, ndomains)
sxx_in_left    = Vector{Matrix{Float64}}(undef, ndomains)
sxx_in_right   = Vector{Matrix{Float64}}(undef, ndomains)
syx_in_left    = Vector{Matrix{Float64}}(undef, ndomains)
syx_in_right   = Vector{Matrix{Float64}}(undef, ndomains)
sxy_in_bottom  = Vector{Matrix{Float64}}(undef, ndomains)
sxy_in_top     = Vector{Matrix{Float64}}(undef, ndomains)
syy_in_bottom  = Vector{Matrix{Float64}}(undef, ndomains)
syy_in_top     = Vector{Matrix{Float64}}(undef, ndomains)
N_ux           = Vector{Int}(undef, ndomains)     
N_sigmax       = Vector{Int}(undef, ndomains)
N_uy           = Vector{Int}(undef, ndomains)     
N_sigmay       = Vector{Int}(undef, ndomains)
Ux_total       = zeros(Nx, Ny)

f0 = 60000.0  
t0 = 1.0 / f0
function ricker_wavelet(t::Float64, f0::Float64, t0::Float64)
    arg = π * f0 * (t - t0)
    return (1.0 - 2.0 * arg^2) * exp(-arg^2)
end
x_s_global = 0.5; y_s_global = 0.5
source_ix = round(Int, Nx / 2)
source_iy = round(Int, Ny / 2) 

for j=1:ndomains
   
    current_domain = subdomains[j]
    x_left=current_domain.x_min
    x_right=current_domain.x_max
    y_left=current_domain.y_min
    y_right=current_domain.y_max

    x1knots = shifted_knots(x_left, x_right, nx_local,order,v1[j])
    x1basis[j] = BSplineBasis(order, x1knots[order:end-order+1])
    x2knots = shifted_knots(x_left, x_right, nx_local,order,v2[j])
    x2basis[j] = BSplineBasis(order, x2knots[order:end-order+1])
    y1knots = shifted_knots(y_left,  y_right, ny_local,order,v1[j])
    y1basis[j] = BSplineBasis(order, y1knots[order:end-order+1])
    y2knots = shifted_knots(y_left,  y_right, ny_local,order,v2[j])
    y2basis[j] = BSplineBasis(order, y2knots[order:end-order+1])
  
    MX1 = galerkin_matrix(x1basis[j], (D0,D0))
    LX1[j] = cholesky(Hermitian(MX1)).L
    MX2 = galerkin_matrix(x2basis[j], (D0,D0)) 
    LX2[j] = cholesky(Hermitian(MX2)).L
    MY1 = galerkin_matrix(y1basis[j], (D0,D0))
    LY1[j] = cholesky(Hermitian(MY1)).L
    MY2 = galerkin_matrix(y2basis[j], (D0,D0))
    LY2[j] = cholesky(Hermitian(MY2)).L

    N_ux[j] = length(x1basis[j])
    N_sigmax[j] =length(x2basis[j])
    N_uy[j] = length(y1basis[j])
    N_sigmay[j] =length(y2basis[j])

    qx1plus[j]  = [evaluate(x1basis[j], i, x_right, D0) for i in 1:N_ux[j]]
    qx1minus[j] = [evaluate(x1basis[j], i, x_left,  D0) for i in 1:N_ux[j]]
    qx2plus[j]  = [evaluate(x2basis[j], i, x_right, D0) for i in 1:N_sigmax[j]]
    qx2minus[j] = [evaluate(x2basis[j], i, x_left,  D0) for i in 1:N_sigmax[j]]
    qy1plus[j]  = [evaluate(y1basis[j], i, y_right, D0) for i in 1:N_uy[j]]
    qy1minus[j] = [evaluate(y1basis[j], i, y_left,  D0) for i in 1:N_uy[j]]
    qy2plus[j]  = [evaluate(y2basis[j], i, y_right, D0) for i in 1:N_sigmay[j]]
    qy2minus[j] = [evaluate(y2basis[j], i, y_left,  D0) for i in 1:N_sigmay[j]]
end
  
for j=1:ndomains

    dom = subdomains[j]
    (left, right, bottom, top) = get_neighbors(dom.ix, dom.iy, nx_dom, ny_dom)
    idx = dom.ix
    idy = dom.iy

    bx1minus[j]=(-1)*(1-alphaxm[idx])*qx1minus[j]
    bx2minus[j]=-alphaxm[idx]*qx2minus[j]
    by1minus[j]=-(1-alphaym[idy])*qy1minus[j]
    by2minus[j]=-alphaym[idy]*qy2minus[j]
    bx1plus[j] = (1 - alphaxp[idx]) * qx1plus[j]
    bx2plus[j] = (alphaxp[idx]) * qx2plus[j]
    by1plus[j] = (1 - alphayp[idy]) * qy1plus[j]
    by2plus[j] = (alphayp[idy]) * qy2plus[j]
    
    Dx0[j] = galerkin_cross(x2basis[j], x1basis[j],(D0,D0); nquad = 24)
    Dy0[j] = galerkin_cross(y2basis[j], y1basis[j],(D0,D0);nquad = 24)
    
    if bottom !== nothing
       Dx1bottom[j] = galerkin_cross(x1basis[j], x2basis[bottom], (D0,D0); nquad = 24)
       Dx2bottom[j] = galerkin_cross(x2basis[j], x1basis[bottom],(D0,D0); nquad = 24)
    else
       Dx1bottom[j] = zeros(length(x1basis[j]), length(x2basis[j]))
       Dx2bottom[j] = zeros(length(x2basis[j]), length(x1basis[j]))  
    end

    if top !== nothing
       Dx1top[j] = galerkin_cross(x1basis[j], x2basis[top], (D0,D0); nquad = 24)
       Dx2top[j] = galerkin_cross(x2basis[j], x1basis[top],(D0,D0); nquad = 24)
    else
       Dx1top[j] = zeros(length(x1basis[j]), length(x2basis[j]))
       Dx2top[j] = zeros(length(x2basis[j]), length(x1basis[j]))
    end

    if left !== nothing
       Dy1left[j] = galerkin_cross(y1basis[j], y2basis[left], (D0,D0); nquad = 24)
       Dy2left[j] = galerkin_cross(y2basis[j], y1basis[left],(D0,D0); nquad = 24)
    else
       Dy1left[j] = zeros(length(y1basis[j]), length(y2basis[j]))
       Dy2left[j] = zeros(length(y2basis[j]), length(y1basis[j]))
    end

    if right !== nothing
        Dy1right[j] = galerkin_cross(y1basis[j], y2basis[right], (D0,D0); nquad = 24)
        Dy2right[j] = galerkin_cross(y2basis[j], y1basis[right],(D0,D0); nquad = 24)
    else
        Dy1right[j] = zeros(length(y1basis[j]), length(y2basis[j]))
        Dy2right[j] = zeros(length(y2basis[j]), length(y1basis[j]))
    end
    
    Kx[j] = galerkin_cross(x2basis[j], x1basis[j],(D0,D1);nquad = 24)
    Ky[j] = galerkin_cross(y2basis[j], y1basis[j],(D0,D1);nquad = 24)

    Qx_plus = (+1.0) .* (qx2plus[j] * qx1plus[j]')
    Qx_minus = (-1.0) .* (qx2minus[j] * qx1minus[j]')
    Qy_plus =  (+1.0) .* (qy2plus[j] * qy1plus[j]')
    Qy_minus =(-1.0) .* (qy2minus[j] * qy1minus[j]')

    Dx1[j] = Kx[j] - alphaxm[idx] * Qx_minus - alphaxp[idx] * Qx_plus
    Dy1[j] = Ky[j] - alphaym[idy] * Qy_minus - alphayp[idy] * Qy_plus
    
    Ux[j] = zeros(N_ux[j], N_uy[j])
    Uy[j] = zeros(N_ux[j], N_uy[j])
    Uxhat[j]=to_hat(Ux[j],LX1[j],LY1[j])
    Uyhat[j]=to_hat(Uy[j],LX1[j],LY1[j])

    source_time_function = ricker_wavelet(dt,f0,t0)
    F_hat_y[j] = spzeros(N_ux[j], N_uy[j])
    if (x_s_global >= dom.x_min && x_s_global <= dom.x_max && y_s_global >= dom.y_min && y_s_global <= dom.y_max)
        F_hat_x[j] = project_point_source(x1basis[j], y1basis[j], LX1[j], LY1[j], x_s_global, y_s_global, source_time_function)
    else
        F_hat_x[j] = spzeros(N_ux[j], N_uy[j])
    end
    Uxhat1[j]=(dt^2/(2*rho)).*(F_hat_x[j])
    Uyhat1[j]= (dt^2/(2*rho)).*(F_hat_y[j])

    Uxhat_old[j]=Uxhat[j]
    Uxhat[j]=Uxhat1[j]
    Uyhat_old[j]=Uyhat[j]
    Uyhat[j]=Uyhat1[j]

    Ux[j]=from_hat(Uxhat[j],LX1[j],LY1[j])
    Uy[j]=from_hat(Uyhat[j],LX1[j],LY1[j])

    Sxx[j] = zeros(N_sigmax[j], N_sigmay[j])
    Sxy[j] = zeros(N_sigmax[j], N_sigmay[j])
    Syx[j] = zeros(N_sigmax[j], N_sigmay[j])
    Syy[j] = zeros(N_sigmax[j], N_sigmay[j])
end

anim = Animation()
frame_rate = 20

for i=1:n
    for j in 1:ndomains
      Iux = I(N_ux[j]);  Iuy = I(N_uy[j])

      ux_out_minus[j] = qx1minus[j]' * Ux[j] * Iuy        
      ux_out_plus[j]  = qx1plus[j]'  * Ux[j] * Iuy       
      uy_out_minus[j] = qx1minus[j]' * Uy[j] * Iuy
      uy_out_plus[j]  = qx1plus[j]'  * Uy[j] * Iuy
      ux_out_bottom[j] = Iux * Ux[j] * reshape(qy1minus[j], :, 1)   
      ux_out_top[j]    = Iux * Ux[j] * reshape(qy1plus[j], :, 1)   
      uy_out_bottom[j] = Iux * Uy[j] * reshape(qy1minus[j], :, 1)
      uy_out_top[j]    = Iux * Uy[j] * reshape(qy1plus[j], :, 1)
    end

    for j in 1:ndomains
      dom = subdomains[j]
      (left, right, bottom, top) = get_neighbors(dom.ix, dom.iy, nx_dom, ny_dom)

      ux_in_left[j]  = zeros(size(ux_out_plus[j]))
      ux_in_right[j] = zeros(size(ux_out_minus[j]))
      uy_in_left[j]  = zeros(size(uy_out_plus[j]))
      uy_in_right[j] = zeros(size(uy_out_minus[j])) 
      ux_in_bottom[j] = zeros(size(ux_out_top[j]))     
      ux_in_top[j]    = zeros(size(ux_out_bottom[j])) 
      uy_in_bottom[j] = zeros(size(uy_out_top[j]))  
      uy_in_top[j]    = zeros(size(uy_out_bottom[j])) 

      if left !== nothing
          ux_in_left[j]  = ux_out_plus[left]     
          uy_in_left[j]  = uy_out_plus[left]
      end
      if right !== nothing
          ux_in_right[j]  = ux_out_minus[right]
          uy_in_right[j]  = uy_out_minus[right]
      end
      if bottom !== nothing
          ux_in_bottom[j]  = ux_out_top[bottom]
          uy_in_bottom[j]  = uy_out_top[bottom]
      end
      if top !== nothing
          ux_in_top[j]  = ux_out_bottom[top]
          uy_in_top[j]  = uy_out_bottom[top]
      end
    end

    for j=1:ndomains
      domain = subdomains[j]
      source_time_function = ricker_wavelet(i*dt,f0,t0)
      F_hat_y[j] = spzeros(N_ux[j], N_uy[j])

      if (x_s_global >= domain.x_min && x_s_global <= domain.x_max &&
         y_s_global >= domain.y_min && y_s_global <= domain.y_max)
         F_hat_x[j] = project_point_source(x1basis[j], y1basis[j], LX1[j], LY1[j], x_s_global, y_s_global, source_time_function)
      else
         F_hat_x[j] = spzeros(N_ux[j], N_uy[j])
      end
    end

    for j=1:ndomains

      Ustarx10 = Dx1[j] * Ux[j] * Dy0[j]' + bx2minus[j] * ux_in_left[j]  * Dy2left[j]' + bx2plus[j]  * ux_in_right[j] * Dy2right[j]'
      
      Ustarx01 = Dx0[j] * Ux[j] * Dy1[j]' + Dx2bottom[j] * ux_in_bottom[j] * by2minus[j]' + Dx2top[j]    * ux_in_top[j]    * by2plus[j]'
      
      Ustary10 = Dx1[j] * Uy[j] * Dy0[j]' + bx2minus[j] * uy_in_left[j]  * Dy2left[j]' + bx2plus[j]  * uy_in_right[j] * Dy2right[j]'
       
      Ustary01 = Dx0[j] * Uy[j] * Dy1[j]' + Dx2bottom[j] * uy_in_bottom[j] * by2minus[j]' + Dx2top[j]    * uy_in_top[j]    * by2plus[j]'
      
      Uhat_x10=to_hat(Ustarx10,LX2[j],LY2[j])
      Uhat_x01=to_hat(Ustarx01,LX2[j],LY2[j])
      Uhat_y10=to_hat(Ustary10,LX2[j],LY2[j])
      Uhat_y01=to_hat(Ustary01,LX2[j],LY2[j])

      σxx_hat=(λ+2*μ)*(Uhat_x10)+λ*(Uhat_y01)
      σyy_hat=(λ)*(Uhat_x10)+(λ+2*μ)*(Uhat_y01)
      σxy_hat=(μ)*(Uhat_x01+Uhat_y10)

      Sxx[j]=from_hat(σxx_hat,LX2[j],LY2[j])
      Sxy[j]=from_hat(σxy_hat,LX2[j],LY2[j])
      Syx[j]=from_hat(σxy_hat,LX2[j],LY2[j])
      Syy[j]=from_hat(σyy_hat,LX2[j],LY2[j])
    end

    for j=1:ndomains
      ISx = I(N_sigmax[j]); ISy = I(N_sigmay[j])

      sxx_out_minus[j] = qx2minus[j]' * Sxx[j] * ISy
      sxx_out_plus[j]  = qx2plus[j]'  * Sxx[j] * ISy
      syx_out_minus[j] = qx2minus[j]' * Syx[j] * ISy
      syx_out_plus[j]  = qx2plus[j]'  * Syx[j] * ISy
      sxy_out_bottom[j] = ISx * Sxy[j] * reshape(qy2minus[j], :, 1)
      sxy_out_top[j]    = ISx * Sxy[j] * reshape(qy2plus[j], :, 1)
      syy_out_bottom[j] = ISx * Syy[j] * reshape(qy2minus[j], :, 1)
      syy_out_top[j]    = ISx * Syy[j] * reshape(qy2plus[j], :, 1)
    end

    for j in 1:ndomains
       dom = subdomains[j]
       (left, right, bottom, top) = get_neighbors(dom.ix, dom.iy, nx_dom, ny_dom)

       sxx_in_left[j] = zeros(size(sxx_out_plus[j]))
       sxx_in_right[j]= zeros(size(sxx_out_minus[j]))
       syx_in_left[j] = zeros(size(syx_out_plus[j]))
       syx_in_right[j]= zeros(size(syx_out_minus[j]))
       sxy_in_bottom[j]= zeros(size(sxy_out_top[j]))
       sxy_in_top[j]   = zeros(size(sxy_out_bottom[j]))
       syy_in_bottom[j]= zeros(size(syy_out_top[j]))
       syy_in_top[j]   = zeros(size(syy_out_bottom[j]))

       if left !== nothing
           sxx_in_left[j] = sxx_out_plus[left]
           syx_in_left[j] = syx_out_plus[left]
       end
       if right !== nothing
           sxx_in_right[j] = sxx_out_minus[right]
           syx_in_right[j] = syx_out_minus[right]
       end
       if bottom !== nothing
           sxy_in_bottom[j] = sxy_out_top[bottom]
           syy_in_bottom[j] = syy_out_top[bottom]
       end
       if top !== nothing
           sxy_in_top[j] = sxy_out_bottom[top]
           syy_in_top[j] = syy_out_bottom[top]
       end
     end

    for j=1:ndomains
      Sxx10_star = -(Dx1[j]' * Sxx[j] * Dy0[j]) + bx1minus[j] * sxx_in_left[j]  * Dy1left[j]' + bx1plus[j]  * sxx_in_right[j] * Dy1right[j]'

      Sxy01_star = -(Dx0[j]' * Sxy[j] * Dy1[j]) + Dx1bottom[j] * sxy_in_bottom[j] * by1minus[j]' + Dx1top[j]    * sxy_in_top[j]    * by1plus[j]'

      Syx10_star = -(Dx1[j]' * Syx[j] * Dy0[j]) + bx1minus[j] * syx_in_left[j]  * Dy1left[j]' + bx1plus[j]  * syx_in_right[j] * Dy1right[j]'

      Syy01_star = -(Dx0[j]' * Syy[j] * Dy1[j]) + Dx1bottom[j] * syy_in_bottom[j] * by1minus[j]' + Dx1top[j]    * syy_in_top[j]    * by1plus[j]'

      Sxx10_hat=to_hat(Sxx10_star,LX1[j],LY1[j])
      Syx10_hat=to_hat(Syx10_star,LX1[j],LY1[j])
      Sxy01_hat=to_hat( Sxy01_star,LX1[j],LY1[j])
      Syy01_hat=to_hat( Syy01_star,LX1[j],LY1[j])

      Uxhat_new[j]=2*Uxhat[j]-Uxhat_old[j]+(dt^2/rho)*(Sxx10_hat+Sxy01_hat+F_hat_x[j])
      Uyhat_new[j]=2*Uyhat[j]-Uyhat_old[j]+(dt^2/rho)*(Syx10_hat+Syy01_hat+F_hat_y[j])

      Uxhat_old[j]=Uxhat[j]
      Uxhat[j]=Uxhat_new[j]
      Uyhat_old[j]=Uyhat[j]
      Uyhat[j]=Uyhat_new[j]

      Ux[j]=from_hat(Uxhat[j],LX1[j],LY1[j])
      Uy[j]=from_hat(Uyhat[j],LX1[j],LY1[j])
    end

    if i % frame_rate == 0
      current_t = i * dt
      Ux_total .= 0.0

        for j in 1:ndomains
          dom = subdomains[j]
          x_local = LinRange(dom.x_min, dom.x_max, nx_local)
          y_local = LinRange(dom.y_min, dom.y_max, ny_local)
            
          B_x_vals = [evaluate(x1basis[j], i, val) for val in x_local, i in 1:N_ux[j]]
          B_y_vals = [evaluate(y1basis[j], i, val) for val in y_local, i in 1:N_uy[j]]
          Ux_grid = B_x_vals * Ux[j] * B_y_vals'

          ix_start = (dom.ix - 1) * (nx_local - 1) + 1
          ix_end   = ix_start + nx_local - 1
          iy_start = (dom.iy - 1) * (ny_local - 1) + 1
          iy_end   = iy_start + ny_local - 1

          Ux_total[ix_start:ix_end, iy_start:iy_end] = Ux_grid
        end

      max_abs_val=5*10^(-13)
      p = heatmap(x, y, Ux_total', aspect_ratio=:equal, c=custom_colors,clim=(-max_abs_val, max_abs_val), title="Horizontal Displacement (Ux) | t=$(round(current_t*1000000, sigdigits=2)) μs",
          xlabel="x (m)",ylabel="y (m)",size=(900, 700))
      scatter!(p, [x[source_ix]],[y[source_iy]], marker=:diamond, markersize=5,  markercolor=:red)
      println("Frame $(i/frame_rate)/$(n/frame_rate) at t=$(round(current_t*1000000, sigdigits=2)) μs captured.")
      savefig(p, "Final_time_plot.png")
      frame(anim, p)
    end
end
gif(anim, "Elastic_Wave_DFDM_2D.gif", fps=10)


