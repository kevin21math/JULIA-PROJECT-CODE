using SparseArrays, BSplineKit, LinearAlgebra, QuadGK, Plots,FastGaussQuadrature
include("C:\\Users\\hp\\Downloads\\VS CODE PRG\\julia project\\DFD_Method\\galerkin_cross.jl")
include("C:\\Users\\hp\\Downloads\\VS CODE PRG\\julia project\\DFD_Method\\shifted_knots.jl")

function to_hat(U::Matrix{Float64}, LX::LowerTriangular, LY::LowerTriangular)  
    return LX\U/LY'
end

function from_hat(U::Matrix{Float64}, LX::LowerTriangular, LY::LowerTriangular)  
    return LX'\U/LY
end

function project_point_source(xbasis, ybasis,LX1::LowerTriangular, LY1::LowerTriangular, x_s, y_s, amplitude)
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
Nx,Ny = 250,250
x  = LinRange(a,b,Nx)
y  = LinRange(a,b,Ny)
dx=(b-a)/(Nx-1)
dy=(b-a)/(Ny-1)
order = 4 

D1 = BSplineKit.Derivative(1)
D0 = BSplineKit.Derivative(0)
alphap=0
alpham=0

x1knots = shifted_knots(a,b,Nx,order,-0.125)
x1basis = BSplineBasis(order, x1knots[order:end-order+1])
x2knots = shifted_knots(a,b,Nx,order,0.125)
x2basis = BSplineBasis(order, x2knots[order:end-order+1])
y1knots = shifted_knots(a,b,Ny,order,-0.1)
y1basis = BSplineBasis(order, y1knots[order:end-order+1])
y2knots = shifted_knots(a,b,Ny,order,0.1)
y2basis = BSplineBasis(order, y2knots[order:end-order+1])

nbasis_x1, nbasis_x2 = length(x1basis), length(x2basis)
nbasis_y1, nbasis_y2 = length(y1basis), length(y2basis)
MX1 = galerkin_matrix(x1basis, (D0,D0))
LX1 = cholesky(Hermitian(MX1)).L
MX2 = galerkin_matrix(x2basis, (D0,D0)) 
LX2 = cholesky(Hermitian(MX2)).L
MY1 = galerkin_matrix(y1basis, (D0,D0))
LY1 = cholesky(Hermitian(MY1)).L
MY2 = galerkin_matrix(y2basis, (D0,D0))
LY2 = cholesky(Hermitian(MY2)).L

Dx0 = galerkin_cross(x2basis, x1basis,(D0,D0); nquad = 24)
Dy0 = galerkin_cross(y2basis, y1basis,(D0,D0);nquad = 24)

Kx = galerkin_cross(x2basis, x1basis,(D0,D1);nquad = 24)
Ky = galerkin_cross(y2basis, y1basis,(D0,D1);nquad = 24)

x1plus = [evaluate(x1basis, j, b, D0) for j in 1:nbasis_x1]
x1minus = [evaluate(x1basis, j, a, D0) for j in 1:nbasis_x1]
x2plus = [evaluate(x2basis, j, b, D0) for j in 1:nbasis_x2]
x2minus = [evaluate(x2basis, j, a, D0) for j in 1:nbasis_x2]
y1plus = [evaluate(y1basis, j, b, D0) for j in 1:nbasis_y1]
y1minus = [evaluate(y1basis, j, a, D0) for j in 1:nbasis_y1]
y2plus = [evaluate(y2basis, j, b, D0) for j in 1:nbasis_y2]
y2minus = [evaluate(y2basis, j, a, D0) for j in 1:nbasis_y2]

Qx_plus = (+1.0) .* (x2plus * x1plus')
Qx_minus = (-1.0) .* (x2minus * x1minus')
Qy_plus =  (+1.0) .* (y2plus * y1plus')
Qy_minus =(-1.0) .* (y2minus * y1minus')

Dx1 = Kx - alpham * Qx_minus - alphap * Qx_plus
Dy1 = Ky - alpham * Qy_minus - alphap * Qy_plus

f0 = 60000.0  
t0 = 1.0 / f0
function ricker_wavelet(t::Float64, f0::Float64, t0::Float64)
    arg = π * f0 * (t - t0)
    return (1.0 - 2.0 * arg^2) * exp(-arg^2)
end

source_ix = round(Int, Nx / 2)
source_iy = round(Int, Ny / 2) 


Ux = zeros(nbasis_x1, nbasis_y1)
Uy = zeros(nbasis_x1, nbasis_y1)

source_time_function = ricker_wavelet(dt,f0,t0)
F_hat_y = spzeros(nbasis_x1, nbasis_y1)
x_s = x[source_ix]; y_s = y[source_iy]
F_hat_x = project_point_source(x1basis, y1basis,LX1, LY1,x_s, y_s, source_time_function)

Uxhat1=Ux+(dt^2/(2*rho))*(F_hat_x)
Uyhat1=Uy+(dt^2/(2*rho))*(F_hat_y)

Uxhat_old=Ux
Uxhat=Uxhat1
Uyhat_old=Uy
Uyhat=Uyhat1

Ux=from_hat(Uxhat,LX1,LY1)
Uy=from_hat(Uyhat,LX1,LY1)

anim = Animation()
frame_rate = 50

for i=1:n
     
    global Uxhat_old, Uxhat, Uyhat_old, Uyhat, Ux, Uy , F_hat_x, F_hat_y , source_time_function , x_s , y_s

    source_time_function = ricker_wavelet((i+1) * dt,f0,t0)
    F_hat_y = spzeros(nbasis_x1, nbasis_y1)
    x_s = x[source_ix]; y_s = y[source_iy]
    F_hat_x = project_point_source(x1basis, y1basis,LX1, LY1, x_s, y_s, source_time_function)

    Ustarx10= Dx1*Ux*Dy0'
    Ustarx01= Dx0*Ux*Dy1'
    Ustary10= Dx1*Uy*Dy0'
    Ustary01= Dx0*Uy*Dy1'

    Uhat_x10=to_hat(Ustarx10,LX2,LY2)
    Uhat_x01=to_hat(Ustarx01,LX2,LY2)
    Uhat_y10=to_hat(Ustary10,LX2,LY2)
    Uhat_y01=to_hat(Ustary01,LX2,LY2)

    Sxx_hat=(λ+2*μ)*(Uhat_x10)+λ*(Uhat_y01)
    Syy_hat=(λ)*(Uhat_x10)+(λ+2*μ)*(Uhat_y01)
    Sxy_hat=(μ)*(Uhat_x01+Uhat_y10)

    Sxx=from_hat(Sxx_hat,LX2,LY2)
    Sxy=from_hat(Sxy_hat,LX2,LY2)
    Syx=from_hat(Sxy_hat,LX2,LY2)
    Syy=from_hat(Syy_hat,LX2,LY2)

    Sxx10_star=-(Dx1'*Sxx*Dy0)
    Sxy01_star=-(Dx0'*Sxy*Dy1)
    Syx10_star=-(Dx1'*Syx*Dy0)
    Syy01_star=-(Dx0'*Syy*Dy1)

    Sxx10_hat=to_hat(Sxx10_star,LX1,LY1)
    Syx10_hat=to_hat(Syx10_star,LX1,LY1)
    Sxy01_hat=to_hat( Sxy01_star,LX1,LY1)
    Syy01_hat=to_hat( Syy01_star,LX1,LY1)

    Uxhat_new=2*Uxhat-Uxhat_old+(dt^2/rho)*(Sxx10_hat+Sxy01_hat+F_hat_x)
    Uyhat_new=2*Uyhat-Uyhat_old+(dt^2/rho)*(Syx10_hat+Syy01_hat+F_hat_y)

    Uxhat_old=Uxhat
    Uxhat=Uxhat_new
    Uyhat_old=Uyhat
    Uyhat=Uyhat_new

    Ux=from_hat(Uxhat,LX1,LY1)
    Uy=from_hat(Uyhat,LX1,LY1)

    if i % frame_rate == 0
       current_t = i * dt
            
       B_x_vals = [evaluate(x1basis, i, val) for val in x, i in 1:nbasis_x1]
       B_y_vals = [evaluate(y1basis, i, val) for val in y, i in 1:nbasis_y1]
       Ux_grid = B_x_vals * Ux * B_y_vals'
       Uy_grid = B_x_vals * Uy * B_y_vals'

       max_abs_val=5*10^(-13)
       p = heatmap(x, y, Ux_grid', aspect_ratio=:equal, c=custom_colors,clim=(-max_abs_val, max_abs_val), title="Horizontal Displacement (Ux) | t=$(round(current_t*1000000, sigdigits=2)) μs",
       xlabel="x (m)",ylabel="y (m)",size=(900, 700))
       scatter!(p, [x[source_ix]],[y[source_iy]], marker=:diamond, markersize=5,  markercolor=:red)
       frame(anim, p)
       println("Frame $(i/frame_rate)/$(n/frame_rate) at t=$(round(current_t*1000000, sigdigits=2)) μs captured.")
    end
end
gif(anim, "Elastic_Wave_DFDM_2D_1Domain.gif", fps=10)

