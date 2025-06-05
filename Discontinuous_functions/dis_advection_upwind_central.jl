# Discontinuos function
# Advection Upwind_Central Scheme
# ut + a*ux = 0
# u(x,0) = 1 if x < 0.5, 0 otherwise

using DelimitedFiles

Nx = 200
L  = 2
a  = 1
T  = 1
CFL = 0.5

dx = L/Nx
dt = CFL * dx / a
Nt = Int(floor(T/dt))

x = LinRange(0, L, Nx)

u0(x) = ifelse.(x .< 0.5, 1.0, 0.0)
u = u0.(x)
u_new = copy(u)

U = zeros(Nt+1, Nx) 
U[1, :] = u

for n in 1:Nt
    for i in 2:Nx-1
        u_new[i] = u[i] - CFL*(u[i+1] - u[i-1])/2 + (dt/dx)*abs(a)*(u[i+1] - 2*u[i] + u[i-1])/2
    end
    u_new[Nx] = 0   
    u_new[1]  = 1   
    u .= u_new
    U[n+1, :] = u
end

writedlm("dis_advection_upwind_central_data.csv", U, ',')
