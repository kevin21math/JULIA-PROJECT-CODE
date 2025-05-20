# Continuous periodic function
# Advection Upwind_Central Scheme
# ut + a*ux = 0
# u(x,0) = sin(2πx)

using DelimitedFiles

Nx = 200 
L  = 2 
a  = 1   
T  = 1  
CFL = 0.5

dx = L/Nx                    
dt = CFL*dx/a              
Nt = Int(floor(T/dt))       

x = LinRange(0, L, Nx) 
u0(x) = sin(2π*x)
u = u0.(x)
u_new = copy(u)

U = zeros(Nt+1, Nx)
U[1, :] = u 

for n in 1:Nt
    for i in 2:Nx-1
        u_new[i] = u[i] - CFL*(u[i+1] - u[i-1])/2 + (dt/dx)*abs(a)*(u[i+1] - 2*u[i] + u[i-1])/2
    end
   
    u_new[Nx] = u[Nx] - CFL*(u[1] - u[Nx-1])/2 + (dt/dx)*abs(a)*(u[1] - 2*u[Nx] + u[Nx-1])/2
    u_new[1]  = u[1]  - CFL*(u[2] - u[Nx])/2   + (dt/dx)*abs(a)*(u[2] - 2*u[1] + u[Nx])/2
   
    u .= u_new
    U[n+1, :] = u

    t = n*dt
    exact_u = u0.(mod.(x .- a*t, L))

    if n == Nt
     diff = abs.(u .- exact_u)
     L1_norm = (dx / 3) * (diff[1] + 4 * sum(diff[2:2:Nx-1]) + 2 * sum(diff[3:2:Nx-2]) + diff[Nx])
     println("L1 norm (Simpson) at time $t: $L1_norm")

     diff_sq = (u .- exact_u).^2
     L2_sq = (dx / 3) * (diff_sq[1] + 4 * sum(diff_sq[2:2:Nx-1]) + 2 * sum(diff_sq[3:2:Nx-2]) + diff_sq[Nx])
     L2_norm = sqrt(L2_sq)
     println("L2 norm (Simpson) at time $t: $L2_norm")

     L_inf_norm = maximum(abs.(u - exact_u))
     println("L∞ norm at time $t: $L_inf_norm")
    end
end

writedlm("advection_upwind_central_data.csv", U, ',')
