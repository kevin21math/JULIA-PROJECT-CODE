# Advection Central Scheme
# ut + a*ux = 0
# u(x,0) = sin(2πx)

using Plots

Nx = 100 
L  = 1   
a  = 1   
T  = 1  
CFL = 0.5

dx = L/Nx                    
dt = CFL*dx/a               
Nt = T/dt                

x = LinRange(0,L,Nx)

u0(x) = sin(2π*x)
u = u0.(x)
u_new = copy(u)

anim = @animate for n in 1:Nt
    for i in 2:Nx-1
        u_new[i] = u[i] - CFL*(u[i+1] - u[i-1]) / 2
    end
    u_new[Nx]= u[Nx] - CFL*(u[1] - u[Nx-1]) / 2
    u_new[1] = u[1] - CFL*(u[2] - u[Nx]) / 2
    u .= u_new
    plot(x,u,label="sine curve",xlabel="x",ylabel="u",title="1D Linear Advection")
end
gif(anim,"advection_central.gif",fps=20)