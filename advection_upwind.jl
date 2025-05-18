# Advection Upwind Scheme
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
    for i in 2:Nx
        u_new[i] = u[i] - CFL*(u[i] - u[i-1]) 
    end
    u_new[1] = u[1] - CFL*(u[1] - u[Nx])
    #u_new[1] = u_new[Nx]
    u .= u_new
    plot(x,u,label="sine curve",xlabel="x",ylabel="u",xlims = (0, L),ylims = (-1, 1), title="1D Linear Advection")
end
gif(anim,"advection_upwind.gif",fps=20)