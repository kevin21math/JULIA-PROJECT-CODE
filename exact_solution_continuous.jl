using DelimitedFiles

Nx = 200
L  = 2
a  = 1
T  = 1
CFL = 0.5

dx = L / Nx
dt = CFL * dx / a
Nt = Int(floor(T / dt))

x = LinRange(0, L, Nx)
U_exact = zeros(Nt+1, Nx)

for n in 0:Nt
    t = n * dt
    U_exact[n+1, :] = sin.(2π * mod.(x .- a*t, 2)) 
end

writedlm("exact_solution.csv", U_exact, ',')
