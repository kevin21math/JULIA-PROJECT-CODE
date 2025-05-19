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


u0(x) = ifelse.(x .< 0.5, 1.0, 0.0)
U_exact[1, :] = u0.(x)

for n in 1:Nt
    t = n * dt
    #shifted_x = mod.(x .- a*t, L)
    U_exact[n+1, :] = u0.(x .- a*t)
end

writedlm("dis_exact_solution.csv", U_exact, ',')
