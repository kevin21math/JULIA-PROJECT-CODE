# basis: Legendre-Gauss-Lobatto
using OrdinaryDiffEqLowStorageRK
using Trixi, LinearAlgebra, Plots
using DelimitedFiles
polydeg = 3 #= polynomial degree =#
basis = LobattoLegendreBasis(polydeg)
nodes = basis.nodes # Gauss-Lobatto nodes in [-1, 1]
D = basis.derivative_matrix
M = diagm(basis.weights) # mass matrix
B = diagm([-1; zeros(polydeg - 1); 1])
plt = plot()
#errors
dx_vals = Float64[]
L1_errors = Float64[]
L2_errors = Float64[]
Linf_errors = Float64[]

# mesh
coordinates_min = -1.0 # minimum coordinate
coordinates_max = 1.0  # maximum coordinate
set = [8,16,32,64,128,256] # number of elements to test
for n_elements in set   # number of elements
    dx = (coordinates_max - coordinates_min) / n_elements # length of one element

    x = Matrix{Float64}(undef, length(nodes), n_elements)
    for element in 1:n_elements
        x_l = -1 + (element - 1) * dx + dx / 2
        for i in eachindex(nodes) # basis points in [-1, 1]
            ξ = nodes[i]
            x[i, element] = x_l + dx / 2 * ξ
        end
    end    
   # initial condition
   initial_condition_sine_wave(x) = 1.0 + 0.5 * sin(pi * x)
   u0 = initial_condition_sine_wave.(x)
   if(n_elements == 16)
     plot!(vec(x), vec(u0), label = "exact solution", legend = :topleft)
   end
   # flux Lax-Friedrichs
   surface_flux = flux_lax_friedrichs

   # rhs! method
   function rhs!(du, u, x, t)
      # reset du
      du .= zero(eltype(du))
      flux_numerical = copy(du)

      # calculate interface and boundary fluxes
      equations = LinearScalarAdvectionEquation1D(1.0)
      for element in 2:(n_elements - 1)
          # left interface
          flux_numerical[1, element] = surface_flux(u[end, element - 1], u[1, element], 1, equations)

          flux_numerical[end, element - 1] = flux_numerical[1, element]
          # right interface
          flux_numerical[end, element] = surface_flux(u[end, element], u[1, element + 1], 1, equations)
        
          flux_numerical[1, element + 1] = flux_numerical[end, element]
      end
      # boundary flux
      flux_numerical[1, 1] = surface_flux(u[end, end], u[1, 1], 1, equations)
      flux_numerical[end, end] = flux_numerical[1, 1]

      # calculate surface integrals
      for element in 1:n_elements
          du[:, element] -= (M \ B) * flux_numerical[:, element]
      end

      # calculate volume integral
      for element in 1:n_elements
          flux = u[:, element]
          du[:, element] += (M \ transpose(D)) * M * flux
      end

      # apply Jacobian from mapping to reference element
      for element in 1:n_elements
          du[:, element] *= 2 / dx
      end

      return nothing
   end

   # create ODE problem
   tspan = (0.0, 2.0)
   ode = ODEProblem(rhs!, u0, tspan, x)

   # solve
   sol = solve(ode, RDPK3SpFSAL49(); abstol = 1.0e-15, reltol = 1.0e-15, ode_default_options()...)

   plot!(vec(x), vec(sol.u[end]), label = "solution at n_elements=$(n_elements)", legend = :topleft, lw = 2)

   diff = abs.(sol.u[end] - initial_condition_sine_wave.(x))
   L1 = (dx / 3) * (diff[1] + 4 * sum(diff[2:2:end-1]) + 2 * sum(diff[3:2:end-2]) + diff[end])
   diff_sq = (sol.u[end] - initial_condition_sine_wave.(x)).^2
   L2_sq = (dx / 3) * (diff_sq[1] + 4 * sum(diff_sq[2:2:end-1]) + 2 * sum(diff_sq[3:2:end-2]) + diff_sq[end])
   L2 = sqrt(L2_sq)
   Linf = maximum(abs.(sol.u[end] - initial_condition_sine_wave.(x)))
   push!(dx_vals, dx)
   push!(L1_errors, L1)
   push!(L2_errors, L2)
   push!(Linf_errors, Linf)
end

function convergence_rate(errors, dx_vals)
    rates = Any[]
    for i in 1:6
        if i == 1
            push!(rates, "" )  
        else
            rate = log2(errors[i-1]/errors[i])
            push!(rates, rate)
        end
    end
    return rates
end

L1_rates = convergence_rate(L1_errors, dx_vals)
L2_rates = convergence_rate(L2_errors, dx_vals)
Linf_rates = convergence_rate(Linf_errors, dx_vals)

println("L1 errors: ", L1_errors)
println("L2 errors: ", L2_errors)
println("L∞ errors: ", Linf_errors)
println("L1 convergence rates: ", L1_rates)
println("L2 convergence rates: ", L2_rates)
println("L∞ convergence rates: ", Linf_rates)
display(plt)
savefig(plt, "dgsem.png")
header = ["dx", "L1_error", "L1_rate","L2_error", "L2_rate","Linf_error", "Linf_rate"]
data = [dx_vals L1_errors L1_rates L2_errors L2_rates Linf_errors Linf_rates]
open("dgsem_data.csv", "w") do io
    writedlm(io, [header], ',')
end
open("dgsem_data.csv", "a") do io
    writedlm(io, data, ',')
end
