using SparseArrays, BSplineKit, LinearAlgebra, QuadGK, Plots

fixzeros(x) = x == 0.0 ? 0.0 : x
fixsmall(x) = abs(x) < 1e-3 ? 0.0 : x

include("C:\\Users\\hp\\Downloads\\VS CODE PRG\\julia project\\DFD_Method\\shifted_knots.jl")
include("C:\\Users\\hp\\Downloads\\VS CODE PRG\\julia project\\DFD_Method\\galerkin_cross.jl")
include("C:\\Users\\hp\\Downloads\\VS CODE PRG\\julia project\\DFD_Method\\galerkin_cross_elementwise.jl")

const Gal = BSplineKit.Galerkin
v = 3000
dt = 2.55*10^(-3)
#dt = 1.0*10^(-4)
rho = 2500
kval = rho*(v^2)
T = 1.0
nt = Int(round(T/dt))

a, b = 0, 6000
N = 157
#N = 751
x  = LinRange(0,6000,N)
p = 4 
αk1 = 0.125
αk2 = -0.125
Nd = 3
dom_nodes = [1, 53, 105, 157]
#dom_nodes = [1, 251, 501, 751]
D1 = BSplineKit.Derivative(1)
D0 = BSplineKit.Derivative(0)
alpha  = [0, 0.5, 0.5, 1]

L1     = Vector{Matrix{Float64}}(undef, Nd)
L2     = Vector{Matrix{Float64}}(undef, Nd)
Dstar  = Vector{Matrix{Float64}}(undef, Nd)
D21    = Vector{Matrix{Float64}}(undef, Nd)
q1_minus = Vector{Vector{Float64}}(undef, Nd)
q1_plus  = Vector{Vector{Float64}}(undef, Nd)
q2_minus = Vector{Vector{Float64}}(undef, Nd)
q2_plus  = Vector{Vector{Float64}}(undef, Nd)
uhat   = Vector{Vector{Float64}}(undef, Nd)
uhat_old = Vector{Vector{Float64}}(undef, Nd)
uhat_new = Vector{Vector{Float64}}(undef, Nd)
σhat   = Vector{Vector{Float64}}(undef, Nd)
B_minus = Vector{Matrix{Float64}}(undef, Nd)
B_plus  = Vector{Matrix{Float64}}(undef, Nd)
Phi_list = Vector{Matrix{Float64}}(undef, Nd)
nbasis1_list = Vector{Int}(undef, Nd) 
nloc_list = Vector{Int}(undef, Nd)
u_global = zeros(N)

for i in 1:Nd
    xloc=dom_nodes[i]:dom_nodes[i+1]
    nloc=length(xloc)
    nloc_list[i] = nloc 
    knots1loc = shifted_knots(x[xloc[1]],x[xloc[end]],nloc,p,αk1)
    basis1loc = BSplineBasis(p, knots1loc[p:end-p+1])
    nbasis1 = length(basis1loc)
    nbasis1_list[i] = nbasis1
    M1loc = galerkin_matrix(basis1loc, (D0,D0))
    L1[i] = cholesky(M1loc).L
    #basis1loc = L1[i]\basis1loc
    #Phi = [evaluate(basis1loc, j, x) for j in 1:nbasis1, x in x[xloc]]
    #Psi = L1[i] \ Phi
    #println(Psi*Psi)
    knots2loc = shifted_knots(x[xloc[1]],x[xloc[end]],nloc,p,αk2)
    basis2loc = BSplineBasis(p, knots2loc[p:end-p+1])
    nbasis2 = length(basis2loc)
    M2loc = galerkin_matrix(basis2loc, (D0,D0))
    L2[i] = cholesky(M2loc).L
    #basis2loc = L2[i]\basis2loc
    
    K1loc = spzeros(nbasis2, nbasis1)
    #K1loc = (-1).*galerkin_cross_elementwise(basis2loc, basis1loc, (D1, D0); nquad=8)
    #K1loc = (-1).*galerkin_cross(basis2loc, basis1loc, (D1, D0); a = x[xloc[1]], b = x[xloc[end]], nquad = 24)
    K1loc = galerkin_cross(basis2loc, basis1loc, (D1, D0); a = x[xloc[1]], b = x[xloc[end]], nquad = 24)
    
    K1loc = fixzeros.(K1loc)
    
    x_left  = x[xloc[1]]
    x_right = x[xloc[end]]
    
    q1_minus[i]  = [evaluate(basis1loc, i,x_left)  for i in 1:nbasis1]            
    q1_plus[i] = [evaluate(basis1loc,i,x_right) for i in 1:nbasis1]
    
    q2_minus[i]  = [evaluate(basis2loc,i, x_left) for i in 1:nbasis2]                
    q2_plus[i] = [evaluate(basis2loc,i, x_right) for i in 1:nbasis2]
    
    Qminus =  (-1) .* (q2_minus[i]  * q1_minus[i]') 
    Qminus = fixzeros.(Qminus) 
    Qplus  =  (+1) .* (q2_plus[i] * q1_plus[i]')
      
    #Dstar[i]  = K1loc+(1 - alpha[i]).*Qminus+(1-alpha[i+1]).*Qplus
    Dstar[i]  = K1loc - alpha[i].*Qminus - (1 - alpha[i+1]).*Qplus
    
    D21[i] = LowerTriangular(L2[i]) \ ((LowerTriangular(L1[i]) \ transpose(Dstar[i]))')
    #D21[i] = fixsmall.(D21[i])
    
    f0 = 10.0
    ξ = @. π*f0*(x[xloc] - 500) / v
    u0_vals = @. (1 - 2*ξ^2) * exp(-ξ^2)
    u0_vals = ifelse.(abs.(u0_vals) .< 1e-10, 0.0, u0_vals)
    
    x_loc_pts = x[xloc]
    
    Phi = zeros(nbasis1, nloc)
    for j in 1:nbasis1
       Phi[j, :] .= evaluate(basis1loc, j, x_loc_pts)
    end
    Phi_list[i] = Phi
       
    dx_local = diff(x_loc_pts)
    w = zeros(nloc)
    w[1] = 0.5
    w[end] = 0.5
    w[2:end-1] .= 1.0
    w.= w .* vcat(dx_local, dx_local[end])  
    r = Phi * (w .* u0_vals)
    
    fB =  UpperTriangular(transpose(L1[i])) \ (LowerTriangular(L1[i]) \ r)
              
    fhat_coeffs = transpose(L1[i]) * fB   
    
    uhat[i] = copy(fhat_coeffs)
    uhat_old[i] = copy(uhat[i])
    uhat_new[i] = copy(uhat[i])
    σhat[i] = zeros(nloc)
    
    #plt = plot(x[xloc], uhat[i], lw=2, xlabel="x", ylabel="u", title="Initial condition in element $i")
    #display(plt)
      
end


for i in 2:Nd
    B_minus[i] = alpha[i] .* (LowerTriangular(L2[i]) \ (q2_minus[i] * transpose(LowerTriangular(L1[i-1]) \ q1_plus[i-1])))
    B_minus[i] = fixsmall.(B_minus[i])
end
for i in 1:Nd-1
    B_plus[i]  = (1 - alpha[i+1]) .*(LowerTriangular(L2[i]) \ (q2_plus[i] * transpose(LowerTriangular(L1[i+1]) \ q1_minus[i+1])))
    B_plus[i] = fixsmall.(B_plus[i])
end

B_minus[1] = zeros(size(B_minus[2]))
B_plus[Nd] = zeros(size(B_plus[Nd-1]))



snapshots = []


for i in 1:nt
    σhat[1] .= kval.*(D21[1]*uhat[1] + B_plus[1]*uhat[2])
    σhat[2] .= kval.*(D21[2]*uhat[2] - B_minus[2]*uhat[1] + B_plus[2]*uhat[3])
    σhat[3] .= kval.*(D21[3]*uhat[3] - B_minus[3]*uhat[2])
    
    if(i == 1)
        uhat_new[1] .= uhat[1]  + ((dt^2)*0.5/rho) .* (-transpose(D21[1]) * σhat[1] + transpose(B_minus[2])*σhat[2])
        uhat_new[2] .= uhat[2]  + ((dt^2)*0.5/rho) .* (-transpose(D21[2]) * σhat[2] - transpose(B_plus[1])*σhat[1] + transpose(B_minus[3])*σhat[3])
        uhat_new[3] .= uhat[3]  + ((dt^2)*0.5/rho) .* (-transpose(D21[3]) * σhat[3] - transpose(B_plus[2])*σhat[2])
    else
        uhat_new[1] .= 2*uhat[1] - uhat_old[1] + ((dt^2)/rho) .* (-transpose(D21[1]) * σhat[1] + transpose(B_minus[2])*σhat[2])
        uhat_new[2] .= 2*uhat[2] - uhat_old[2] + ((dt^2)/rho) .* (-transpose(D21[2]) * σhat[2] - transpose(B_plus[1])*σhat[1] + transpose(B_minus[3])*σhat[3])
        uhat_new[3] .= 2*uhat[3] - uhat_old[3] + ((dt^2)/rho) .* (-transpose(D21[3]) * σhat[3] - transpose(B_plus[2])*σhat[2])
    end
     
    #if i % 5 == 0
        for i in 1:Nd
            fB = transpose(L1[i]) \ uhat[i]
            Phi = Phi_list[i]                  
            u_phys = transpose(Phi) * fB 
            rng = dom_nodes[i]:dom_nodes[i+1] 
            u_global[rng] = u_phys
        end
        push!(snapshots, copy(u_global))
    #end

    for j in 1:Nd
        uhat_old[j] = copy(uhat[j])
        uhat[j] = copy(uhat_new[j])
    end
end

#figure=plot(x, snapshots[end], lw=2, xlabel="x", ylabel="u", title="Wavefield at final time")
#display(figure)


anim = @animate for k in 1:nt
    plot(x, snapshots[k], ylim = (-0.5, 1.0), xlim = (0, 6000), xlabel = "x", ylabel = "u", title = "time step $(k)", legend = false, lw = 2)
end

gif(anim, "wavefield.gif", fps = 10)

