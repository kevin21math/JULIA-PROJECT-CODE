using LinearAlgebra,Plots

nx, nz = 601, 601
dx, dz = 25.0, 25.0
nt = 500
dt = 0.004  
T = collect(0:dt:(nt-1)*dt)
X = collect(-(nx-2)*dx/200:dx/100:(nx-2)*dx/200)
Z = collect(0:dz/100:(nz-1)*dz/100)

ρ = 2500.0
vp = 4000.0
nu = 0.25
vs = vp * sqrt((0.5 - nu) / (1 - nu))
μ = ρ * vs^2
λ = ρ * (vp^2 - 2 * vs^2)


v_x = zeros(nx, nz)
v_z = zeros(nx-1, nz-1)
τ_xx = zeros(nx-1, nz)
τ_zz = zeros(nx-1, nz)
τ_xz = zeros(nx, nz-1)


isrc, jsrc = div(nx, 2), 1
f(t) = exp.(-200*(t^2))

ir1, jr1 = isrc + Int(1500/dx),3
ir2, jr2 = isrc + Int(2000/dx),3
ir3, jr3 = isrc + Int(2500/dx),3
ir4, jr4 = isrc + Int(3000/dx),3
u_r1 = zeros(nt)
u_r2 = zeros(nt)
u_r3 = zeros(nt)
u_r4 = zeros(nt)
for it = 1:nt
    t = it * dt

    τ_xx[isrc, jsrc] += f(t - 0.1)
    τ_zz[isrc, jsrc] += f(t - 0.1)

    
    for i in 2:nx-1
        for j in 2:nz-1
            v_x[i, j] += dt / ρ * ((τ_xx[i, j] - τ_xx[i-1, j]) / dx + (τ_xz[i, j] - τ_xz[i, j-1]) / dz)
        end
    end

    for i in 1:nx-1
        for j in 1:nz-1
            v_z[i, j] += dt / ρ * ((τ_xz[i+1, j] - τ_xz[i, j]) / dx + (τ_zz[i, j+1] - τ_zz[i, j]) / dz)
        end
    end

    
    for j in 2:nz
        v_x[end, j] = v_x[end, j] - vp * (dt/dx) * (v_x[end, j] - v_x[end-1, j])
        v_x[1, j] = v_x[1, j] + vp * (dt/dx) * (v_x[2, j] - v_x[1, j])
    end
   
    for i in 1:nx
        v_x[i, end] = v_x[i, end] - vp * (dt/dz) * (v_x[i, end] - v_x[i, end-1])
    end

    
    for i in 1:nx-1
        for j in 2:nz-1
            dvx_dx = (v_x[i+1, j] - v_x[i, j]) / dx
            dvz_dz = (v_z[i, j] - v_z[i, j-1]) / dz
            τ_xx[i, j] += dt * ((λ + 2μ) * dvx_dx + λ * dvz_dz)
            τ_zz[i, j] += dt * ((λ + 2μ) * dvz_dz + λ * dvx_dx)
        end
    end
    
    for i in 1:nx-1
         τ_xx[i, end] = τ_xx[i, end] - vp * (dt / dz) * (τ_xx[i, end] - τ_xx[i, end-1])
         τ_zz[i, end] = τ_zz[i, end] - vp * (dt / dz) * (τ_zz[i, end] - τ_zz[i, end-1])
    end

    for i in 2:nx-1
        for j in 1:nz-1
            dvx_dz = (v_x[i, j+1] - v_x[i, j]) / dz
            dvz_dx = (v_z[i, j] - v_z[i-1, j]) / dx
            τ_xz[i, j] += dt * μ * (dvx_dz + dvz_dx)
        end
    end
    
    for j in 1:nz-1
        τ_xz[end, j] = τ_xz[end, j] - vp * (dt / dx) * (τ_xz[end, j] - τ_xz[end-1, j])
        τ_xz[1, j] = τ_xz[1, j] + vp * (dt / dx) * (τ_xz[2, j] - τ_xz[1, j])
    end

    τ_zz[:,1] .= 0.0
    τ_xz[:,1] .= -τ_xz[:,2]

    if it % 125 == 0
        plot = heatmap(X,Z,τ_xx',xlims =(-75, 75), clims=(-1e-3, 1e-3),color=:seismic, xlabel="x", ylabel="z", aspect_ratio=1) 
        display(plot)
        savefig("free surface frame_t$(round(it*dt, digits=3)).png") 
    end

    u_r1[it] = v_x[ir1, jr1]
    u_r2[it] = v_x[ir2, jr2]
    u_r3[it] = v_x[ir3, jr3]
    u_r4[it] = v_x[ir4, jr4]
end
u_disp1 = cumsum(u_r1) .* dt
u_disp2 = cumsum(u_r2) .* dt
u_disp3 = cumsum(u_r3) .* dt
u_disp4 = cumsum(u_r4) .* dt

u_disp1 = u_disp1 / maximum(abs.(u_disp1))
u_disp2 = u_disp2 / maximum(abs.(u_disp2))
u_disp3 = u_disp3 / maximum(abs.(u_disp3))
u_disp4 = u_disp4 / maximum(abs.(u_disp4))

sol1 = plot(T, u_disp1, lw=2, ylim=(-1, 1), xlabel="Time (s)", ylabel="Normalized displacement", title="Horizontal displacement component at 1500m", legend=false)
display(sol1)
savefig("Displacement_1500m.png")

sol2 = plot(T, u_disp2, lw=2, ylim=(-1, 1), xlabel="Time (s)", ylabel="Normalized displacement", title="Horizontal displacement component at 2000m", legend=false)
display(sol2)
savefig("Displacement_2000m.png")

sol3 = plot(T, u_disp3, lw=2, ylim=(-1, 1), xlabel="Time (s)", ylabel="Normalized displacement", title="Horizontal displacement component at 2500m", legend=false)
display(sol3)
savefig("Displacement_2500m.png")

sol4 = plot(T, u_disp4, lw=2, ylim=(-1, 1), xlabel="Time (s)", ylabel="Normalized displacement", title="Horizontal displacement component at 3000m", legend=false)
display(sol4)
savefig("Displacement_3000m.png")

