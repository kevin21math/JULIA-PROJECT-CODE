using Plots

nx, nz = 401, 401
dx, dz = 15.0, 15.0
nt = 1000
dt = 0.001  

ρ = 2500.0
vp = 4000.0
nu = 0.25
vs = vp * sqrt((0.5 - nu) / (1 - nu))
μ = ρ * vs^2
λ = ρ * (vp^2 - 2 * vs^2)


v_x = zeros(Float64, nx, nz)
v_z = zeros(Float64, nx-1, nz-1)
τ_xx = zeros(Float64, nx-1, nz)
τ_zz = zeros(Float64, nx-1, nz)
τ_xz = zeros(Float64, nx, nz-1)


isrc, jsrc = div(nx, 2), div(nz, 2)
ricker(t) = (1 - 2 * (π*15*t)^2) * exp(-(π*15*t)^2)

for it = 1:nt
    t = it * dt

    τ_xx[isrc, jsrc] += ricker(t - 0.1)
    τ_zz[isrc, jsrc] += ricker(t - 0.1)

    
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

    
    for j in 1:nz
        v_x[end, j] = v_x[end, j] - vp * (dt/dx) * (v_x[end, j] - v_x[end-1, j])
        v_x[1, j] = v_x[1, j] + vp * (dt/dx) * (v_x[2, j] - v_x[1, j])
    end
   
    for i in 1:nx
        v_x[i, end] = v_x[i, end] - vp * (dt/dz) * (v_x[i, end] - v_x[i, end-1])
        v_x[i, 1] = v_x[i, 1] + vp * (dt/dz) * (v_x[i, 2] - v_x[i, 1])
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
         τ_xx[i, 1] = τ_xx[i, 1] + vp * (dt / dz) * (τ_xx[i, 2] - τ_xx[i, 1])
         τ_zz[i, end] = τ_zz[i, end] - vp * (dt / dz) * (τ_zz[i, end] - τ_zz[i, end-1])
         τ_zz[i, 1] = τ_zz[i, 1] + vp * (dt / dz) * (τ_zz[i, 2] - τ_zz[i, 1])
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
    
    if it % 250 == 0
        plot = heatmap(τ_xx', clims=(-1e-3, 1e-3),xlabel="x", ylabel="z", title="wavefront of τ_xx at time t =$(round(it*dt, digits=3)) s ", aspect_ratio=1)
        display(plot)
        savefig("Derivative frame_t$(round(it*dt, digits=3)).png") 
    end  
end
