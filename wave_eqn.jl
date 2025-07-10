using Plots


nx, nz = 201, 201
dx, dz = 10.0, 10.0
nt = 1000
dt = 0.001  
T = collect(0:dt:(nt-1)*dt)

ρ = 2500.0
vp = 4000.0
nu = 0.2
vs = vp * sqrt((0.5 - nu) / (1 - nu))
μ = ρ * vs^2
λ = ρ * (vp^2 - 2 * vs^2)


v_x = zeros(Float64, nx, nz)
v_z = zeros(Float64, nx-1, nz-1)
σ_xx = zeros(Float64, nx-1, nz)
σ_zz = zeros(Float64, nx-1, nz)
σ_xz = zeros(Float64, nx, nz-1)


isrc, jsrc = div(nx, 2), div(nz, 2)
ricker(t) = (1 - 2 * (π*15*t)^2) * exp(-(π*15*t)^2)
g(t) = -2*40*t*exp(-40*t*t)

ir, jr = isrc + 10, jsrc
u_r = zeros(Float64, nt)
for it = 1:nt
    t = it * dt
    
    σ_xx[isrc, jsrc] += ricker(t - 0.1)
    σ_zz[isrc, jsrc] += ricker(t - 0.1)

    
    for i in 2:nx-1
        for j in 2:nz-1
            v_x[i, j] += dt / ρ * ((σ_xx[i, j] - σ_xx[i-1, j]) / dx + (σ_xz[i, j] - σ_xz[i, j-1]) / dz)
        end
    end

    for i in 1:nx-1
        for j in 1:nz-1
            v_z[i, j] += dt / ρ * ((σ_xz[i+1, j] - σ_xz[i, j]) / dx + (σ_zz[i, j+1] - σ_zz[i, j]) / dz)
        end
    end

    
    for j in 1:nz
        v_x[end, j] = v_x[end, j] + vp*(dt / dx) * (v_x[end, j] - v_x[end-1, j])
        v_x[1, j] = v_x[1, j] - vp*(dt / dx) * (v_x[2, j] - v_x[1, j])
    end
   
    for i in 1:nx
        v_x[i, end] = v_x[i, end] + vp*(dt / dz) * (v_x[i, end] - v_x[i, end-1])
        v_x[i, 1] = v_x[i, 1] - vp*(dt / dz) * (v_x[i, 2] - v_x[i, 1])
    end

    
    for i in 1:nx-1
        for j in 2:nz-1
            dvx_dx = (v_x[i+1, j] - v_x[i, j]) / dx
            dvz_dz = (v_z[i, j] - v_z[i, j-1]) / dz
            σ_xx[i, j] += dt * ((λ + 2μ) * dvx_dx + λ * dvz_dz)
            σ_zz[i, j] += dt * ((λ + 2μ) * dvz_dz + λ * dvx_dx)
        end
    end
    
    for i in 1:nx-1
         σ_xx[i, end] = σ_xx[i, end] + vp*(dt / dz) * (σ_xx[i, end] - σ_xx[i, end-1])
         σ_xx[i, 1] = σ_xx[i, 1] - vp*(dt / dz) * (σ_xx[i, 2] - σ_xx[i, 1])
         σ_zz[i, end] = σ_zz[i, end] + vp*(dt / dz) * (σ_zz[i, end] - σ_zz[i, end-1])
         σ_zz[i, 1] = σ_zz[i, 1] - vp*(dt / dz) * (σ_zz[i, 2] - σ_zz[i, 1])
    end

    for i in 2:nx-1
        for j in 1:nz-1
            dvx_dz = (v_x[i, j+1] - v_x[i, j]) / dz
            dvz_dx = (v_z[i, j] - v_z[i-1, j]) / dx
            σ_xz[i, j] += dt * μ * (dvx_dz + dvz_dx)
        end
    end
    
    for j in 1:nz-1
        σ_xz[end, j] = σ_xz[end, j] + vp*(dt / dx) * (σ_xz[end, j] - σ_xz[end-1, j])
        σ_xz[1, j] = σ_xz[1, j] - vp*(dt / dx) * (σ_xz[2, j] - σ_xz[1, j])
    end
    
    if it % 50 == 0
        plot = heatmap(σ_xx', clims=(-1e-3, 1e-3), title="Time step: $it", aspect_ratio=1)
        display(plot)
    end

    u_r[it] = v_x[ir, jr]
end
u_disp = cumsum(u_r) .* dt
plot(T, u_disp ./ maximum(abs.(u_disp)), lw=2, xlabel="Time (s)", ylabel="Normalized displacement", title="2D Elastic Wave - Displacement Seismogram", legend=false)
