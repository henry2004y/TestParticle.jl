import DisplayAs #hide
using TestParticle, OrdinaryDiffEq
using TestParticle: mᵢ
using LinearAlgebra
using Statistics: mean
using Printf
using CairoMakie
CairoMakie.activate!(type = "png") #hide

const γ = 5/3

"Set initial conditions."
function prob_func(prob, i, repeat)
   v₀ = sample(vdf₁, 1)
   r₀ = [5_000e3, 0.0, 0.0]

   prob = remake(prob; u0 = [r₀..., v₀...])
end

# MHD states in SI units
n₁ = 1.0e6
T₁ = 720471.8506664868
Pi₁ = 0.0049735919716217296 * 1e-9
Pe₁ = 0.0049735919716217296 * 1e-9
Pth₁ = Pi₁ + Pe₁
V₁ = [-545.1484121835928, 0.0, 0.0] .* 1e3
B₁ = [0.0, 0.0, 5.0] .* 1e-9

n₂ = 3.1662479035540008e6
T₂ = 5.957947703036288e6
Pi₂ = 0.13022511153415584 * 1e-9
Pe₂ = 0.13022511153415584 * 1e-9
Pth₂ = Pi₂ + Pe₂
V₂ = [-172.17489874108816, 0.0, 0.0] .* 1e3
B₂ = [0.0, 0.0, 15.831239517770003] .* 1e-9;

E₁ = B₁ × V₁
E₂ = B₂ × V₂
# Thermal speed used to generate the Maxwellian distribution
Uth₁ = √(γ * Pth₁ / (n₁ * mᵢ))
Uth₂ = √(γ * Pth₂ / (n₂ * mᵢ))
# Shock normal direction range
x = range(-5_000e3, 5_000e3, length=100)
B = repeat(B₁, 1, length(x))
E = repeat(E₁, 1, length(x))
# Index for the shock location
mid_ = length(x) ÷ 2

B[:, 1:mid_] .= B₂
E[:, 1:mid_] .= E₂

const vdf₁ = Maxwellian(V₁, Uth₁)

trajectories = 100
weight₁ = n₁ / trajectories # relation between test particle and real particles

prob = let
   # BC type 3 is Flat
   param = prepare(x, E, B; species=Proton, bc=3);
   stateinit = zeros(6) # particle position and velocity to be modified
   tspan = (0.0, 30.0)
   ODEProblem(trace!, stateinit, tspan, param)
end
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)

sols = solve(ensemble_prob, Vern9(), EnsembleSerial(); trajectories);

function plot_traj(sols)
   f = Figure(fontsize=18)
   ax = Axis3(f[1, 1],
      title = "Shock particles",
      xlabel = "x [km]",
      ylabel = "y [km]",
      zlabel = "z [km]",
      aspect = :data,
   )

   invL = 1 / 1e3

   for i in eachindex(sols)
      lines!(ax, sols[i], idxs=(1,2,3), label="$i")
      ##TODO: wait for https://github.com/MakieOrg/Makie.jl/issues/3623 to be fixed!
      ax.scene.plots[9+2*i-1].color = Makie.wong_colors()[mod(i-1, 7)+1]
      scale!(ax.scene.plots[9+2*i-1], invL, invL, invL)
   end

   f
end

f = plot_traj(sols[1:4])
f = DisplayAs.PNG(f) #hide

function plot_dist(x, sols; nxchunks::Int=2, ntchunks::Int=20)
   trange = range(sols[1].prob.tspan..., length=ntchunks)

   xrange = range(x[1], x[end], length=nxchunks+1)
   dx = (x[end] - x[1]) / nxchunks
   xmid = range(x[1] + 0.5dx, x[end] - 0.5dx, length=nxchunks) ./ 1e3

   vx = [Float64[] for _ in 1:nxchunks]

   for sol in sols
      for t in trange
         xv = sol(t)
         for i in 1:nxchunks
            if xrange[i] < xv[1] ≤ xrange[i+1]
               push!(vx[i], xv[4] / 1e3)
            end
         end
      end
   end

   for i in eachindex(vx)
      if isempty(vx[i])
         push!(vx[i], 0.0)
      end
   end

   f = Figure(size = (1200, 600), fontsize=18)
   ax = Axis(f[1, 1],
      limits = (nothing, nothing, -750, 400),
      title = "Phase space distributions at different spatial locations",
      xlabel = "Location [km]",
      ylabel = "Vx [km/s]",
      xminorticksvisible = true)

   for i in 1:nxchunks
      hist!(ax, vx[i], normalization = :pdf, bins = 50,
         scale_to=-5000/nxchunks, offset=xmid[i], direction=:x)
   end

   means_str = [@sprintf "%d [km/s]" mean(vx[i]) for i in eachindex(vx)]
   text!(Point.(xmid, 300.0), text = means_str, align = (:right, :center),
      offset = (-60, 0), color = :black, fontsize=24)

   f
end

f = plot_dist(x, sols; nxchunks=4, ntchunks=100)
f = DisplayAs.PNG(f) #hide

# MHD states in SI units
n₁ = 1.0e6
T₁ = 720471.8506664868
Pi₁ = 0.0049735919716217296 * 1e-9
Pe₁ = 0.0049735919716217296 * 1e-9
Pth₁ = Pi₁ + Pe₁
V₁ = [-545.1484121835928, 0.0, 0.0] .* 1e3
B₁ = [5.0, 0.0, 0.0] .* 1e-9

n₂ = 3.6363636363636362e6
T₂ = 7.380333520264822e6
Pi₂ = 0.18526630094290936 * 1e-9
Pe₂ = 0.18526630094290936 * 1e-9
Pth₂ = Pi₂ + Pe₂
V₂ = [-149.91581335048804, 0.0, 0.0] .* 1e3
B₂ = [5.0, 0.0, 0.0] .* 1e-9;

E₁ = B₁ × V₁
E₂ = B₂ × V₂

# Thermal speed used to generate the Maxwellian distribution
Uth₁ = √(γ * Pth₁ / (n₁ * mᵢ))
Uth₂ = √(γ * Pth₂ / (n₂ * mᵢ))
# Shock normal direction range
x = range(-5_000e3, 5_000e3, length=100)
B = repeat(B₁, 1, length(x))
E = repeat(E₁, 1, length(x))
# Index for the shock location
mid_ = length(x) ÷ 2

B[:, 1:mid_] .= B₂
E[:, 1:mid_] .= E₂

vdf₁ = Maxwellian(V₁, Uth₁)
vdf₂ = Maxwellian(V₂, Uth₂)

trajectories = 2
weight₁ = n₁ / trajectories

prob = let
   # BC type 3 is Flat
   param = prepare(x, E, B; species=Proton, bc=3);
   stateinit = zeros(6) # particle position and velocity to be modified
   tspan = (0.0, 30.0)
   ODEProblem(trace!, stateinit, tspan, param)
end
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)

sols = solve(ensemble_prob, Vern9(), EnsembleSerial(); trajectories);

f = plot_traj(sols)
f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl