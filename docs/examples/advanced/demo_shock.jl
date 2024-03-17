# ---
# title: Shock
# id: demo_shock
# date: 2024-02-21
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.10.2
# description: Tracing charged particle in an MHD plane shock
# ---

# This example shows how to trace protons of a certain energy in MHD plane shocks.

import DisplayAs #hide
using TestParticle, OrdinaryDiffEq
using TestParticle: mᵢ, kB
using LinearAlgebra
using Statistics: mean, std
using Printf
using Random
using CairoMakie, PairPlots
CairoMakie.activate!(type = "png") #hide

## For reproducible results
Random.seed!(1234)

"Set initial conditions."
function prob_func(prob, i, repeat)
   v₀ = sample(vdf₁)
   r₀ = [5_000e3, 0.0, 0.0]

   prob = remake(prob; u0 = [r₀..., v₀...])
end

# Perpendicular shock is a special shock in which both the upstream and downstream plasma flows are perpendicular to the magnetic field, as well as the shock front.
## MHD states in SI units
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

# In ideal MHD, the electric field simply contains the convection term.
# For perpendicular shocks, the electric field across the shock is continuous.
# In the upstream, particles follow straight lines because it's purely an ExB drift; across the shock, ExB drift changes to the downstream bulk velocity.

E₁ = B₁ × V₁
E₂ = B₂ × V₂

## Shock normal direction range
x = range(-5_000e3, 5_000e3, length=100)
B = repeat(B₁, 1, length(x))
E = repeat(E₁, 1, length(x))
## Index for the shock location
mid_ = length(x) ÷ 2

B[:, 1:mid_] .= B₂
E[:, 1:mid_] .= E₂

const vdf₁ = Maxwellian(V₁, Pth₁, n₁; m=mᵢ)
vdf₂ = Maxwellian(V₂, Pth₂, n₂; m=mᵢ)

trajectories = 400
weight₁ = n₁ / trajectories # relation between test particle and real particles

prob = let
   ## BC type 3 is Flat
   param = prepare(x, E, B; species=Proton, bc=3);
   stateinit = zeros(6) # particle position and velocity to be modified
   tspan = (0.0, 30.0)
   ODEProblem(trace!, stateinit, tspan, param)
end
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)

sols = solve(ensemble_prob, Vern9(), EnsembleSerial(); trajectories);

# Sample particle trajectories

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

# Phase space distributions

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

   v̄x = mean.(vx)
   vth = [std(vx[i]; corrected=false, mean=v̄x[i]) for i in 1:nxchunks]
   means_str = [@sprintf "Vx: %d [km/s]" v̄x[i] for i in eachindex(v̄x)]
   std_str = [@sprintf "Vth: %d [km/s]" vth[i] for i in eachindex(vth)]
   text!(Point.(xmid.+400, 300.0), text = means_str, align = (:right, :center),
      offset = (-60, 0), color = :black, fontsize=24)
   text!(Point.(xmid.+400, -700.0), text = std_str, align = (:right, :center),
      offset = (-60, 0), color = :black, fontsize=24)

   f
end

f = plot_dist(x, sols; nxchunks=4, ntchunks=100)
f = DisplayAs.PNG(f) #hide

# Even with 400 particles, we are still able to statistically approximate the velocity moment downstream of the perpendicular shock.
# While the upstream thermal speed is close to what we set, the downstream thermal speed is higher than our precalculated value. (Why?)

println("Vth₁ = ", round(vdf₁.vth / 1e3, digits=2), " km/s") #hide
println("Vth₂ = ", round(vdf₂.vth / 1e3, digits=2), " km/s") #hide

# A nice way to present the 3d distributions in 2d is via the pair plots:

function plot_dist_pairplot(x, sols; ntchunks::Int=20)
   nxchunks = 2
   trange = range(sols[1].prob.tspan..., length=ntchunks)

   xrange = range(x[1], x[end], length=nxchunks+1)
   dx = (x[end] - x[1]) / nxchunks
   xmid = range(x[1] + 0.5dx, x[end] - 0.5dx, length=nxchunks) ./ 1e3

   table = [(;
   vx = Float64[],
   vy = Float64[],
   vz = Float64[],
   ) for _ in 1:nxchunks]

   for sol in sols
      for t in trange
         xv = sol(t)
         for i in 1:nxchunks
            if xrange[i] < xv[1] ≤ xrange[i+1]
               push!(table[i].vx, xv[4] / 1e3)
               push!(table[i].vy, xv[5] / 1e3)
               push!(table[i].vz, xv[6] / 1e3)
            end
         end
      end
   end

   for i in eachindex(table)
      if isempty(v[i])
         push!(table[i].vx, 0.0)
         push!(table[i].vy, 0.0)
         push!(table[i].vz, 0.0)
      end
   end

   f = Figure(size = (800, 800), fontsize=18)

   c1 = Makie.wong_colors(0.5)[1]
   c2 = Makie.wong_colors(0.5)[2]

   l1 = @sprintf "x: %d [km] downstream" xmid[1]
   l2 = @sprintf "x: %d [km] upstream" xmid[2]

   pairplot(f[1,1],
       PairPlots.Series(table[1], label=l1, color=c1, strokecolor=c1),
       PairPlots.Series(table[2], label=l2, color=c2, strokecolor=c2),
   )

   f
end

f = plot_dist_pairplot(x, sols; ntchunks=20)
f = DisplayAs.PNG(f) #hide

# For parallel shocks, we have a different scenario.
# Parallel shock is another special shock in which both the upstream and downstream plasma flows are parallel to the magnetic field, as well as perpendicular to the shock front.
## MHD states in SI units
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

# There is no convection electric field in the parallel shock:

E₁ = B₁ × V₁
E₂ = B₂ × V₂

# Therefore, when we trace particles, there is no deceleration across the shock:

## Shock normal direction range
x = range(-5_000e3, 5_000e3, length=100)
B = repeat(B₁, 1, length(x))
E = repeat(E₁, 1, length(x))
## Index for the shock location
mid_ = length(x) ÷ 2

B[:, 1:mid_] .= B₂
E[:, 1:mid_] .= E₂

vdf₁ = Maxwellian(V₁, Pth₁, n₁; m=mᵢ)
vdf₂ = Maxwellian(V₂, Pth₂, n₂; m=mᵢ)

trajectories = 2
weight₁ = n₁ / trajectories

prob = let
   ## BC type 3 is Flat
   param = prepare(x, E, B; species=Proton, bc=3);
   stateinit = zeros(6) # particle position and velocity to be modified
   tspan = (0.0, 30.0)
   ODEProblem(trace!, stateinit, tspan, param)
end
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)

sols = solve(ensemble_prob, Vern9(), EnsembleSerial(); trajectories);

f = plot_traj(sols)
f = DisplayAs.PNG(f) #hide

# Clearly, test particle tracing in MHD parallel shocks fails to recover physics.
# MHD parallel shocks are essentially hydrodynamic shocks where magnetic field plays no role.
# Due to the lack of collision and other diffusion processes, we are unable to capture the correct microscopic scenario here.