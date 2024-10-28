# ---
# title: Electron Fermi Acceleration Inside Foreshock Transient Cores
# id: demo_electron_acceleration_foreshock_transient
# date: 2024-10-27
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.11.1
# description: Fermi acceleration demonstration
# ---

# This example demonstrates electron acceleration via Fermi acceleration using a simple 1D model. It follows the setup described in [Fermi acceleration of electrons inside foreshock transient cores](https://doi.org/10.1002/2017JA024480).
# The simulation domain is $2\,R_E$ in x, where the initial magnetosheath locates at $x \in [0, 0.5] R_E$ with a magnetic field $B_z = 20\,\mathrm{nT}$ and the bow shock at $x = 0.5\,R_E$.
# The foreshock transient boundary is at $x = 1.5\,R_E$ initially and moves toward the bow shock at a speed $U=100\,\mathrm{km/s}$.
# The region beyond $x=1.5\,R_E$ is the foreshock transient sheath in which the magnetic field is $B_z=10\,\mathrm{nT}$. A convection electric field $E_y = -1\,\mathrm{mV/m}$, consistent with the velocity U, is introduced in the foreshock transient sheath.
# Between the two boundaries, $0.5\,R_E < x < 1.5\,R_E$, is the foreshock transient's core region.
#
# In the first case, the magnetic field in the core region is set to zero.
# In the second case, a magnetic fluctuation is imposed:
#
# ```math
# \begin{aligned}
# \delta B_{x,y,z} = \sum_{N = N_0}^{N_1} \delta B_N \cos\left( \frac{2\pi N x}{L_0} + \phi_{x,y,z}^N \right) \\
# \delta B_N = \tilde{B} (N/N_0)^{-1.2},\quad N = N_0, N_0 + 1, ..., N_1
# \end{aligned}
# ```
# Here $L_0 = 1\,R_E$ is the initial length of the core in the x direction.
# We choose $N_0 = 100, N_1 = 1000, \tilde{B}=0.2\,\mathrm{nT}$, and $\phi_{x,y,z}^N$ as the random phases of various modes between 0 and 2π (independently different in the x, y, and z directions).
# As the low-frequency wave speed ($\sim \mathcal{O}(10)\mathrm{km/s}$) is much smaller than the electron speed ($\sim \mathcal{O}(10^3)\mathrm{km/s}$), we do not include wave propagation in this 1D model.

#import DisplayAs #hide
using TestParticle
using TestParticle: kB, mₑ, Rₑ, qₑ
using OrdinaryDiffEq
using Random
using StaticArrays
using FHist
using CairoMakie, Printf
CairoMakie.activate!(type = "png") #hide

## For reproducible results
Random.seed!(1234)

## Analytic EM fields

function B(xu, t)
   Bz =
      if xu[1] < 0.5Rₑ
         20e-9
      elseif xu[1] > 1.5Rₑ - U*t
         10e-9
      else
         0.0
      end

   SA[0.0, 0.0, Bz]
end

function E(xu, t)
   Ey =
      if xu[1] > 1.5Rₑ - U*t
         -1e-3
      else
         0.0
      end

   SA[0.0, Ey, 0.0]
end

function B_field_perturb(xu, t)
   L₀, N₀, N₁, B̃ = 1Rₑ, 100, 1000, 0.2e-9

   δBx, δBy, δBz = 0.0, 0.0, 0.0
   for N in N₀:N₁
      δBn = B̃(N / N0)^-1.2
      ϕx, ϕy, ϕz = rand(3) .* 2π
      δBx += δBn * cos(2π * N * x / L₀ + ϕx)
      δBy += δBn * cos(2π * N * x / L₀ + ϕy)
      δBz += δBn * cos(2π * N * x / L₀ + ϕz)
   end

   SA[δBx, δBy, δBz]
end

function isoutofdomain(xv, p, t)
   if xv[1] < 0 || xv[1] > 2Rₑ
      return true
   else
      return false
   end
end

function prob_func(prob, i, repeat)
   x0 = [(0.5 + rand())*Rₑ, 0.0, 0.0] # launched in the core region
   u0 = [0.0, 0.0, 0.0]
   T₀ = 10 # [eV]
   vth = √(2T₀*abs(qₑ)/mₑ) # [m/s]
   vdf = Maxwellian(u0, vth)
   v0 = TestParticle.sample(vdf)

   prob = @views remake(prob, u0=[x0..., v0...])
end

"Kinetic energy."
get_kinetic_energy(dx, dy, dz) = 1 // 2 * (dx^2 + dy^2 + dz^2)

function plot_multiple(sol)
   energy = map(x -> get_kinetic_energy(x[4:6]...), sol.u) .* mₑ ./ abs(qₑ)
   ## Obtain the EM fields along the particle trajectory
   ## [mV/m]
   E = [sol.prob.p[2](sol[:,istep], sol.t[istep]).*1e3 for istep in eachindex(sol)]
   ## [nT]
   B = [sol.prob.p[3](sol[:,istep], sol.t[istep]).*1e9 for istep in eachindex(sol)]

   Ex = [e[1] for e in E]
   Ey = [e[2] for e in E]
   Ez = [e[3] for e in E]
   Bx = [b[1] for b in B]
   By = [b[2] for b in B]
   Bz = [b[3] for b in B]

   t = sol.t
   x = sol[1,:] ./ Rₑ
   y = sol[2,:] ./ Rₑ
   z = sol[3,:] ./ Rₑ

   fig = Figure(size = (900, 600), fontsize = 20)

   xlabels = ("", "", "", "t [s]")
   ylabels = ("KE [eV]", "Locations [RE]", "E [mV/m]", "B [nT]")
   limits = (
      (nothing, (nothing, nothing)),
      (nothing, (nothing, nothing)),
      (nothing, (nothing, nothing)),
      (nothing, (nothing, nothing)))
   
   axs = [Axis(fig[row, col], xlabel=xlabels[row], ylabel=ylabels[row], limits=limits[row])
      for row in eachindex(xlabels), col in 1:1]

   linkxaxes!(axs...)

   lines!(axs[1], t, energy)
   lines!(axs[2], t, x, label="x")
   lines!(axs[2], t, y, label="y")
   lines!(axs[2], t, z, label="z")
   lines!(axs[3], t, Ex, label="x")
   lines!(axs[3], t, Ey, label="y")
   lines!(axs[3], t, Ez, label="z")
   lines!(axs[4], t, Bx, label="x")
   lines!(axs[4], t, By, label="y")
   lines!(axs[4], t, Bz, label="z")
   
   for ax in @view axs[2:4]
      axislegend(ax, framevisible = false, orientation = :horizontal)
   end
   
   fig
end

function plot_dist_pairplots(sols; t=0, case=1)
   ##TODO: Optimization
   vx = Vector{eltype(sols[1].u[1])}(undef, 0)
   vy = similar(vx)
   vz = similar(vx)
   n = 0
   for sol in sols
      if (sol.t[end] ≥ t) && (1.5Rₑ - U*sol.t[end] > sol[1,end] > 0)
         n += 1
         v = sol(t)[4:6] ./ 1e3
         append!(vx, v[1])
         append!(vy, v[2])
         append!(vz, v[3])
      end
   end

   f = Figure(size=(650, 600), fontsize=18)
   h2d = Hist2D((vx, vy); nbins=(40, 40))
   _, _heatmap = plot(f[1,1], h2d;
      axis=(title="t = $t, case = $case, particle count = $(length(vx))",
      xlabel=L"V_x [km/s]", ylabel=L"V_y [km/s]", aspect=1, limits=(-1e4, 1e4, -1e4, 1e4)))

   Colorbar(f[1,2], _heatmap)

   f
end

function find_max_acceleration_index(sols; countall=true, tend=40)
   if countall
      ratio = [get_kinetic_energy(sol[4:6,end]...) / get_kinetic_energy(sol[4:6,1]...)
         for sol in sols]
   else
      ## only count the particles that are still trapped at t=tend
      ratio = [get_kinetic_energy(sol[4:6,end]...) / get_kinetic_energy(sol[4:6,1]...)
         for sol in sols if sol.t[end] > tend-0.1]
   end
   imax = argmax(ratio)

   energy_init = get_kinetic_energy(sols[imax][4:6,1]...) .* mₑ ./ abs(qₑ)
   energy_final = get_kinetic_energy(sols[imax][4:6,end]...) .* mₑ ./ abs(qₑ)
   @printf "Initial energy [eV]: %.2f " energy_init
   @printf "Final energy [eV]: %.2f " energy_final
   @printf "Kinetic energy change ratio: %.2f\n" ratio[imax]

   imax
end

const U = -100e3 # [m/s]

## Initial condition
stateinit = zeros(6)

## Time span [s]
tspan = (0, 40)

trajectories = 5000;

# Case 1: 0 core field

param = prepare(E, B; species=Electron);
prob = ODEProblem(trace!, stateinit, tspan, param) 
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)
## dt=2e-4, adaptive=false)
sols = solve(ensemble_prob, Vern9(), EnsembleThreads();
   isoutofdomain, trajectories, verbose=false)

## maximum acceleration ratio particle index
imax = find_max_acceleration_index(sols)

f = plot_multiple(sols[imax])
f = DisplayAs.PNG(f) #hide

f = plot_dist_pairplots(sols, t=tspan[1], case=1)
f = DisplayAs.PNG(f) #hide

f = plot_dist_pairplots(sols, t=tspan[2], case=1)
f = DisplayAs.PNG(f) #hide

# Case 2: B fluctuation core field

param = prepare(E, B_field_perturb; species=Electron);
prob = ODEProblem(trace!, stateinit, tspan, param) 
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)
## dt=2e-4, adaptive=false)
sols = solve(ensemble_prob, Vern9(), EnsembleThreads();
   isoutofdomain, trajectories, verbose=false)

## maximum acceleration ratio particle index
imax = find_max_acceleration_index(sols)

f = plot_multiple(sols[imax])
f = DisplayAs.PNG(f) #hide

f = plot_dist_pairplots(sols, t=tspan[1], case=1)
f = DisplayAs.PNG(f) #hide

f = plot_dist_pairplots(sols, t=tspan[2], case=1)
f = DisplayAs.PNG(f) #hide