import DisplayAs #hide
using TestParticle
using TestParticle: getB_dipole, getE_dipole, sph2cart, mᵢ, qᵢ, c, Rₑ
using OrdinaryDiffEq
using StaticArrays
using FieldTracer
using CairoMakie
CairoMakie.activate!(type = "png")

function getB_superposition(xu)
   getB_dipole(xu) + SA[0.0, 0.0, -10e-9]
end

"Boundary condition check."
function isoutofdomain(u, p, t)
   rout = 18Rₑ
   if hypot(u[1], u[2], u[3]) < 1.1Rₑ ||
      abs(u[1]) > rout || abs(u[2]) > rout || abs(u[3]) > rout
      return true
   else
      return false
   end
end

"Set initial conditions."
function prob_func(prob, i, repeat)
   # initial particle energy
   Ek = 5e3 # [eV]
   # initial velocity, [m/s]
   v₀ = sph2cart(c*sqrt(1-1/(1+Ek*qᵢ/(mᵢ*c^2))^2), 0.0, π/4)
   # initial position, [m]
   r₀ = sph2cart(13*Rₑ, π*i, π/2)
   prob.u0 .= [r₀..., v₀...]

   prob
end

# obtain field
param = prepare(getE_dipole, getB_superposition)
stateinit = zeros(6) # particle position and velocity to be modified
tspan = (0.0, 2000.0)
trajectories = 2

prob = ODEProblem(trace!, stateinit, tspan, param)
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)

# See https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/
# for the solver options
sols = solve(ensemble_prob, Vern9(), EnsembleSerial(); reltol=1e-5,
   trajectories, isoutofdomain, dense=true, save_on=true)

### Visualization

f = Figure()
ax = Axis3(f[1, 1],
   title = "5 keV Protons in a vacuum superposition magnetosphere",
   xlabel = "x [Re]",
   ylabel = "y [Re]",
   zlabel = "z [Re]",
   aspect = :data,
)

invRE = 1 / Rₑ
for sol in sols
   l = lines!(ax, sol)
   scale!(l, invRE, invRE, invRE)
end

# Field lines
function get_numerical_field(x, y, z)
   bx = zeros(length(x), length(y), length(z))
   by = similar(bx)
   bz = similar(bx)

   for i in CartesianIndices(bx)
      pos = [x[i[1]], y[i[2]], z[i[3]]]
      bx[i], by[i], bz[i] = getB_superposition(pos)
   end

   bx, by, bz
end

function trace_field!(ax, x, y, z, unitscale)
   bx, by, bz = get_numerical_field(x, y, z)

   zs = 0.0
   nr, nϕ = 8, 4
   dϕ = 2π / nϕ
   for r in range(8Rₑ, 16Rₑ, length=nr), ϕ in range(0, 2π-dϕ, length=nϕ)
      xs = r * cos(ϕ)
      ys = r * sin(ϕ)

      x1, y1, z1 = FieldTracer.trace(bx, by, bz, xs, ys, zs, x, y, z; ds=0.1, maxstep=10000)

      lines!(ax, x1.*unitscale, y1.*unitscale, z1.*unitscale, color=:gray)
   end
end

x = range(-18Rₑ, 18Rₑ, length=50)
y = range(-18Rₑ, 18Rₑ, length=50)
z = range(-18Rₑ, 18Rₑ, length=50)

trace_field!(ax, x, y, z, invRE)

f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl