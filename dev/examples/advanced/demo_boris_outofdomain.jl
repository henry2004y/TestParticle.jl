import DisplayAs #hide
using TestParticle, OrdinaryDiffEqVerner, StaticArrays
using TestParticle: qᵢ, mᵢ
using CairoMakie
CairoMakie.activate!(type = "png") #hide

"""
Set initial states.
"""
function prob_func(prob, i, repeat)
   prob = @views remake(prob; u0 = [prob.u0[1:3]..., 10.0 - i*2.0, prob.u0[5:6]...])
end

isoutofdomain(xv, p, t) =
   if isnan(xv[1])
      return true
   else
      return false
   end

# Number of cells for the field along each dimension
nx, ny = 4, 6
# Unit conversion factors between SI and dimensionless units
B₀ = 10e-9            # [T]
Ω = abs(qᵢ) * B₀ / mᵢ # [1/s]
t₀ = 1 / Ω            # [s]
U₀ = 1.0              # [m/s]
l₀ = U₀ * t₀          # [m]
E₀ = U₀*B₀            # [V/m]

x = range(0, 11, length = nx) # [l₀]
y = range(-21, 0, length = ny) # [l₀]

B = fill(0.0, 3, nx, ny) # [B₀]
B[3, :, :] .= 1.0

E(x) = SA[0.0, 0.0, 0.0] # [E₀]

# If bc == 1, we set a NaN value outside the domain (default);
# If bc == 2, we set periodic boundary conditions.
param = prepare(x, y, E, B; species = User, bc = 1);

# Initial conditions to be modified in prob_func
x0 = [0.0, 0.0, 0.0] # initial position [l₀]
u0 = [0.0, 0.0, 0.0] # initial velocity [v₀], will be overwritten in prob_func
stateinit = [x0..., u0...]
tspan = (0.0, 1.5π) # 3/4 gyroperiod

dt = 0.1
savestepinterval = 1
trajectories = 2
prob = TraceProblem(stateinit, tspan, param; prob_func)

sols = TestParticle.solve(prob; dt, savestepinterval, isoutofdomain, trajectories)

f = Figure(fontsize = 18)
ax = Axis(f[1, 1],
   title = "Proton trajectory",
   xlabel = "X",
   ylabel = "Y",
   limits = (-10.1, 10.1, -20.1, 0.1),
   aspect = DataAspect()
)

for i in eachindex(sols)
   lines!(ax, sols[i]; idxs = (1, 2), label = string(i), color = Makie.wong_colors()[i])
end

axislegend(position = :lt, framevisible = false)

f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
