import DisplayAs #hide
using TestParticle
using TestParticle: qᵢ, mᵢ
using StaticArrays
using OrdinaryDiffEq
using CairoMakie
CairoMakie.activate!(type = "png")

uniform_B(x) = SA[0.0, 0.0, 0.01]
uniform_E(x) = SA[0.0, 0.0, 0.0]

"Set initial states."
function prob_func(prob, i, repeat)
   prob.u0[4] = 10.0 - i*2.0

   prob
end

function isoutofdomain(xv)
   if isnan(xv[1])
      return true
   else
      return false
   end
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

x = range(0, 11, length=nx) # [l₀]
y = range(-21, 0, length=ny) # [l₀]

B = fill(0.0, 3, nx, ny) # [B₀]
B[3,:,:] .= 1.0

E(x) = SA[0.0, 0.0, 0.0] # [E₀]

# If bc == 1, we set a NaN value outside the domain (default);
# If bc == 2, we set periodic boundary conditions.
param = prepare(x, y, E, B; species=User, bc=1);

# Initial conditions to be modified in prob_func
x0 = [0.0, 0.0, 0.0] # initial position [l₀]
u0 = [0.0, 0.0, 0.0] # initial velocity [v₀]
stateinit = [x0..., u0...]
tspan = (0.0, 1.5π) # 3/4 gyroperiod

paramBoris = BorisMethod(param)
dt = 0.1
savestepinterval = 1
trajectories = 2
prob = TraceProblem(stateinit, tspan, dt, paramBoris; prob_func)

sols = trace_trajectory(prob; savestepinterval, isoutofdomain, trajectories)

f = Figure(fontsize = 18)
ax = Axis(f[1, 1],
   title = "Proton trajectory",
   xlabel = "X",
   ylabel = "Y",
   limits = (-10.1, 10.1, -20.1, 0.1),
   aspect = DataAspect()
)

for i in eachindex(sols)
   @views lines!(ax, sols[i].u[1,:], sols[i].u[2,:], label=string(i))
end

axislegend()

f = DisplayAs.PNG(f) #hide

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
