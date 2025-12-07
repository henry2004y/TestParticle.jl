# # Reconnection X-line
#
# This example traces protons and electrons near a magnetic reconnection X-line.
# The magnetic field is modeled as a Harris current sheet with a superposed
# X-type perturbation.

import DisplayAs #hide
using TestParticle, OrdinaryDiffEq, StaticArrays
import TestParticle as TP
using TestParticle: Rₑ
using LinearAlgebra: norm
using CairoMakie
CairoMakie.activate!(type = "png") #hide

### Define the Field

## Parameters
const B₀ = 20e-9  # Asymptotic magnetic field [T]
const Bₙ = 1e-9   # Normal magnetic field gradient parameter [T]
const L = 0.5Rₑ   # Current sheet thickness [m]
const Eᵣ = 0.1e-3 # Reconnection electric field [V/m]

function getB(xu)
    x, y, z = xu[1], xu[2], xu[3]
    return SVector{3}(B₀ * tanh(z / L), 0.0, Bₙ * x / L)
end

function getE(xu)
    return SVector{3}(0.0, Eᵣ, 0.0)
end

### Trace Particles

## Initialize proton
x0_p = [-1.0Rₑ, 0.0, 0.1Rₑ]
Ek_p = 1e3 # eV
v_p = TP.get_velocity(Ek_p, TP.mᵢ, 1) # magnitude
v0_p = [v_p, 0.0, 0.0]

stateinit_p = [x0_p..., v0_p...]
param_p = prepare(getE, getB, species=TP.Proton)
tspan_p = (0.0, 20.0) # seconds

prob_p = ODEProblem(TP.trace!, stateinit_p, tspan_p, param_p)
sol_p = solve(prob_p, Vern9())

## Initialize electron
x0_e = [-1.0Rₑ, 0.0, 0.1Rₑ]
Ek_e = 1000 # eV
v_e = TP.get_velocity(Ek_e, TP.mₑ, 1) # magnitude
v0_e = [v_e, 0.0, 0.0]

stateinit_e = [x0_e..., v0_e...]
param_e = prepare(getE, getB, species=TP.Electron)
tspan_e = (0.0, 2.0) # seconds, electron is much faster

prob_e = ODEProblem(TP.trace!, stateinit_e, tspan_e, param_e)
sol_e = solve(prob_e, Vern9())

### Visualization

f = Figure(size=(1000, 500), fontsize = 18)

# Plot Proton
ax1 = Axis3(f[1, 1],
   title = "Proton trajectory",
   xlabel = "x [Re]",
   ylabel = "y [Re]",
   zlabel = "z [Re]",
   aspect = :data
)

n = 2000 # number of timepoints
ts_p = range(tspan_p..., length = n)
x_p = sol_p(ts_p, idxs = 1) ./ Rₑ |> Vector
y_p = sol_p(ts_p, idxs = 2) ./ Rₑ |> Vector
z_p = sol_p(ts_p, idxs = 3) ./ Rₑ |> Vector

lines!(ax1, x_p, y_p, z_p, label = "1 keV Proton", color = :red)
axislegend(ax1)

# Plot Electron
ax2 = Axis3(f[1, 2],
   title = "Electron trajectory",
   xlabel = "x [Re]",
   ylabel = "y [Re]",
   zlabel = "z [Re]",
   aspect = :data
)

ts_e = range(tspan_e..., length = n)
x_e = sol_e(ts_e, idxs = 1) ./ Rₑ |> Vector
y_e = sol_e(ts_e, idxs = 2) ./ Rₑ |> Vector
z_e = sol_e(ts_e, idxs = 3) ./ Rₑ |> Vector

lines!(ax2, x_e, y_e, z_e, label = "1 keV Electron", color = :blue)
axislegend(ax2)

# Add B field arrows to proton plot for context
function plot_B!(ax)
   xrange = range(-2, 2, length = 5)
   yrange = range(-1, 1, length = 3)
   zrange = range(-1, 1, length = 5)

   ps = [Point3f(x, y, z) for x in xrange for y in yrange for z in zrange]
   B = map(p -> Vec3f(getB(p .* Rₑ) ./ B₀), ps)
   Bmag = norm.(B)

   arrows!(ax, ps, B, fxaa = true, color = Bmag, lengthscale = 0.4, arrowsize = 0.05)
end

plot_B!(ax1)
plot_B!(ax2)

f = DisplayAs.PNG(f) #hide
