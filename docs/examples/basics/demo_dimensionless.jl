# ---
# title: Dimensionless Units
# id: demo_dimensionless
# date: 2023-12-14
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.4
# description: Tracing charged particle in dimensionless units
# ---

# This example shows how to trace charged particles in dimensionless units.
# Let ``B_0`` be the reference magnetic field strength and ``U_0`` the reference velocity. The Lorentz equation without the electric field can be written as
#
# ```math
# \frac{\mathrm{d}\mathbf{u}}{\mathrm{d}t} = \Omega \mathbf{u}\times\mathbf{b}
# ```
#
# where ``\Omega = qB_0/m`` is the reference gyrofrequency and  ``\mathbf{b} = \mathbf{B}/B_0`` is the normalized magnetic field.
#
# Then we can normalize the time scale to the gyroperiod ``T=2\pi/\Omega``. We use the gyroradius for the spatial scale, which is proportional to the velocity scale
#
# ```math
# r_L = \frac{U_0}{\Omega}
# ```
#
# For instance, if the magnetic field is homogeneous and the initial perpendicular velocity is 4 (in the normalized units), then the gyroradius is 4.
# Tracing in dimensionless units is beneficial for many scenarios. For example, MHD simulations do not have intrinsic scales. Therefore, we can do dimensionless particle tracing in MHD fields, and then convert to any scale we would like.
#
# Now let's demonstrate this with `trace_normalized!`.

import DisplayAs #hide
using TestParticle
using OrdinaryDiffEq
using Meshes

using CairoMakie
CairoMakie.activate!(type = "png")

x = range(-10, 10, length=15)
y = range(-10, 10, length=20)
z = range(-10, 10, length=25)
B = fill(0.0, 3, length(x), length(y), length(z)) # [B₀]
E = fill(0.0, 3, length(x), length(y), length(z)) # [E₀]

B₀ = 10e-9

B[3,:,:,:] .= 1.0

Δx = x[2] - x[1]
Δy = y[2] - y[1]
Δz = z[2] - z[1]

mesh = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
   (x[1], y[1], z[1]),
   (Δx, Δy, Δz))

param = prepare(mesh, E, B, B₀; species=Proton)

Ω = param[1]
U₀ = 1.0

x0 = [0.0, 0.0, 0.0] # initial position [l₀]
u0 = [4*Ω*U₀, 0.0, 0.0] # initial velocity [v₀]
stateinit = [x0..., u0...]
tspan = (0.0, 2π/Ω) # one gyroperiod

prob = ODEProblem(trace_normalized!, stateinit, tspan, param)
sol = solve(prob, Vern9())

### Visualization
f = Figure(fontsize = 18)
ax = Axis(f[1, 1],
   title = "Proton trajectory",
   xlabel = "X",
   ylabel = "Y",
   limits = (-4.1, 4.1, -8.1, 0.1),
   aspect = DataAspect()
)

lines!(ax, sol, vars=(1,2))

f = DisplayAs.PNG(f) #hide