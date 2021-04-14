# Tracing charged particle in Tokamak.
# This example shows how to trace protons in a stationary magnetic field that
# corresponds to a Tokamak.
#
# Hongyang Zhou, hyzhou@umich.edu

using TestParticle
import TestParticle.MagneticConfinement: getB_Tokamak
using DifferentialEquations
using PyPlot
using Statistics: mean
using Printf

"""
    set_axes_equal(ax)

Set 3D plot axes to equal scale.
Make axes of 3D plot have equal scale so that spheres appear as spheres and
cubes as cubes. Required since `ax.axis('equal')` and `ax.set_aspect('equal')`
don't work on 3D.
"""
function set_axes_equal(ax)
   limits = zeros(2,3)
   limits[:,1] .= ax.get_xlim3d()
   limits[:,2] .= ax.get_ylim3d()
   limits[:,3] .= ax.get_zlim3d()
   origin = mean(limits, dims=1)
   radius = @. 0.5 * max(abs(limits[2,:] - limits[1,:]))
   _set_axes_radius(ax, origin, radius)
end

function _set_axes_radius(ax, origin, radius)
   x, y, z = origin
   ax.set_xlim3d([x - radius[1], x + radius[1]])
   ax.set_ylim3d([y - radius[2], y + radius[2]])
   ax.set_zlim3d([z - radius[3], z + radius[3]])
end

## Obtain field

# Magnetic bottle parameters in SI units
const ICoil = 80. # current in the coil
const N = 15000 # number of windings
const IPlasma = 1e6 # current in the plasma
const a = 1.5 # radius of each coil
const b = 0.8 # radius of central region

function getB(xu)
   getB_Tokamak(xu[1], xu[2], xu[3], a, b, ICoil*N, IPlasma)
end

function getE(xu)
   [0.0, 0.0, 0.0]
end

## Initialize particles
m = TestParticle.mᵢ
q = TestParticle.qᵢ
c = TestParticle.c

# initial velocity, [m/s]
v₀ = [-0.1, -0.15, 0.0] .* c
# initial position, [m]
r₀ = [2.3, 0.0, 0.0]
stateinit = [r₀..., v₀...]

param = prepare(getE, getB; species=Proton)
tspan = (0.0, 1e-6)

prob = ODEProblem(trace_analytic_relativistic!, stateinit, tspan, param)

@printf "Speed = %6.4f %s\n" √(v₀[1]^2+v₀[2]^2+v₀[3]^2)/c*100 "% speed of light"
@printf "Energy = %6.4f MeV\n" (1/√(1-(v₀[1]/c)^2-(v₀[2]/c)^2-(v₀[3]/c)^2)-1)*m*c^2/abs(q)/1e6

# Default Tsit5() does not work in this case!
sol = solve(prob, AB3(); dt=2e-11, save_idxs=[1,2,3])

## Visualization

using3D()
fig = plt.figure(figsize=(10,6))
ax = fig.gca(projection="3d")

ax.plot(sol[1,:], sol[2,:], sol[3,:], label="proton")

# Plot coils
θ = range(0, 2π, length=201)
y = a * cos.(θ)
z = a * sin.(θ)
for i in 0:17
   ϕ = i*π/9
   ax.plot(y*sin(ϕ) .+ (a+b)*sin(ϕ), y*cos(ϕ) .+ (a+b)*cos(ϕ), z, color="k")
end

# Plot Tokamak
u = range(0, 2π, length=100)
v = range(0, 2π, length=100)

U = [y for _ in u, y in v]
V = [x for x in u, _ in v]

X = @. (a + b + (a - 0.05)*cos(U)) * cos(V)
Y = @. (a + b + (a - 0.05)*cos(U)) * sin(V)
Z = @. (a - 0.05) * sin(U)

ax.plot_surface(X, Y, Z, alpha=0.1)

ax.legend()
title("Particle trajectory in Tokamak")
xlabel("x [m]")
ylabel("y [m]")
zlabel("z [m]")

xMin, xMax = -4, 4
yMin, yMax = -4, 4
zMin, zMax = -4, 4

xlim(xMin, xMax)
ylim(yMin, yMax)
zlim(zMin, zMax)

set_axes_equal(ax)

#savefig("test.png", bbox_inches="tight")