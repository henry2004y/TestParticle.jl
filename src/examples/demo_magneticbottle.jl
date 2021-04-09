# Tracing charged particle in the magnetic bottle.
# This example shows how to trace protons in a stationary magnetic field that
# corresponds to the magnetic bottle.
# Reference: https://en.wikipedia.org/wiki/Magnetic_mirror#Magnetic_bottles
#
# Hongyang Zhou, hyzhou@umich.edu

using TestParticle
import TestParticle.MagneticConfinement: getB_Bottle
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
const I1 = 80. # current in the solenoid
const N1 = 1500 # number of windings
const I2 = 100. # current in the central solenoid
const N2 = 8000 # number of windings
const distance = 10. # distance between solenoids
const a = 1.0 # radius of each coil
const b = 4.0 # radius of central coil


function getB(xu)
   getB_Bottle(xu[1], xu[2], xu[3], distance, a, b, I1*N1, I2*N2)
end

function getE(xu)
   [0.0, 0.0, 0.0]
end

## Initialize particles
m = TestParticle.mₑ
q = TestParticle.qₑ
c = TestParticle.c

# initial velocity, [m/s]
v₀ = [0.8, 0.48, 0.3595] .* c
# initial position, [m]
r₀ = [0.0, 0.0, -0.8]
stateinit = [r₀..., v₀...]

# Theoretically we can take advantage of the fact that magnetic field does not
# accelerate particles, so that γ remains constant. However, we are not doing
# that here since it is not generally true in the EM field.

param = prepare(getE, getB; species="electron")
tspan = (0.0, 3e-7)

prob = ODEProblem(trace_analytic_relativistic!, stateinit, tspan, param)

@printf "Speed = %6.4f %s\n" √(v₀[1]^2+v₀[2]^2+v₀[3]^2)/c*100 "% speed of light"
@printf "Energy = %6.4f MeV\n" (1/√(1-(v₀[1]/c)^2-(v₀[2]/c)^2-(v₀[3]/c)^2)-1)*m*c^2/abs(q)/1e6

# Default Tsit5() does not work in this case!
sol = solve(prob, AB3(); dt=1.5e-12, save_idxs=[1,2,3])

## Visualization

fig = plt.figure(figsize=(10,6))
ax = fig.gca(projection="3d")

ax.plot(sol[1,:], sol[2,:], sol[3,:], label="electron")

# Plot coils
θ = range(0, 2π, length=100)
x = a.*cos.(θ)
y = a.*sin.(θ)
z = fill(distance/2, size(x))
ax.plot(x, y, z, color="r", label="coil")

z = fill(-distance/2, size(x))
ax.plot(x, y, z, color="r")

x = b.*cos.(θ)
y = b.*sin.(θ)
z = fill(0.0, size(x))
ax.plot(x, y, z, color="r")

ax.legend()
title("Electron trajectory in the magnetic bottle")
xlabel("x [m]")
ylabel("y [m]")
zlabel("z [m]")

xMin, xMax = -5, 5
yMin, yMax = -5, 5
zMin, zMax = -5.2, 5.2

xlim(xMin, xMax)
ylim(yMin, yMax)
zlim(zMin, zMax)

set_axes_equal(ax)

#savefig("test.png", bbox_inches="tight")