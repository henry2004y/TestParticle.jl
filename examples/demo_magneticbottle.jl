# Tracing charged particle in the magnetic bottle.
# This example shows how to trace protons in a stationary magnetic field that
# corresponds to the magnetic bottle.
# Reference: https://en.wikipedia.org/wiki/Magnetic_mirror#Magnetic_bottles
#
# Hongyang Zhou, hyzhou@umich.edu

using TestParticle
using TestParticle.MagneticConfinement: getB_Bottle
using OrdinaryDiffEq
using StaticArrays
using PyPlot
using Statistics: mean
using Printf

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
   SVector{3}(getB_Bottle(xu[1], xu[2], xu[3], distance, a, b, I1*N1, I2*N2))
end

function getE(xu)
   SA[0.0, 0.0, 0.0]
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

param = prepare(getE, getB; species=Electron)
tspan = (0.0, 3e-7)

prob = ODEProblem(trace_analytic_relativistic!, stateinit, tspan, param)

@printf "Speed = %6.4f %s\n" √(v₀[1]^2+v₀[2]^2+v₀[3]^2)/c*100 "% speed of light"
@printf "Energy = %6.4f MeV\n" (1/√(1-(v₀[1]/c)^2-(v₀[2]/c)^2-(v₀[3]/c)^2)-1)*m*c^2/abs(q)/1e6

# Default Tsit5() and many solvers does not work in this case!
sol = solve(prob, AB3(); dt=1.5e-12, save_idxs=[1,2,3])

## Visualization

using3D()
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

TestParticle.set_axes_equal(ax)

#savefig("test.png", bbox_inches="tight")