# Tracing charged particle in Tokamak.
# This example shows how to trace protons in a stationary magnetic field that
# corresponds to a Tokamak.
#
# Hongyang Zhou, hyzhou@umich.edu

using TestParticle
using TestParticle: getB_tokamak
using OrdinaryDiffEq
using StaticArrays
using PyPlot
using Statistics: mean
using Printf

## Obtain field

# Magnetic bottle parameters in SI units
const ICoil = 80. # current in the coil
const N = 15000 # number of windings
const IPlasma = 1e6 # current in the plasma
const a = 1.5 # radius of each coil
const b = 0.8 # radius of central region

function getB(xu)
   SVector{3}(getB_tokamak(xu[1], xu[2], xu[3], a, b, ICoil*N, IPlasma))
end

function getE(xu)
   SA[0.0, 0.0, 0.0]
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

# Default Tsit5() alone does not work in this case! You need to set a maximum
# timestep to maintain stability, or choose a different algorithm as well.
# The sample figure in the gallery is generated with AB3() and dt=2e-11.
sol = solve(prob, Tsit5(); dt=2e-11, save_idxs=[1,2,3])

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

TestParticle.set_axes_equal(ax)

#savefig("test.png", bbox_inches="tight")