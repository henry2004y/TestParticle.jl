# Tracing charged particle in the Harris current sheet.
# This example shows how to trace protons in a stationary magnetic field that
# corresponds to the 1D Harris current sheet.
# Reference: https://en.wikipedia.org/wiki/Current_sheet
#
# Hongyang Zhou, hyzhou@umich.edu

using TestParticle
using DifferentialEquations
using PyPlot
using Statistics: mean

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


m = TestParticle.mᵢ
q = TestParticle.qᵢ
c = TestParticle.c
Rₑ = TestParticle.Rₑ

## Obtain field

# Harris current sheet parameters in SI units
const B₀, L = 20e-9, 0.4Rₑ

include("../utility/current_sheet.jl")

function getB(xu)
   getB_CS_harris(xu[1:3], B₀, L)
end

function getE(xu)
   [0.0, 0.0, 0.0]
end

## Initialize particles

Ek = 5e7 # [eV]

# initial velocity, [m/s]
v₀ = [c*sqrt(1-1/(1+Ek*q/(m*c^2))^2), 0.0, 0.0]
# initial position, [m]
r₀ = [-5.0Rₑ, 0.0, 0.0]
stateinit = [r₀..., v₀...]

param = prepare(getE, getB)
tspan = (0.0,10.0)

prob = ODEProblem(trace_analytic!, stateinit, tspan, param)

sol = solve(prob; save_idxs=[1,2,3], alg_hints=[:nonstiff])

## Visualization

fig = plt.figure(figsize=(10,6))
ax = fig.gca(projection="3d")

ax.plot(sol[1,:]./Rₑ,sol[2,:]./Rₑ, sol[3,:]./Rₑ, label="50 MeV proton, B0 = 20 nT")

ax.legend()
title("Particle trajectory near the Harris current sheet")
xlabel("x [Re]")
ylabel("y [Re]")
zlabel("z [Re]")

x = range(-8, 8, length=20)

X = collect(Float32, x)
Y = zeros(Float32, size(X)...)
Z = zeros(Float32, size(X)...)

B = zeros(Float32, 3, size(X)...)

i = 1
for (x,y) in zip(X,Y)
   B[1+3*(i-1):3*i] = getB_CS_harris([x,0.0,0.0], 4e-2, 1.0)
   global i += 1
end

for s = 1:3
   ax.quiver(X, Y, Z, vec(B[1,:,:]), vec(B[2,:,:]), vec(B[3,:,:]),
      color="k", alpha=0.6)
   @. Y -= 15
end

ax.set_box_aspect([2,5,1])
set_axes_equal(ax)

savefig("test.png", bbox_inches="tight")