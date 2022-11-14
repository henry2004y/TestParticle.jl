# Tracing charged particle in the Harris current sheet.
# This example shows how to trace protons in a stationary magnetic field that
# corresponds to the 1D Harris current sheet.
# Reference: https://en.wikipedia.org/wiki/Current_sheet
#
# Hongyang Zhou, hyzhou@umich.edu

using TestParticle
using TestParticle: getB_CS_harris
using OrdinaryDiffEq
using StaticArrays
using PyPlot
using Statistics: mean

## Obtain field

# Harris current sheet parameters in SI units
const B₀, L = 20e-9, 0.4TestParticle.Rₑ

function getB(xu)
   SVector{3}(getB_CS_harris(xu[1:3], B₀, L))
end

function getE(xu)
   SA[0.0, 0.0, 0.0]
end

## Initialize particles

m = TestParticle.mᵢ
q = TestParticle.qᵢ
c = TestParticle.c
Rₑ = TestParticle.Rₑ

Ek = 5e7 # [eV]

# initial velocity, [m/s]
v₀ = [c*√(1-1/(1+Ek*q/(m*c^2))^2), 0.0, 0.0]
# initial position, [m]
r₀ = [-5.0Rₑ, 0.0, 0.0]
stateinit = [r₀..., v₀...]

param = prepare(getE, getB)
tspan = (0.0, 10.0)

prob = ODEProblem(trace!, stateinit, tspan, param)

sol = solve(prob, Tsit5(); save_idxs=[1,2,3])

## Visualization

using3D()
fig = plt.figure(figsize=(10,6))
ax = fig.gca(projection="3d")

n = 100 # number of timepoints
ts = range(tspan..., length=n)
ax.plot(sol(ts, idxs=1)./Rₑ,sol(ts, idxs=2)./Rₑ, sol(ts, idxs=3)./Rₑ,
   label="50 MeV proton, B0 = 20 nT")

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
TestParticle.set_axes_equal(ax)

#savefig("test.png", bbox_inches="tight")