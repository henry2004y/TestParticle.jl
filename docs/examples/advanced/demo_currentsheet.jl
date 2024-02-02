# ---
# title: Current sheet
# id: demo_currentsheet
# date: 2023-04-20
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.0
# description: Tracing charged particle in the Harris current sheet
# ---

# This example shows how to trace protons in a stationary magnetic field that corresponds to
# the 1D Harris current sheet defined by a reference strength and width.
# Reference: https://en.wikipedia.org/wiki/Current_sheet

import DisplayAs #hide
using TestParticle
using TestParticle: getB_CS_harris
using OrdinaryDiffEq
using StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png")

### Obtain field

## Harris current sheet parameters in SI units
const B₀, L = 20e-9, 0.4TestParticle.Rₑ

function getB(xu)
   SVector{3}(getB_CS_harris(xu[1:3], B₀, L))
end

function getE(xu)
   SA[0.0, 0.0, 0.0]
end

### Initialize particles

m = TestParticle.mᵢ
q = TestParticle.qᵢ
c = TestParticle.c
Rₑ = TestParticle.Rₑ
## initial particle energy, [eV]
Ek = 5e7
## initial velocity, [m/s]
v₀ = [c*√(1-1/(1+Ek*q/(m*c^2))^2), 0.0, 0.0]
## initial position, [m]
r₀ = [-5.0Rₑ, 0.0, 0.0]
stateinit = [r₀..., v₀...]

param = prepare(getE, getB)
tspan = (0.0, 10.0)

prob = ODEProblem(trace!, stateinit, tspan, param)

sol = solve(prob, Tsit5(); save_idxs=[1,2,3])

### Visualization

f = Figure()
ax = Axis3(f[1, 1],
   title = "Particle trajectory near the Harris current sheet",
   xlabel = "x [Re]",
   ylabel = "y [Re]",
   zlabel = "z [Re]",
   aspect = :data,
)

n = 100 # number of timepoints
ts = range(tspan..., length=n)
x = sol(ts, idxs=1)./Rₑ |> Vector
y = sol(ts, idxs=2)./Rₑ |> Vector
z = sol(ts, idxs=3)./Rₑ |> Vector

l = lines!(ax, x, y, z, label="50 MeV proton, B0 = 20 nT")
axislegend()

X, Y, Z = let xrange = range(-8, 8, length=20)
   X = collect(Float32, xrange)
   Y = zeros(Float32, size(X)...)
   Z = zeros(Float32, size(X)...)
   X, Y, Z
end

B = zeros(Float32, 3, size(X)...)

i = 1
for (x,y) in zip(X,Y)
   B[1+3*(i-1):3*i] = getB_CS_harris([x,0.0,0.0], 4e-2, 1.0)
   global i += 1
end

for s = 1:3
   quiver!(ax, X, Y, Z, vec(B[1,:,:]), vec(B[2,:,:]), vec(B[3,:,:]),
      color=:black, alpha=0.6, lengthscale=100.0)
   @. Y -= 15
end

f = DisplayAs.PNG(f) #hide