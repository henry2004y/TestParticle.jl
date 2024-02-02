# ---
# title: Magnetic mirror
# id: demo_magnetic_mirror
# date: 2023-04-20
# author: "[Tiancheng Liu](https://github.com/TCLiuu); [Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.9.0
# description: Charged particle in the magnetic mirror
# ---

# This example demonstrates the particle motion trajectory in a magnetic mirror
# and also illustrates the conservation of magnetic moment.
# From the third figure on the right, it can be seen that the zero-order quantity of
# magnetic moment is conserved, but its high-order part is not conserved under this definition
# of magnetic moment and oscillates rapidly. We can observe from this the oscillation
# characteristics of magnetic moments at different levels.
#
# This example is based on `demo_magneticbottle.jl`.

import DisplayAs #hide
using TestParticle
using TestParticle: getB_mirror
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
using Printf
using CairoMakie
CairoMakie.activate!(type = "png")

### Obtain field

## Magnetic mirror parameters in SI units
const I = 20. # current in the solenoid
const N = 45 # number of windings
const distance = 10. # distance between solenoids
const a = 4.0 # radius of each coil

function getB(xu)
   SVector{3}(getB_mirror(xu[1], xu[2], xu[3], distance, a, I*N))
end

function getE(xu)
   SA[0.0, 0.0, 0.0]
end

## velocity in the direction perpendicular to the magnetic field
function v_perp(xu)
   vu = @view xu[4:6]
   B = getB(xu)
   b = normalize(B)
   v_pa = (vu⋅b).*b
   return norm(vu - v_pa)
end

## magnetic field
absB(xu) = hypot(getB(xu)...)

## μ, magnetic moment
function mu(xu)
   return v_perp(xu)^2/hypot(getB(xu)...)
end

Et(xu) = hypot(xu[4:6]...)

### Initialize particles
m = TestParticle.mₑ
q = TestParticle.qₑ
c = TestParticle.c

## initial velocity, [m/s]
v₀ = [2.75, 2.5, 1.3] .* 0.001c  # confined
##v₀ = [0.25, 0.25, 5.9595] .* 0.01c  # escaped
## initial position, [m]
r₀ = [0.8, 0.8, 0.0]  # confined
##r₀ = [1.5, 1.5, 2.4]  # escaped
stateinit = [r₀..., v₀...]

param = prepare(getE, getB; species=Electron)
tspan = (0.0, 1e-4)

prob = ODEProblem(trace!, stateinit, tspan, param)

## Default Tsit5() and many solvers does not work in this case!
## AB4() has better performance in maintaining magnetic moment conservation compared to AB3().
sol_non = solve(prob, AB4(); dt=3e-9)

### Visualization

f = monitor(sol_non, vars=[absB, v_perp, mu])

## Plot coils
θ = range(0, 2π, length=100)
x = a.*cos.(θ)
y = a.*sin.(θ)
z = fill(distance/2, size(x))
ax = f[1:3,1:3]
lines!(ax, x, y, z, color=:red)
z = fill(-distance/2, size(x))
lines!(ax, x, y, z, color=:red)

## # The distribution of magnetic field along the z-axis or x-axis
## Bz(z) = hypot(getB(SA[0.0, 0.0, z])...)
## Bx(x) = hypot(getB(SA[x, 0.0, 0.5*distance])...)
## z = collect(-10:0.01:10)
## x = collect(-0.99*a:0.01:0.99*a)
## # Ba = Bz.(z)
## Ba = Bx.(x)
## # lines(z, Ba, color=:red)
## lines(x, Ba, color=:red)

f = DisplayAs.PNG(f) #hide