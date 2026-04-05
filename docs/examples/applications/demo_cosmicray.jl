import DisplayAs #hide
using TestParticle, OrdinaryDiffEqVerner, StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png") #hide

## Number of cells for the field along each dimension
nx, ny, nz = 64, 64, 1
## Unit conversion factors for length
rL0 = 1.0
L = 4.0
## Set length scales
x = range(-L / 2, L / 2, length = nx) # [rL0]
y = range(-L / 2, L / 2, length = ny) # [rL0]
z = range(-1.0, 1.0, length = nz) # [rL0]

B = fill(0.0, 3, nx, ny, nz) # [B0]
B[3, :, :, :] .= 1.0
E(x) = SA[0.0, 0.0, 0.0] # [E₀]
## periodic bc = 2
param = prepare(x, y, z, E, B; m = 1, q = 1, bc = 2)

## Initial condition
stateinit = let
    x0 = [0.0, 0.0, 0.0] # initial position [l₀]
    u0 = [1.0, 0.0, 0.0] # initial velocity [v₀] -> r = 1 * rL0
    [x0..., u0...]
end
## Time span
tspan = (0.0, 100.0) # half gyroperiod

prob = ODEProblem(trace_normalized!, stateinit, tspan, param)
sol = solve(prob, Vern9())

### Visualization
f = Figure(fontsize = 18)
ax = Axis(
    f[1, 1],
    title = "Proton trajectory",
    xlabel = "X",
    ylabel = "Y",
    limits = (x[1], x[end], y[1], y[end]),
    aspect = DataAspect()
)

lines!(ax, sol, idxs = (1, 2))

xgrid = [i for i in x, _ in y]
ygrid = [j for _ in x, j in y]
scatter!(ax, xgrid[:], ygrid[:], color = :tomato)

f = DisplayAs.PNG(f) #hide