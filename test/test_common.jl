# Shared fixtures for the test suite.
#
# Every test file runs inside its own module so that it can be included on its
# own, and reaches the definitions collected here through `using ..test_common`.
# Keeping a single definition per field is not only tidier: `prepare` specializes
# on the type of the field function, so a fresh closure at every call site makes
# Julia compile another full set of methods.

module test_common

using TestParticle
import TestParticle as TP
using StaticArrays
using Meshes: CartesianGrid
using LinearAlgebra: norm

export uniform_B, uniform_B_strong, uniform_Ex, zero_E, zero_field,
    unit_Bz, unit_Ey, constant_E, gradient_B, time_varying_B, curved_B, curved_B_t,
    bottle_B, bottle_B0, bottle_α, sheared_B,
    uniform_grid, cartesian_grid, max_rel_diff, prob_func_vy_by_id

## Analytic fields ############################################################

"Uniform 10 nT magnetic field along +z."
uniform_B(x) = SA[0.0, 0.0, 1.0e-8]

"Uniform 0.01 T magnetic field along +z; a gyroperiod short enough to resolve."
uniform_B_strong(x) = SA[0.0, 0.0, 0.01]

"Uniform 1.0e-9 V/m electric field along +x."
uniform_Ex(x) = SA[1.0e-9, 0.0, 0.0]

"""
Zero electric field as an ordinary function.

Distinct from [`zero_field`](@ref) on purpose: the two take different dispatch
paths inside `Field`.
"""
zero_E(x) = SA[0.0, 0.0, 0.0]

"Zero field object."
const zero_field = ZeroField()

"Unit magnetic field along +z, with the time-dependent `(x, t)` signature."
unit_Bz(x, t) = SA[0.0, 0.0, 1.0]

"Unit electric field along +y, with the time-dependent `(x, t)` signature."
unit_Ey(x, t) = SA[0.0, 1.0, 0.0]

"Constant 1.0e5 V/m electric field along +y."
constant_E(x, t) = SA[0.0, 1.0e5, 0.0]

"Magnetic field along +z that grows linearly in x: Bz = 0.01 (1 + x)."
gradient_B(x, t) = SA[0.0, 0.0, 0.01 * (1.0 + x[1])]

"""
    time_varying_B(x, t)

10 nT ahead of a front advancing at 100 m/s and behind x = -10 m, zero in
between.
"""
function time_varying_B(x, t)
    Bz = if x[1] > 100t
        1.0e-8
    elseif x[1] < -10.0
        1.0e-8
    else
        0.0
    end
    return SA[0.0, 0.0, Bz]
end

"""
    curved_B(x)
    curved_B_t(x, t)

Divergence-free field whose lines are circles around x = -3 m; B_θ = 1/r, so
∂B_θ/∂θ = 0 keeps ∇ ⋅ B = 0. The `_t` spelling carries the `(x, t)` signature
that the derivative helpers require; keeping the two apart keeps `curved_B`
time-independent as far as `Field` is concerned.
"""
function curved_B(x)
    θ = atan(x[3] / (x[1] + 3))
    r = sqrt((x[1] + 3)^2 + x[3]^2)
    return SA[-1.0e-6 * sin(θ) / r, 0.0, 1.0e-6 * cos(θ) / r]
end

curved_B_t(x, t) = curved_B(x)

const bottle_B0 = 1.0e-4  # [T]
const bottle_α = 1.0e-2   # [m⁻²]

"""
    bottle_B(x, t)

Magnetic bottle carrying both curvature and a gradient, strong enough to push a
hybrid solver back and forth between the guiding center and the full orbit.
"""
function bottle_B(x, t)
    Bz = bottle_B0 * (1 + bottle_α * x[3]^2)
    Bx = -bottle_B0 * bottle_α * x[1] * x[3]
    By = -bottle_B0 * bottle_α * x[2] * x[3]
    return SA[Bx, By, Bz]
end

const sheared_B0 = 0.01   # [T]
const sheared_k = 100.0   # [m⁻¹]

"""
    sheared_B(x, t)

Field rotating rapidly in x, which keeps the adiabaticity parameter large
everywhere and forces the non-adiabatic branch.
"""
sheared_B(x, t) = SA[
    sheared_B0 * cos(sheared_k * x[1]),
    sheared_B0 * sin(sheared_k * x[1]),
    0.0,
]

## Grids ######################################################################

"Node coordinates of a uniform grid over [-10, 10] in every direction."
uniform_grid(nx::Integer, ny::Integer, nz::Integer) = (
    range(-10, 10; length = nx),
    range(-10, 10; length = ny),
    range(-10, 10; length = nz),
)

"`CartesianGrid` spanning the nodes `x`, `y` and `z`."
function cartesian_grid(x, y, z)
    return CartesianGrid(
        (first(x), first(y), first(z)), (last(x), last(y), last(z));
        dims = (length(x) - 1, length(y) - 1, length(z) - 1)
    )
end

"`CartesianGrid` spanning the nodes `x` and `y`."
function cartesian_grid(x, y)
    return CartesianGrid(
        (first(x), first(y)), (last(x), last(y));
        dims = (length(x) - 1, length(y) - 1)
    )
end

## Ensemble fixtures #########################################################

"""
    prob_func_vy_by_id(prob, ctx)

Ensemble initial state that scales v_y by the trajectory index, so every
trajectory differs. This lives here rather than in a single test file because
`EnsembleDistributed` serializes the function through its module, which the
workers can only resolve if they load this file too.
"""
function prob_func_vy_by_id(prob, ctx)
    return remake(
        prob; u0 = SA[
            prob.u0[1], prob.u0[2], prob.u0[3],
            prob.u0[4], ctx.sim_id * 1.0e5, prob.u0[6],
        ]
    )
end

## Comparison helpers #########################################################

"""
    max_rel_diff(a, b)

The largest difference between two states, relative to the state itself.
"""
function max_rel_diff(a, b)
    return maximum(zip(a, b)) do (u, v)
        return norm(collect(u) - collect(v)) / max(norm(collect(v)), one(eltype(v)))
    end
end

end # module test_common
