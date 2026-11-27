# Phase-space density reconstruction from traced trajectories. SI units throughout:
# velocities [m/s], phase-space densities [s³/m⁶], 2-D projections [s²/m⁵].

"""
    sample_velocity_ball([rng], radius; center = SA[0.0, 0.0, 0.0]) -> SVector{3}

Velocity drawn uniformly from the ball of `radius` around `center`, so that each of `N`
draws stands for the velocity volume `4πR³ / (3N)`.
"""
function sample_velocity_ball(rng::AbstractRNG, radius::Real; center = SA[0.0, 0.0, 0.0])
    return center .+ radius .* cbrt(rand(rng)) .* sample_unit_sphere(rng)
end

function sample_velocity_ball(radius::Real; center = SA[0.0, 0.0, 0.0])
    return sample_velocity_ball(default_rng(), radius; center)
end

"""
    bin_velocity_space(vs, weights, edges; average = false) -> CartesianGridVDF

Phase-space density on the bin centers of `edges` from the velocities `vs` with
`weights`; `edges` is a range shared by the three axes or one range per axis. Samples
outside the bins are dropped; empty bins are zero. Throws an `ArgumentError` when no
weight lands in a bin.

  - `average = false` sums the weights: for weights `f ΔV / Δv³` the result is the
    histogram estimate of the bin-averaged `f`, with counting noise `∝ 1/√n` for `n`
    samples per bin.
  - `average = true` averages them: a ratio estimate of `f` without the counting noise,
    for weights `f`.

For crossings of a detector plane, `ΔV` is the source velocity volume mapped onto the
detector, which brings in the flux ratio described in [`vdf_forward_trace`](@ref).
"""
function bin_velocity_space(vs, weights, edges; average = false)
    es = _triple(edges)
    A = _histogram(vs, weights, es, average)
    all(iszero, A) && throw(ArgumentError("no weight lands in the velocity bins"))
    return CartesianGridVDF(map(_bin_centers, es)..., A)
end

function _histogram(vs, weights, es, average)
    nb = map(e -> length(e) - 1, es)
    A = zeros(float(eltype(weights)), nb)
    counts = average ? zeros(Int, nb) : nothing

    @inbounds for i in eachindex(vs, weights)
        v = vs[i]
        I = map((e, x) -> floor(Int, (x - first(e)) / step(e)) + 1, es, (v[1], v[2], v[3]))
        all(map((k, n) -> 1 ≤ k ≤ n, I, nb)) || continue
        A[I...] += weights[i]
        counts === nothing || (counts[I...] += 1)
    end

    if counts !== nothing
        @inbounds for i in eachindex(A)
            counts[i] > 0 && (A[i] /= counts[i])
        end
    end
    return A
end

"""
    project_vdf(d::CartesianGridVDF) -> (f_xy, f_xz, f_yz)

Projections `∫ f dv_k` of `d` on each pair of velocity axes, by the midpoint rule on its
uniform grids.
"""
function project_vdf(d::CartesianGridVDF)
    proj(k, g) = dropdims(sum(d.f; dims = k); dims = k) .* step(g)
    return proj(3, d.vz), proj(2, d.vy), proj(1, d.vx)
end

"""
    refine_vdf_window(f, vx, vy, vz, v0, dv; margin = 0.0, relthresh = 1.0e-5,
        bounds = nothing)

Refine a coarse velocity grid down to the spacing `dv`, keeping only the populated part of
`f`. This is the second pass of the adaptive backward Liouville method: a cheap coarse
pass locates the region of velocity space that the detector can actually see, and only
that region is retraced on the fine grid.

Returns `(vx_new, vy_new, vz_new)` in the same unit as the inputs: the smallest box
enclosing every grid point where `f` exceeds `relthresh` times its maximum, expanded by
`margin` and clamped to `bounds`. The box bounds are snapped onto `v0 + k·dv` so that the
refined grid shares bin centers with the full grid starting at `v0`. `v0`, `dv` and
`margin` may be scalars or 3-tuples for per-axis values.

`bounds` is a `(lo, hi)` pair, or one such pair per axis. It defaults to the extent of the
coarse grid, which is the range the coarse pass actually probed; pass the extent of the
full display grid instead when the refined result has to live on that grid.
"""
function refine_vdf_window(
        f, vx, vy, vz, v0, dv; margin = 0.0, relthresh = 1.0e-5, bounds = nothing
    )
    kept = findall(f .> maximum(f) * relthresh)
    isempty(kept) && return (vx, vy, vz)

    vs = (vx, vy, vz)
    v0s = _triple(v0)
    dvs = _triple(dv)
    margins = _triple(margin)
    axis_bounds = _bounds(bounds, vs)

    return ntuple(3) do d
        lo = minimum(p -> vs[d][p[d]], kept) - margins[d]
        hi = maximum(p -> vs[d][p[d]], kept) + margins[d]
        _refine_range(lo, hi, axis_bounds[d][1], axis_bounds[d][2], v0s[d], dvs[d])
    end
end

_bounds(::Nothing, vs) = map(g -> (first(g), last(g)), vs)
_bounds(b::Tuple{<:Tuple, <:Tuple, <:Tuple}, _) = b
_bounds(b::Tuple, _) = (b, b, b)

function _refine_range(lo, hi, bound_lo, bound_hi, v0, dv)
    snap(x, r) = v0 + r((x - v0) / dv) * dv

    return range(max(bound_lo, snap(lo, floor)), min(bound_hi, snap(hi, ceil)); step = dv)
end

# `f_src` at the first crossing of `source` by the trajectory traced backward from `x0`
# with each velocity of the grid, zero for those that never cross
function _run_backward_pass(
        vx, vy, vz, x0, param, tspan, dt, source, f_src, isoutside, alg, ensemblealg
    )
    grid = CartesianIndices((length(vx), length(vy), length(vz)))

    function prob_func(prob, ctx)
        ix, iy, iz = Tuple(grid[ctx.sim_id])
        return remake(prob; u0 = SA[x0[1], x0[2], x0[3], vx[ix], vy[iy], vz[iz]])
    end

    function output_func(sol, i)
        st = get_first_crossing(sol, source)
        return (any(isnan, st) ? 0.0 : f_src(st[SA[4, 5, 6]]), false)
    end

    prob = TraceProblem(SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], tspan, param; prob_func)
    ensemble_prob = EnsembleProblem(prob; prob_func, output_func, safetycopy = false)
    solve_kwargs = isoutside === nothing ? (;) : (; isoutside)
    sols = solve(
        ensemble_prob, alg, ensemblealg;
        dt = -abs(dt), trajectories = length(grid), solve_kwargs...
    )
    return reshape(sols.u, size(grid))
end

"""
    vdf_backward_trace(param, detector, source, f_src; dt, alg, kw...) -> CartesianGridVDF

Phase-space density at the point `detector` by backward Liouville tracing: each velocity
of the grid of bin centers of `-v_range:dv:v_range` is traced backward over
`tspan` to its first crossing of the surface `source`, where `f_src` gives `f`. Grid
points that never reach `source` get `f = 0`.

`detector` is an x coordinate, a 3-vector or a Meshes `Point`. `f_src` maps a velocity to
`f`, or is a PlasmaPopulations `Population`.

With `adaptive = true`, a first pass at spacing `dv_coarse` locates where `f` exceeds
`relthresh` times its maximum; only that box, widened by `margin` and clamped to `bounds`
(see [`refine_vdf_window`](@ref)), is retraced at `dv`.

# Keywords

  - `v_range = 1000.0e3` (or one per axis), `dv = 20.0e3`
  - `tspan = (0.0, -20.0)`
  - `adaptive = true`, `dv_coarse = 3dv`, `margin = 3dv`, `relthresh = 1.0e-5`,
    `bounds = nothing` (the full grid)
  - `isoutside = nothing`, `ensemblealg = EnsembleThreads()`: passed to `solve`
"""
function vdf_backward_trace(
        param, detector, source, f_src;
        v_range = 1000.0e3,
        dv::Real = 20.0e3, dt::Real, tspan = (0.0, -20.0),
        adaptive::Bool = true, dv_coarse::Real = 3 * dv, margin::Real = 3 * dv,
        relthresh::Real = 1.0e-5, bounds = nothing, isoutside = nothing,
        alg, ensemblealg = EnsembleThreads()
    )
    x0 = _detector_point(detector)
    f = _density(f_src)
    v0 = map(r -> -r + dv / 2, _triple(v_range))  # first bin center per axis
    pass(vx, vy, vz) = _run_backward_pass(
        vx, vy, vz, x0, param, tspan, dt, source, f, isoutside, alg, ensemblealg
    )

    if adaptive
        coarse = map(r -> range(-r, r; step = dv_coarse), _triple(v_range))
        f_coarse = pass(coarse...)
        all(iszero, f_coarse) && throw(_unreached("source"))

        bnds = bounds === nothing ? map(a -> (a, -a), v0) : bounds
        grids = refine_vdf_window(
            f_coarse, coarse..., v0, dv; margin, relthresh, bounds = bnds
        )
    else
        grids = map(a -> range(a, -a; step = dv), v0)
    end
    f3d = pass(grids...)
    all(iszero, f3d) && throw(_unreached("source"))
    return CartesianGridVDF(grids..., f3d)
end

"""
    vdf_forward_trace(param, source, detector, f_src; vradius, dt, alg, kw...)

Phase-space density at `detector` by forward Liouville tracking, as a `CartesianGridVDF`
on the bin centers of `-v_range:dv:v_range`.

`nparticles` velocities drawn uniformly from the ball of `vradius` around `center` are
launched from the point `source` (as `detector` in [`vdf_backward_trace`](@ref)) over
`tspan`. By Liouville's theorem a crossing at `v` carries `f = f_src(v₀)`; it deposits
`f ΔV / Δv³` into its bin, with `ΔV` the ball volume per particle mapped onto the
detector, `ΔV = 4πR³/(3N) · |v₀ₓ| / |vₓ|`, which holds for a steady state that varies only
along x. `f_src` is as in [`vdf_backward_trace`](@ref).

# Keywords

  - `nparticles = 50000`, `center = SA[0.0, 0.0, 0.0]`, `seed = 42`
  - `v_range = 1000.0e3` (or one per axis), `dv = 20.0e3`
  - `tspan = (0.0, 20.0)`
  - `isoutside = nothing`, `ensemblealg = EnsembleThreads()`: passed to `solve`
"""
function vdf_forward_trace(
        param, source, detector, f_src;
        nparticles::Int = 50000, vradius::Real, center = SA[0.0, 0.0, 0.0],
        v_range = 1000.0e3,
        dv::Real = 20.0e3, tspan = (0.0, 20.0), dt::Real,
        isoutside = nothing, alg, ensemblealg = EnsembleThreads(), seed = 42
    )
    x_source = _detector_point(source)
    v_center = SVector{3, Float64}(center)
    f = _density(f_src)

    function prob_func(prob, ctx)
        v = sample_velocity_ball(ctx.rng, vradius; center = v_center)
        return remake(prob; u0 = vcat(x_source, v))
    end

    prob = TraceProblem(SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], tspan, param; prob_func)
    solve_kwargs = isoutside === nothing ? (;) : (; isoutside)
    sols = solve(
        prob, alg, ensemblealg;
        dt, trajectories = nparticles, seed, solve_kwargs...
    )

    ws0 = [f(s.u[1][SA[4, 5, 6]]) for s in sols.u]
    vxi = [s.u[1][4] for s in sols.u]
    vs, (ws, ws_vxi) = get_particle_crossings(sols, detector, (ws0, vxi))

    ΔV = 4π * vradius^3 / (3 * nparticles)
    w = ws .* (ΔV / dv^3) .* abs.(ws_vxi) ./ abs.(first.(vs))
    edges = map(r -> range(-r, r; step = dv), _triple(v_range))
    f3d = _histogram(vs, w, edges, false)
    all(iszero, f3d) && throw(_unreached("detector"))
    return CartesianGridVDF(map(_bin_centers, edges)..., f3d)
end

_density(f) = f
_density(p::Population) = let n = p.n, g = pdf(p.vdf)
    v -> n * g(v)
end

_unreached(what) = ArgumentError(
    "no trajectory reached `$what` with f > 0; widen `tspan` or the velocity range"
)

_bin_centers(e) = range(first(e) + step(e) / 2; step = step(e), length = length(e) - 1)

_triple(x::Tuple) = x
_triple(x) = (x, x, x)

_detector_point(x::Real) = SA[Float64(x), 0.0, 0.0]
_detector_point(p) = SA[Float64(p[1]), Float64(p[2]), Float64(p[3])]
