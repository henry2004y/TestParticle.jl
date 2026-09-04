# Phase-space (VDF) reconstruction utilities.
#
# These are the building blocks shared by the three ways of reconstructing a distribution
# function from traced trajectories: forward Monte-Carlo binning, forward Liouville
# tracking and backward Liouville tracing. They deliberately depend on no histogram
# package, so that all methods return the same `(vi, vj, f)` tuples and can be compared
# bin by bin.
#
# Convention: velocities are accepted in [m/s] and returned in [km/s]. Phase-space
# densities are returned in [s³/km⁶] on 3-D grids and in [s²/km⁵] for 2-D projections.
# Source distributions are passed as callables taking a velocity in [m/s] and returning
# `f` in [s³/m⁶], e.g. `v -> n * pdf(vdf, v)`.

"""
    sample_velocity_ball([rng], radius; center = SA[0.0, 0.0, 0.0]) -> SVector{3}

Sample a velocity [m/s] uniformly from the ball of `radius` [m/s] centered on `center`
[m/s].

Uniformity is what makes the sample usable as a quadrature rule in velocity space: each
draw stands for the same volume `4πR³ / (3N)`, which forward Liouville tracking needs in
order to convert a sampled `f` into a contribution to a detector bin.
"""
function sample_velocity_ball(rng::AbstractRNG, radius::Real; center = SA[0.0, 0.0, 0.0])
    return center .+ radius .* cbrt(rand(rng)) .* sample_unit_sphere(rng)
end

function sample_velocity_ball(radius::Real; center = SA[0.0, 0.0, 0.0])
    return sample_velocity_ball(default_rng(), radius; center)
end

"""
    bin_centers(edges)

Centers of the 1-D bins delimited by `edges`, in the same unit as `edges`.
"""
bin_centers(edges) = (edges[1:(end - 1)] .+ edges[2:end]) ./ 2

"""
    bin_velocity_space(vs, weights, edges; vx_source = nothing, average = false)

Accumulate velocity-space samples into a 3-D grid.

`vs` are sample velocities in [m/s] as returned by [`get_particle_crossings`](@ref) and
`weights` the weight of each sample. `edges` are the bin edges in [km/s], shared by all
three velocity axes. Samples outside the binned range are discarded, and the returned
array carries the unit of `weights`. Empty bins are zero.

Two estimators are available, and they differ in how the samples that land in a bin are
combined:

  * `average = false` (default) *sums* the weights, which is the histogram estimator of
    the bin-averaged `f`: feeding `f·ΔV/dv³` with `f` in [s³/km⁶] returns `f` in
    [s³/km⁶], and feeding the constant `n/(N·dv³)` turns binned counts into a
    Monte-Carlo estimate of `f`. Its error carries the counting noise `∝ 1/√n` of the
    samples per bin.
  * `average = true` takes the *mean* of the weights, a ratio estimator of the local `f`
    that divides out the counting noise. It is the better choice when the sampling spreads
    only a handful of samples over each bin, e.g. forward Liouville tracking of a
    uniformly sampled velocity sphere onto a fine detector grid.

Either way [`project_vdf`](@ref) turns the result into 2-D projections in [s²/km⁵].

If `vx_source` is given, each sample is rescaled by `|vx_source| / |vx|`. A trajectory
maps a source velocity-space volume onto a detector volume that differs by exactly that
ratio, so the rescale turns a source-referred weight into the detector-referred one. This
is what makes a flux-weighted launch ensemble consistent with the crossing flux recorded
at the detector.
"""
function bin_velocity_space(
        vs, weights, edges::AbstractRange; vx_source = nothing, average = false
    )
    nb = length(edges) - 1
    lo, dv = first(edges), step(edges)
    A = zeros(float(eltype(weights)), nb, nb, nb)
    n = average ? zeros(Int, nb, nb, nb) : nothing

    @inbounds for i in eachindex(vs)
        v = vs[i]
        ix = floor(Int, (v[1] * 1.0e-3 - lo) / dv) + 1
        iy = floor(Int, (v[2] * 1.0e-3 - lo) / dv) + 1
        iz = floor(Int, (v[3] * 1.0e-3 - lo) / dv) + 1
        (1 ≤ ix ≤ nb && 1 ≤ iy ≤ nb && 1 ≤ iz ≤ nb) || continue
        w = weights[i]
        vx_source === nothing || (w *= abs(vx_source[i]) / abs(v[1]))
        A[ix, iy, iz] += w
        n === nothing || (n[ix, iy, iz] += 1)
    end

    if n !== nothing
        @inbounds for i in eachindex(A)
            n[i] > 0 && (A[i] /= n[i])
        end
    end

    return A
end

"""
    project_vdf(f3d, dv) -> (f_xy, f_xz, f_yz)

Integrate a 3-D phase-space density `f3d` in [s³/km⁶] over each velocity axis in turn,
returning the projections `f(vᵢ, vⱼ) = ∫ f dvₖ` in [s²/km⁵]. `dv` is the bin width in
[km/s]. The array axes are ordered `(v₁, v₂, v₃)`, so the projections come out as
`(v₁, v₂)`, `(v₁, v₃)` and `(v₂, v₃)`.
"""
function project_vdf(f3d::AbstractArray{T, 3}, dv::Real) where {T}
    f_xy = dropdims(sum(f3d; dims = 3); dims = 3) .* dv
    f_xz = dropdims(sum(f3d; dims = 2); dims = 2) .* dv
    f_yz = dropdims(sum(f3d; dims = 1); dims = 1) .* dv

    return f_xy, f_xz, f_yz
end

"""
    analytic_projection(f, k, vi, vj, vk, dvk) -> Matrix

Analytic 2-D projection `f(vᵢ, vⱼ) = ∫ f dvₖ` in [s²/km⁵] of a 3-D phase-space density,
evaluated on the velocity grids `vi` and `vj` in [km/s] while the `k`-th axis is
integrated over `vk` in [km/s] with the midpoint rule at spacing `dvk`. The two remaining
axes follow from `k` as `i < j`: `k = 3` projects onto Vx–Vy, `k = 2` onto Vx–Vz and
`k = 1` onto Vy–Vz.

`f` takes a velocity in [m/s] and returns `f` in [s³/m⁶].
"""
function analytic_projection(f::F, k::Integer, vi, vj, vk, dvk::Real) where {F}
    k in (1, 2, 3) || throw(ArgumentError("integration axis k must be 1, 2 or 3, got $k"))
    T = float(promote_type(eltype(vi), eltype(vj), eltype(vk)))
    M = zeros(T, length(vi), length(vj))
    i, j = k == 1 ? (2, 3) : k == 2 ? (1, 3) : (1, 2)
    v = MVector{3, T}(0.0, 0.0, 0.0)

    @inbounds for (bi, a) in enumerate(vi), (bj, b) in enumerate(vj)
        v[i] = a
        v[j] = b
        s = zero(T)
        for c in vk
            v[k] = c
            s += f(SVector(v[1], v[2], v[3]) * 1.0e3)
        end
        M[bi, bj] = s * 1.0e18 * dvk
    end

    return M
end

"""
    relative_l2(rec, ref; thresh = 1.0e-4)

Relative L2 difference `‖rec − ref‖ / ‖ref‖`, evaluated over the entries where `ref` is
above `thresh` times its own maximum and `rec` is finite. Restricting to the populated
cells keeps the metric from being dominated by the empty tails of a reconstruction.
"""
function relative_l2(rec, ref; thresh = 1.0e-4)
    m = (ref .> maximum(ref) * thresh) .& isfinite.(rec)

    return norm(rec[m] .- ref[m]) / norm(ref[m])
end

"""
    velocity_moments(v, A; dv)
    velocity_moments(proj::Tuple; dv)

Density `n` [m⁻³] and momentum density `n·V` [m⁻³·km/s] of a 2-D velocity distribution
`A` in [s²/km⁵], sampled on the first-axis bin centers `v` in [km/s] with bin width `dv`
in [km/s]. Since the momentum is taken along the first axis, pass the Vx–Vy or the Vx–Vz
projection to obtain `n·Vx`.
"""
function velocity_moments(v, A; dv::Real)
    Af = Float64.(A)
    c = float(dv)^2 * 1.0e-9 # [s²/km⁵] · [km/s]² → [m⁻³]

    return (n = sum(Af) * c, nV = sum(Af .* v) * c)
end

velocity_moments(proj::Tuple; dv::Real) = velocity_moments(proj[1], proj[3]; dv)

"""
    vdf_grid_problem(vx, vy, vz, x0, param, tspan) -> TraceProblem

Build a `TraceProblem` that traces a velocity-space grid: trajectory `i` starts at
position `x0` [m] with the grid velocity `(vx[ix], vy[iy], vz[iz])` [m/s], the `z` index
varying fastest. This is the setup of the backward Liouville method, where the grid sits
at the detector and `tspan` runs backward in time.
"""
function vdf_grid_problem(vx, vy, vz, x0, param, tspan)
    nx, ny, nz = length(vx), length(vy), length(vz)

    function prob_func(prob, ctx)
        iz = (ctx.sim_id - 1) % nz + 1
        iy = ((ctx.sim_id - 1) ÷ nz) % ny + 1
        ix = ((ctx.sim_id - 1) ÷ (nz * ny)) % nx + 1
        return remake(prob; u0 = SA[x0[1], x0[2], x0[3], vx[ix], vy[iy], vz[iz]])
    end

    return TraceProblem(SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], tspan, param; prob_func)
end

_solutions(sols::EnsembleSolution) = sols.u
_solutions(sols) = sols

"""
    vdf_backward(sols, source, f_src, dims) -> Array{T, 3}

Phase-space density [s³/km⁶] on a velocity grid of size `dims = (nx, ny, nz)`, obtained
by tracing that grid backward in time and evaluating the source distribution on the
traced-back states, as prescribed by Liouville's theorem.

Each trajectory of `sols` is followed back to the surface `source` with
[`get_first_crossing`](@ref) and `f_src` is evaluated on the velocity there; grid points
whose trajectory never reaches `source` are left at zero. `f_src` takes a velocity in
[m/s] and returns `f` in [s³/m⁶].
"""
function vdf_backward(sols, source, f_src, dims)
    nx, ny, nz = dims
    states = _solutions(sols)
    T = float(eltype(first(states).u[1]))
    f = zeros(T, nx, ny, nz)

    for (i, sol) in enumerate(states)
        st = get_first_crossing(sol, source)
        any(isnan, st) && continue
        iz = (i - 1) % nz + 1
        iy = ((i - 1) ÷ nz) % ny + 1
        ix = ((i - 1) ÷ (nz * ny)) % nx + 1
        f[ix, iy, iz] = f_src(st[SA[4, 5, 6]]) * 1.0e18
    end

    return f
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

_triple(x::Tuple) = x
_triple(x) = (x, x, x)

_bounds(::Nothing, vs) = map(g -> (first(g), last(g)), vs)
_bounds(b::Tuple{<:Tuple, <:Tuple, <:Tuple}, _) = b
_bounds(b::Tuple, _) = (b, b, b)

function _refine_range(lo, hi, bound_lo, bound_hi, v0, dv)
    snap(x, r) = v0 + r((x - v0) / dv) * dv

    return range(max(bound_lo, snap(lo, floor)), min(bound_hi, snap(hi, ceil)); step = dv)
end
