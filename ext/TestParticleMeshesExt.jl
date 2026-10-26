module TestParticleMeshesExt

using TestParticle
import TestParticle: makegrid, get_cell_centers, prepare, build_interpolator,
    get_particle_crossings, get_first_crossing, get_particle_flux
using Meshes: coords, spacing, paramdim, CartesianGrid, RectilinearGrid, StructuredGrid,
    Plane, Disk, Point, normal, Sphere, area, Vec
using StaticArrays: SVector
using SciMLBase: EnsembleSolution
using LinearAlgebra: norm, ⋅
using PrecompileTools: @setup_workload, @compile_workload
using ChunkSplitters: index_chunks

# Grid build_interpolator forwarding
TestParticle.build_interpolator(::Type{<:CartesianGrid}, args...; kwargs...) =
    TestParticle.build_interpolator(TestParticle.CartesianGrid, args...; kwargs...)

TestParticle.build_interpolator(::Type{<:RectilinearGrid}, args...; kwargs...) =
    TestParticle.build_interpolator(TestParticle.RectilinearGrid, args...; kwargs...)

TestParticle.build_interpolator(::Type{<:StructuredGrid}, args...; kwargs...) =
    TestParticle.build_interpolator(TestParticle.StructuredGrid, args...; kwargs...)

"""
Return uniform range from 2D/3D CartesianGrid.
"""
function makegrid(grid::CartesianGrid)
    gridmin = coords(minimum(grid))
    gridmax = coords(maximum(grid))
    Δx = spacing(grid)
    dim = paramdim(grid)

    gridx = range(gridmin.x.val, gridmax.x.val, step = Δx[1].val)
    gridy = range(gridmin.y.val, gridmax.y.val, step = Δx[2].val)
    if dim == 3
        gridz = range(gridmin.z.val, gridmax.z.val, step = Δx[3].val)
        return gridx, gridy, gridz
    elseif dim == 2
        return gridx, gridy
    elseif dim == 1
        return (gridx,)
    end
end

"""
Return ranges from 2D/3D Meshes.jl RectilinearGrid.
"""
function makegrid(grid::RectilinearGrid)
    if paramdim(grid) == 3
        return grid.xyz[1], grid.xyz[2], grid.xyz[3]
    elseif paramdim(grid) == 2
        return grid.xyz[1], grid.xyz[2]
    end
end

"""
Return cell center coordinates from 2D/3D CartesianGrid.
"""
function get_cell_centers(grid::CartesianGrid)
    grid_coords = makegrid(grid)
    return map(
        r -> range(first(r) + step(r) / 2, step = step(r), length = length(r) - 1),
        grid_coords
    )
end

"""
Return cell center coordinates from 2D/3D RectilinearGrid.
"""
function get_cell_centers(grid::RectilinearGrid)
    grid_coords = makegrid(grid)
    if grid_coords === nothing
        return nothing
    end
    return map(c -> (c[1:(end - 1)] .+ c[2:end]) ./ 2, grid_coords)
end

_gridtype(::CartesianGrid) = TestParticle.CartesianGrid
_gridtype(::RectilinearGrid) = TestParticle.RectilinearGrid

function prepare(
        grid::Union{CartesianGrid, RectilinearGrid}; E = nothing, B = nothing, F = nothing,
        order = 1, bc = TestParticle.FillExtrap(NaN), kw...
    )
    return TestParticle._prepare(
        E, B, F, makegrid(grid)...; gridtype = _gridtype(grid), order, bc, kw...
    )
end

# Compatibility form, removed with those in TestParticle's `deprecated.jl`.
prepare(grid::Union{CartesianGrid, RectilinearGrid}, E, B, F = nothing; kw...) =
    prepare(grid; E, B, F, kw...)

# Virtual detectors

const Surface = Union{Disk, Plane, Sphere}
const AreaSurface = Union{Disk, Sphere}
const Trajectories = Union{AbstractVector, Tuple, EnsembleSolution}

_trajectories(sols::EnsembleSolution) = sols.u
_trajectories(sols) = sols

_weight_at(w::Number, i) = w
_weight_at(w::Tuple, i) = map(wj -> wj[i], w)
_weight_at(w, i) = w[i]

_weight_store(w::Number) = typeof(w)[]
_weight_store(w::Tuple) = map(wj -> eltype(wj)[], w)
_weight_store(w) = eltype(w)[]

_push_weight!(ws::AbstractVector, w) = push!(ws, w)
_push_weight!(ws::Tuple, w::Tuple) = foreach(push!, ws, w)

_vcat_weights(parts::AbstractVector{<:AbstractVector}) = reduce(vcat, parts)
_vcat_weights(parts::AbstractVector{<:Tuple}) =
    ntuple(j -> reduce(vcat, [p[j] for p in parts]), length(first(parts)))

_velocity_type(sol) = SVector{3, float(eltype(first(sol.u)))}

_point(u) = Point(u[1], u[2], u[3])

"""
    _foreach_crossing(f, sol, surface)

Call `f(θ, u1, u2)` for every saved step `u1 → u2` of `sol` that crosses `surface`, the
crossing lying at the fraction `θ` of the step. Stop early when `f` returns `true`.
"""
@inline function _foreach_crossing(f::F, sol, surface) where {F}
    u = sol.u
    u1 = first(u)
    p1 = _point(u1)
    s1 = _signed_distance(p1, surface)
    @inbounds for i in 2:length(u)
        u2 = u[i]
        p2 = _point(u2)
        s2 = _signed_distance(p2, surface)
        if s1 * s2 < 0 || (s1 != 0 && s2 == 0)
            θ = s1 / (s1 - s2)
            if _is_valid_intersection(p1 + θ * (p2 - p1), surface)
                # The callback is not inlined without the call-site hint.
                (@inline f(θ, u1, u2)) === true && return
            end
        end
        u1, p1, s1 = u2, p2, s2
    end
    return
end

@inline _crossing_velocity(::Type{V}, θ, u1, u2) where {V} = V(
    muladd(θ, u2[4] - u1[4], u1[4]),
    muladd(θ, u2[5] - u1[5], u1[5]),
    muladd(θ, u2[6] - u1[6], u1[6])
)

function _push_crossings!(vs, ws, sol, surface, w)
    _foreach_crossing(sol, surface) do θ, u1, u2
        push!(vs, _crossing_velocity(eltype(vs), θ, u1, u2))
        _push_weight!(ws, w)
        return
    end
    return
end

function get_particle_crossings(sol, surface::Surface, weight = 1.0)
    vs = _velocity_type(sol)[]
    ws = _weight_store(weight)
    _push_crossings!(vs, ws, sol, surface, weight)
    return vs, ws
end

function get_particle_crossings(sols::Trajectories, surface::Surface, weights = 1.0)
    trajs = _trajectories(sols)
    V = _velocity_type(first(trajs))
    nchunks = length(trajs) > 200 ? Threads.nthreads() : 1
    chunks = collect(enumerate(index_chunks(eachindex(trajs); n = nchunks)))
    parts = [(V[], _weight_store(weights)) for _ in chunks]

    Threads.@threads for (k, irange) in chunks
        vs, ws = parts[k]
        for i in irange
            _push_crossings!(vs, ws, trajs[i], surface, _weight_at(weights, i))
        end
    end

    return reduce(vcat, first.(parts)), _vcat_weights(last.(parts))
end

function get_particle_crossings(sols, surfaces::AbstractVector{<:Surface}, weights = 1.0)
    crossings = map(surface -> get_particle_crossings(sols, surface, weights), surfaces)
    return first.(crossings), last.(crossings)
end

function get_first_crossing(sol, surface::Surface)
    T = float(eltype(first(sol.u)))
    u1 = first(sol.u)
    _signed_distance(_point(u1), surface) == 0 && return SVector{6, T}(u1)

    state = Ref(fill(T(NaN), SVector{6, T}))
    _foreach_crossing(sol, surface) do θ, u1, u2
        state[] = SVector{6, T}(muladd.(θ, u2 - u1, u1))
        return true
    end
    return state[]
end

function get_first_crossing(sols::Trajectories, surface::Surface)
    trajs = _trajectories(sols)
    T = float(eltype(first(trajs).u[1]))
    states = Vector{SVector{6, T}}(undef, length(trajs))
    Threads.@threads for i in eachindex(trajs)
        states[i] = get_first_crossing(trajs[i], surface)
    end
    return states
end

function _add_flux(flux, sol, surface, w)
    V = typeof(last(flux))
    acc = Ref(flux)
    _foreach_crossing(sol, surface) do θ, u1, u2
        n, v = acc[]
        acc[] = (n + w, v + _crossing_velocity(V, θ, u1, u2) * w)
        return
    end
    return acc[]
end

get_particle_flux(sol, surface::AreaSurface, weight = 1.0) =
    get_particle_flux((sol,), surface, weight)

function get_particle_flux(sols::Trajectories, surface::AreaSurface, weights = 1.0)
    trajs = _trajectories(sols)
    flux = (zero(eltype(weights)), zero(_velocity_type(first(trajs))))
    for i in eachindex(trajs)
        flux = _add_flux(flux, trajs[i], surface, _weight_at(weights, i))
    end
    inv_area = inv(area(surface).val)
    return first(flux) * inv_area, last(flux) * inv_area
end

function get_particle_flux(sols, surfaces::AbstractVector{<:AreaSurface}, weights = 1.0)
    fluxes = map(surface -> get_particle_flux(sols, surface, weights), surfaces)
    return first.(fluxes), last.(fluxes)
end

@inline function _signed_distance(p::Point, surface::Disk)
    n = normal(surface.plane)
    v = p - surface.plane.p
    return (v ⋅ n).val
end

@inline function _signed_distance(p::Point, surface::Plane)
    n = normal(surface)
    v = p - surface.p
    return (v ⋅ n).val
end

@inline function _signed_distance(p::Point, surface::Sphere)
    v = p - surface.center
    return (norm(v) - surface.radius).val
end

@inline function _is_valid_intersection(p::Point, surface::Disk)
    center = surface.plane.p
    return (p - center) ⋅ (p - center) <= surface.radius^2
end

@inline function _is_valid_intersection(p::Point, surface::Union{Plane, Sphere})
    return true
end

@setup_workload begin
    @compile_workload begin
        x = range(-10, 10, length = 4)
        y = range(-10, 10, length = 6)
        z = range(-10, 10, length = 8)
        B = fill(0.0, 3, length(x), length(y), length(z))
        E = fill(0.0, 3, length(x), length(y), length(z))
        B[3, :, :, :] .= 10.0e-9
        E[3, :, :, :] .= 5.0e-10

        mesh = CartesianGrid(
            (first(x), first(y), first(z)), (last(x), last(y), last(z));
            dims = (length(x) - 1, length(y) - 1, length(z) - 1)
        )
        param = prepare(mesh; E, B)
        get_cell_centers(mesh)

        rect_mesh = RectilinearGrid(x, y, z)
        makegrid(rect_mesh)
        get_cell_centers(rect_mesh)
        param_rect = prepare(rect_mesh; E, B)

        # Mock a simple solution for crossing tests
        t_array = collect(0.0:1.0:2.0)
        u_array = [SVector{6, Float64}(-1.0 + t, 0.0, 0.0, 1.0, 0.0, 0.0) for t in t_array]
        struct PrecompileMockSol
            t::Vector{Float64}
            u::Vector{SVector{6, Float64}}
        end
        # Make the mock struct callable
        (s::PrecompileMockSol)(t) = s.u[1]

        sol = PrecompileMockSol(t_array, u_array)

        det = Disk(Plane(Point(0.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0)), 1.0)
        get_particle_flux([sol], det)
        get_particle_flux([sol], [det])
    end
end

end
