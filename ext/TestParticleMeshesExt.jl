module TestParticleMeshesExt

using TestParticle
import TestParticle: makegrid, get_cell_centers, prepare, build_interpolator,
    get_particle_crossings, get_first_crossing, get_particle_flux, get_particle_fluxes
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

function prepare(
        grid::CartesianGrid, E, B, F = TestParticle.ZeroField();
        order = 1, bc = TestParticle.FillExtrap(NaN), kw...
    )
    return TestParticle._prepare(
        E, B, F, makegrid(grid)...;
        gridtype = TestParticle.CartesianGrid, order, bc, kw...
    )
end

function prepare(
        grid::RectilinearGrid, E, B, F = TestParticle.ZeroField();
        order = 1, bc = TestParticle.FillExtrap(NaN), kw...
    )
    return TestParticle._prepare(
        E, B, F, makegrid(grid)...;
        gridtype = TestParticle.RectilinearGrid, order, bc, kw...
    )
end

# Virtual detectors
function get_particle_crossings(sol, surface::Union{Disk, Plane, Sphere}, weight = 1.0)
    u = sol.u
    T = float(eltype(u[1]))
    velocities = SVector{3, T}[]
    weights = typeof(weight)[]

    u1 = u[1]
    p1 = Point(u1[1], u1[2], u1[3])
    s1 = _signed_distance(p1, surface)

    @inbounds for i in 1:(length(u) - 1)
        u2 = u[i + 1]
        p2 = Point(u2[1], u2[2], u2[3])
        s2 = _signed_distance(p2, surface)

        _check_intersection!(
            velocities, weights, s1, s2, surface, u1, u2, weight, p1, p2
        )
        s1 = s2
        p1 = p2
        u1 = u2
    end

    return velocities, weights
end

function get_first_crossing(sol, surface::Union{Disk, Plane, Sphere})
    u = sol.u
    T = float(eltype(u[1]))

    u1 = u[1]
    p1 = Point(u1[1], u1[2], u1[3])
    s1 = _signed_distance(p1, surface)

    if s1 == 0
        return SVector{6, T}(u1)
    end

    @inbounds for i in 1:(length(u) - 1)
        u2 = u[i + 1]
        p2 = Point(u2[1], u2[2], u2[3])
        s2 = _signed_distance(p2, surface)

        if s1 * s2 < 0 || (s1 != 0 && s2 == 0)
            f = s1 / (s1 - s2)
            pcross = p1 + f * (p2 - p1)
            if _is_valid_intersection(pcross, surface)
                return muladd.(f, u2 - u1, u1)
            end
        end
        s1 = s2
        p1 = p2
        u1 = u2
    end

    return fill(T(NaN), SVector{6, T})
end

function get_first_crossing(sols::EnsembleSolution, surface::Union{Disk, Plane, Sphere})
    return get_first_crossing(sols.u, surface)
end

function get_first_crossing(
        sols::Union{AbstractVector, Tuple}, surface::Union{Disk, Plane, Sphere}
    )
    nsols = length(sols)
    T = float(eltype(first(sols).u[1]))
    res = Vector{SVector{6, T}}(undef, nsols)
    Threads.@threads for i in 1:nsols
        res[i] = get_first_crossing(sols[i], surface)
    end
    return res
end

function get_particle_flux(sol, surface::Union{Disk, Sphere}, weight = 1.0)
    vs, ws = get_particle_crossings(sol, surface, weight)

    inv_area = inv(area(surface).val)
    if isempty(ws)
        T = float(eltype(sol.u[1]))
        return zero(eltype(ws)) * inv_area, zero(SVector{3, T}) * inv_area
    end

    number_flux_density = sum(ws) * inv_area
    velocity_flux_density = sum(vs .* ws) * inv_area

    return number_flux_density, velocity_flux_density
end

function get_particle_fluxes(
        sols::Union{AbstractVector, Tuple}, surface::Union{Disk, Sphere},
        weights::Number = 1.0
    )
    return get_particle_fluxes(sols, surface, Base.Iterators.repeated(weights))
end

function get_particle_fluxes(
        sols::EnsembleSolution, surface::Union{Disk, Sphere},
        weights::Number = 1.0
    )
    return get_particle_fluxes(sols.u, surface, weights)
end

function get_particle_fluxes(sols::EnsembleSolution, surface::Union{Disk, Sphere}, weights)
    return get_particle_fluxes(sols.u, surface, weights)
end

function get_particle_fluxes(
        sols::Union{AbstractVector, Tuple}, surface::Union{Disk, Sphere}, weights
    )
    T = float(eltype(first(sols).u[1]))
    W = eltype(weights)
    total_n_flux = zero(W)
    total_v_flux = zero(SVector{3, T})

    @inbounds for (sol, w) in zip(sols, weights)
        total_n_flux, total_v_flux = _get_particle_flux_single_sum!(
            total_n_flux, total_v_flux, sol, surface, w
        )
    end

    inv_area = inv(area(surface).val)
    return total_n_flux * inv_area, total_v_flux * inv_area
end

function _get_particle_flux_single_sum!(total_n_flux, total_v_flux, sol, surface, w)
    u = sol.u
    u1 = u[1]
    p1 = Point(u1[1], u1[2], u1[3])
    s1 = _signed_distance(p1, surface)
    T = float(eltype(u1))

    for i in 1:(length(u) - 1)
        u2 = u[i + 1]
        p2 = Point(u2[1], u2[2], u2[3])
        s2 = _signed_distance(p2, surface)

        if s1 * s2 < 0 || (s1 != 0 && s2 == 0)
            f = s1 / (s1 - s2)
            pcross = p1 + f * (p2 - p1)
            if _is_valid_intersection(pcross, surface)
                vcross = SVector{3, T}(
                    muladd(f, u2[4] - u1[4], u1[4]),
                    muladd(f, u2[5] - u1[5], u1[5]),
                    muladd(f, u2[6] - u1[6], u1[6])
                )
                total_n_flux += w
                total_v_flux += vcross * w
            end
        end
        s1 = s2
        p1 = p2
        u1 = u2
    end

    return total_n_flux, total_v_flux
end

function get_particle_fluxes(
        sols::Union{AbstractVector, Tuple}, surfaces::AbstractVector{T},
        weights::Number = 1.0
    ) where {T <: Union{Disk, Sphere}}
    return get_particle_fluxes(sols, surfaces, Base.Iterators.repeated(weights))
end

function get_particle_fluxes(
        sols::EnsembleSolution, surfaces::AbstractVector{<:Union{Disk, Sphere}},
        weights::Number = 1.0
    )
    return get_particle_fluxes(sols.u, surfaces, weights)
end

function get_particle_fluxes(
        sols::EnsembleSolution, surfaces::AbstractVector{<:Union{Disk, Sphere}}, weights
    )
    return get_particle_fluxes(sols.u, surfaces, weights)
end

function get_particle_fluxes(
        sols::Union{AbstractVector, Tuple}, surfaces::AbstractVector{D}, weights
    ) where {D <: Union{Disk, Sphere}}
    nsurfaces = length(surfaces)
    T = float(eltype(first(sols).u[1]))

    total_n_fluxes = zeros(eltype(weights), nsurfaces)
    total_v_fluxes = zeros(SVector{3, T}, nsurfaces)

    s1_buf = Vector{T}(undef, nsurfaces)

    for (sol, w) in zip(sols, weights)
        _get_particle_fluxes_single_sum!(total_n_fluxes, total_v_fluxes, s1_buf, sol, surfaces, w)
    end

    @inbounds for j in 1:nsurfaces
        inv_area = inv(area(surfaces[j]).val)
        total_n_fluxes[j] *= inv_area
        total_v_fluxes[j] *= inv_area
    end

    return total_n_fluxes, total_v_fluxes
end

function _get_particle_fluxes_single_sum!(
        total_n_fluxes::AbstractVector{W},
        total_v_fluxes::AbstractVector{SVector{3, T}},
        s1s::AbstractVector{T},
        sol::S,
        surfaces::AbstractVector{D},
        w::W
    ) where {T, W, S, D}
    u = sol.u
    u1 = u[1]
    p1 = Point(u1[1], u1[2], u1[3])
    @inbounds for j in eachindex(surfaces)
        s1s[j] = _signed_distance(p1, surfaces[j])
    end

    @inbounds for i in 1:(length(u) - 1)
        u2 = u[i + 1]
        p2 = Point(u2[1], u2[2], u2[3])

        for j in eachindex(surfaces)
            surface = surfaces[j]
            s2 = _signed_distance(p2, surface)

            if s1s[j] * s2 < 0 || (s1s[j] != 0 && s2 == 0)
                f = s1s[j] / (s1s[j] - s2)
                pcross = p1 + f * (p2 - p1)
                if _is_valid_intersection(pcross, surface)
                    vcross = SVector{3, T}(
                        muladd(f, u2[4] - u1[4], u1[4]),
                        muladd(f, u2[5] - u1[5], u1[5]),
                        muladd(f, u2[6] - u1[6], u1[6])
                    )
                    total_n_fluxes[j] += w
                    total_v_fluxes[j] += vcross * w
                end
            end
            s1s[j] = s2
        end
        p1 = p2
        u1 = u2
    end

    return
end

function get_particle_crossings(
        sols::Union{AbstractVector, Tuple}, surface::Union{Disk, Plane, Sphere},
        weights::Number = 1.0
    )
    return get_particle_crossings(sols, surface, Base.Iterators.repeated(weights))
end

function get_particle_crossings(
        sols::EnsembleSolution, surface::Union{Disk, Plane, Sphere},
        weights::Number = 1.0
    )
    return get_particle_crossings(sols.u, surface, weights)
end

function get_particle_crossings(sols::EnsembleSolution, surface::Union{Disk, Plane, Sphere}, weights)
    return get_particle_crossings(sols.u, surface, weights)
end

function get_particle_crossings(
        sols::EnsembleSolution, surface::Union{Disk, Plane, Sphere},
        weights::Tuple{Vararg{AbstractVector}}
    )
    return get_particle_crossings(sols.u, surface, weights)
end

function get_particle_crossings(
        sols::Union{AbstractVector, Tuple}, surface::Union{Disk, Plane, Sphere},
        weights::Tuple{Vararg{AbstractVector}}
    )
    nsols = length(sols)
    T = float(eltype(first(sols).u[1]))
    nweights = length(weights)
    nthreads = Threads.nthreads()

    if nthreads > 1 && nsols > 200
        chunks = index_chunks(1:nsols; n = nthreads)
        th_vels = [SVector{3, T}[] for _ in 1:length(chunks)]
        th_counts = [ntuple(j -> eltype(weights[j])[], nweights) for _ in 1:length(chunks)]

        Threads.@threads for (tid, irange) in collect(enumerate(chunks))
            for i in irange
                w_tuple = ntuple(j -> weights[j][i], nweights)
                get_particle_crossings_single!(
                    th_vels[tid], th_counts[tid], sols[i], surface, w_tuple
                )
            end
        end

        velocities = reduce(vcat, th_vels)
        counts = ntuple(
            j -> reduce(vcat, [th_counts[t][j] for t in 1:length(chunks)]),
            nweights
        )
        return velocities, counts
    else
        velocities = SVector{3, T}[]
        sizehint!(velocities, nsols)
        counts = ntuple(j -> eltype(weights[j])[], nweights)
        for c in counts
            sizehint!(c, nsols)
        end
        w_iter = zip(weights...)
        @inbounds for (sol, w) in zip(sols, w_iter)
            get_particle_crossings_single!(velocities, counts, sol, surface, w)
        end
        return velocities, counts
    end
end

function get_particle_crossings(
        sols::Union{AbstractVector, Tuple}, surface::Union{Disk, Plane, Sphere}, weights
    )
    nsols = length(sols)
    T = float(eltype(first(sols).u[1]))
    nthreads = Threads.nthreads()

    if nthreads > 1 && nsols > 200 && weights isa AbstractVector
        chunks = index_chunks(1:nsols; n = nthreads)
        th_vels = [SVector{3, T}[] for _ in 1:length(chunks)]
        th_counts = [eltype(weights)[] for _ in 1:length(chunks)]

        Threads.@threads for (tid, irange) in collect(enumerate(chunks))
            for i in irange
                get_particle_crossings_single!(
                    th_vels[tid], th_counts[tid], sols[i], surface, weights[i]
                )
            end
        end

        velocities = reduce(vcat, th_vels)
        counts = reduce(vcat, th_counts)
        return velocities, counts
    else
        velocities = SVector{3, T}[]
        sizehint!(velocities, nsols)
        counts = eltype(weights)[]
        sizehint!(counts, nsols)

        @inbounds for (sol, w) in zip(sols, weights)
            get_particle_crossings_single!(velocities, counts, sol, surface, w)
        end

        return velocities, counts
    end
end

function get_particle_crossings_single!(velocities, weights, sol, surface, w)
    u = sol.u
    u1 = u[1]
    p1 = Point(u1[1], u1[2], u1[3])
    s1 = _signed_distance(p1, surface)

    @inbounds for i in 1:(length(u) - 1)
        u2 = u[i + 1]
        p2 = Point(u2[1], u2[2], u2[3])
        s2 = _signed_distance(p2, surface)

        _check_intersection!(
            velocities, weights, s1, s2, surface, u1, u2, w, p1, p2
        )
        s1 = s2
        p1 = p2
        u1 = u2
    end
    return
end

function get_particle_crossings(
        sols::Union{AbstractVector, Tuple},
        surfaces::AbstractVector{D}, weights::Number = 1.0
    ) where {D <: Union{Disk, Plane, Sphere}}
    return get_particle_crossings(sols, surfaces, Base.Iterators.repeated(weights))
end

function get_particle_crossings(
        sols::EnsembleSolution, surfaces::AbstractVector{<:Union{Disk, Plane, Sphere}},
        weights::Number = 1.0
    )
    return get_particle_crossings(sols.u, surfaces, weights)
end

function get_particle_crossings(
        sols::EnsembleSolution, surfaces::AbstractVector{<:Union{Disk, Plane, Sphere}},
        weights
    )
    return get_particle_crossings(sols.u, surfaces, weights)
end

function get_particle_crossings(
        sols::Union{AbstractVector, Tuple}, surfaces::AbstractVector{D},
        weights
    ) where {D <: Union{Disk, Plane, Sphere}}
    nsurfaces = length(surfaces)
    T = float(eltype(first(sols).u[1]))
    results_v = [SVector{3, T}[] for _ in 1:nsurfaces]
    results_w = [eltype(weights)[] for _ in 1:nsurfaces]

    s1_buf = Vector{T}(undef, nsurfaces)

    @inbounds for (sol, w) in zip(sols, weights)
        _get_particle_crossings_single!(results_v, results_w, s1_buf, sol, surfaces, w)
    end

    return results_v, results_w
end

function _get_particle_crossings_single!(
        results_v::Vector{Vector{SVector{3, T}}},
        results_w,
        s1s::AbstractVector{T},
        sol::S,
        surfaces::AbstractVector{D},
        w
    ) where {T, S, D}
    u = sol.u
    u1 = u[1]
    p1 = Point(u1[1], u1[2], u1[3])
    @inbounds for j in eachindex(surfaces)
        s1s[j] = _signed_distance(p1, surfaces[j])
    end

    @inbounds for i in 1:(length(u) - 1)
        u2 = u[i + 1]
        p2 = Point(u2[1], u2[2], u2[3])

        for j in eachindex(surfaces)
            surface = surfaces[j]
            s2 = _signed_distance(p2, surface)
            _check_intersection!(
                results_v[j], results_w[j], s1s[j], s2, surface, u1, u2, w, p1, p2
            )
            s1s[j] = s2
        end
        p1 = p2
        u1 = u2
    end

    return
end

@inline function _push_weight!(weights::AbstractVector, w)
    push!(weights, w)
    return
end

@inline function _push_weight!(weights::Tuple, w::Tuple)
    for j in eachindex(weights)
        push!(weights[j], w[j])
    end
    return
end

@inline function _check_intersection!(
        velocities::AbstractVector{SVector{3, T}},
        weights,
        s1, s2, surface::D, u1, u2, weight, p1::Point, p2::Point
    ) where {T, D}
    if s1 * s2 < 0 || (s1 != 0 && s2 == 0)
        f = s1 / (s1 - s2)
        pcross = p1 + f * (p2 - p1)
        if _is_valid_intersection(pcross, surface)
            vcross = SVector{3, T}(
                muladd(f, u2[4] - u1[4], u1[4]),
                muladd(f, u2[5] - u1[5], u1[5]),
                muladd(f, u2[6] - u1[6], u1[6])
            )
            push!(velocities, vcross)
            _push_weight!(weights, weight)
        end
    end

    return
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
        param = prepare(mesh, E, B)
        get_cell_centers(mesh)

        rect_mesh = RectilinearGrid(x, y, z)
        makegrid(rect_mesh)
        get_cell_centers(rect_mesh)
        param_rect = prepare(rect_mesh, E, B)

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
        get_particle_fluxes([sol], det)
        get_particle_fluxes([sol], [det])
    end
end

end
