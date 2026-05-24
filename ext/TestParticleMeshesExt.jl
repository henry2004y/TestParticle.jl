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
    t, u = sol.t, sol.u
    T = float(eltype(u[1]))
    velocities = SVector{3, T}[]
    weights = typeof(weight)[]

    u1 = u[1]
    p1 = Point(u1[1], u1[2], u1[3])
    s1 = _signed_distance(p1, surface)

    @inbounds for i in 1:(length(t) - 1)
        u2 = u[i + 1]
        p2 = Point(u2[1], u2[2], u2[3])
        s2 = _signed_distance(p2, surface)

        _check_intersection!(
            velocities, weights, s1, s2, surface, sol, t[i], t[i + 1], weight, p1, p2
        )
        s1 = s2
        p1 = p2
    end

    return velocities, weights
end

function get_first_crossing(sol, surface::Union{Disk, Plane, Sphere})
    t, u = sol.t, sol.u
    T = float(eltype(u[1]))

    u1 = u[1]
    p1 = Point(u1[1], u1[2], u1[3])
    s1 = _signed_distance(p1, surface)

    if s1 == 0
        return SVector{6, T}(u1)
    end

    @inbounds for i in 1:(length(t) - 1)
        u2 = u[i + 1]
        p2 = Point(u2[1], u2[2], u2[3])
        s2 = _signed_distance(p2, surface)

        if s1 * s2 < 0 || (s1 != 0 && s2 == 0)
            f = s1 / (s1 - s2)
            pcross = p1 + f * (p2 - p1)
            if _is_valid_intersection(pcross, surface)
                tcross = muladd(f, t[i + 1] - t[i], t[i])
                ucross = sol(tcross)
                return SVector{6, T}(ucross)
            end
        end
        s1 = s2
        p1 = p2
    end

    return fill(T(NaN), SVector{6, T})
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
    t, u = sol.t, sol.u
    p1 = Point(u[1][1], u[1][2], u[1][3])
    s1 = _signed_distance(p1, surface)

    for i in 1:(length(t) - 1)
        u2 = u[i + 1]
        p2 = Point(u2[1], u2[2], u2[3])
        s2 = _signed_distance(p2, surface)

        if s1 * s2 < 0 || (s1 != 0 && s2 == 0)
            f = s1 / (s1 - s2)
            pcross = p1 + f * (p2 - p1)
            if _is_valid_intersection(pcross, surface)
                tcross = muladd(f, t[i + 1] - t[i], t[i])
                ucross = sol(tcross)
                vcross = SVector(ucross[4], ucross[5], ucross[6])
                total_n_flux += w
                total_v_flux += vcross * w
            end
        end
        s1 = s2
        p1 = p2
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
    t, u = sol.t, sol.u
    u1 = u[1]
    p1 = Point(u1[1], u1[2], u1[3])
    @inbounds for j in eachindex(surfaces)
        s1s[j] = _signed_distance(p1, surfaces[j])
    end

    @inbounds for i in 1:(length(t) - 1)
        u2 = u[i + 1]
        p2 = Point(u2[1], u2[2], u2[3])
        tl, tr = t[i], t[i + 1]

        for j in eachindex(surfaces)
            surface = surfaces[j]
            s2 = _signed_distance(p2, surface)

            if s1s[j] * s2 < 0 || (s1s[j] != 0 && s2 == 0)
                f = s1s[j] / (s1s[j] - s2)
                pcross = p1 + f * (p2 - p1)
                if _is_valid_intersection(pcross, surface)
                    tcross = muladd(f, tr - tl, tl)
                    ucross = sol(tcross)
                    vcross = SVector(ucross[4], ucross[5], ucross[6])
                    total_n_fluxes[j] += w
                    total_v_fluxes[j] += vcross * w
                end
            end
            s1s[j] = s2
        end
        p1 = p2
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
        sols::Union{AbstractVector, Tuple}, surface::Union{Disk, Plane, Sphere}, weights
    )
    T = float(eltype(first(sols).u[1]))
    velocities = SVector{3, T}[]
    sizehint!(velocities, length(sols))
    counts = eltype(weights)[]
    sizehint!(counts, length(sols))

    @inbounds for (sol, w) in zip(sols, weights)
        get_particle_crossings_single!(velocities, counts, sol, surface, w)
    end

    return velocities, counts
end

function get_particle_crossings_single!(velocities, weights, sol, surface, w)
    t, u = sol.t, sol.u
    u1 = u[1]
    p1 = Point(u1[1], u1[2], u1[3])
    s1 = _signed_distance(p1, surface)

    @inbounds for i in 1:(length(t) - 1)
        u2 = u[i + 1]
        p2 = Point(u2[1], u2[2], u2[3])
        s2 = _signed_distance(p2, surface)

        _check_intersection!(
            velocities, weights, s1, s2, surface, sol, t[i], t[i + 1], w, p1, p2
        )
        s1 = s2
        p1 = p2
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
        results_w::Vector{Vector{W}},
        s1s::AbstractVector{T},
        sol::S,
        surfaces::AbstractVector{D},
        w::W
    ) where {T, W, S, D}
    t, u = sol.t, sol.u
    u1 = u[1]
    p1 = Point(u1[1], u1[2], u1[3])
    @inbounds for j in eachindex(surfaces)
        s1s[j] = _signed_distance(p1, surfaces[j])
    end

    @inbounds for i in 1:(length(t) - 1)
        u2 = u[i + 1]
        p2 = Point(u2[1], u2[2], u2[3])
        tl, tr = t[i], t[i + 1]

        for j in eachindex(surfaces)
            surface = surfaces[j]
            s2 = _signed_distance(p2, surface)
            _check_intersection!(
                results_v[j], results_w[j], s1s[j], s2, surface, sol, tl, tr, w, p1, p2
            )
            s1s[j] = s2
        end
        p1 = p2
    end

    return
end

@inline function _check_intersection!(
        velocities::AbstractVector{SVector{3, T}},
        weights::AbstractVector{W},
        s1, s2, surface::D, sol::S, tl, tr, weight::W, p1::Point, p2::Point
    ) where {T, D, S, W}
    if s1 * s2 < 0 || (s1 != 0 && s2 == 0)
        f = s1 / (s1 - s2)
        pcross = p1 + f * (p2 - p1)
        if _is_valid_intersection(pcross, surface)
            tcross = muladd(f, tr - tl, tl)
            ucross = sol(tcross)
            vcross = SVector(ucross[4], ucross[5], ucross[6])
            push!(velocities, vcross)
            push!(weights, weight)
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
