# Phase-space (VDF) reconstruction utilities

module test_phasespace

using TestParticle
using TestParticle: get_gyroperiod, get_particle_crossings
using TestParticle: bin_velocity_space, project_vdf, refine_vdf_window, sample_velocity_ball
import TestParticle as TP
using Meshes
using BorisPushers
using StaticArrays
using Random
using LinearAlgebra: norm, ⋅
using Statistics: mean
using PlasmaPopulations: Oriented, BiMaxwellian, CartesianGridVDF, Population, pdf,
    thermal_speed
using Test

# ## Steady E×B drift of a bi-Maxwellian
#
# With `E·B = 0` there is no parallel acceleration, so a bi-Maxwellian centered on the
# E×B drift velocity `u_E = E×B/B²` is an exact steady solution of the Vlasov equation.
# The gyration only rotates the perpendicular velocity, preserving `|v_⊥ − u_E|` and
# `v_∥`, hence `f_det(v) = f_src(v)` and every reconstruction can be checked against the
# analytic source VDF.

const B_mag = 10.0e-9     # uniform magnetic field magnitude [T]
const B_vec = SA[0.0, 0.0, B_mag]  # along +z
const V_drift = -400.0e3  # desired E×B drift speed [m/s], along -x
const E_vec = SA[0.0, B_mag * V_drift, 0.0]  # [V/m]

get_E(r, t = 0.0) = E_vec
get_B(r, t = 0.0) = B_vec

const n0 = 3.0e6     # number density [m⁻³]
const T_par = 15.0   # parallel temperature [eV]
const T_perp = 45.0  # perpendicular temperature [eV]
const vdf = Oriented(
    BiMaxwellian(Proton; T_perp, T_para = T_par), B_vec; u = SA[V_drift, 0.0, 0.0]
)
## Source phase-space density [s³/m⁶]; `pdf` returns a normalized VDF.
f_src(v) = n0 * pdf(vdf, v)
const population = Population(Proton, vdf; n = n0)

const vth_perp = thermal_speed(Proton, T_perp)
const vradius = 3 * vth_perp     # velocity-space sampling radius [m/s]

const x_source = SA[300.0e3, 0.0, 0.0]  # source plane [m]
const x_detector = -200.0e3             # detector plane [m]
const tspan = (0.0, 4.0)                # [s]; > transport time 500 km / 400 km/s
const dt = get_gyroperiod(B_mag) / 40   # [s]
const param = prepare(; E = get_E, B = get_B, species = Proton)
@test TP._detector_point(Meshes.Point(1.0, 2.0, 3.0)) == SA[1.0, 2.0, 3.0]

const detector = Meshes.Plane(Meshes.Point(x_detector, 0.0, 0.0), Meshes.Vec(1.0, 0.0, 0.0))
const source_plane = Meshes.Plane(Meshes.Point(x_source...), Meshes.Vec(1.0, 0.0, 0.0))

const dv = 50.0e3                         # velocity bin width [m/s]
const v_edges = -1.0e6:dv:1.0e6           # [m/s]
const v_centers = (v_edges[1:(end - 1)] .+ v_edges[2:end]) ./ 2

"`‖rec − ref‖ / ‖ref‖` over the cells where `ref` exceeds `1e-4` of its maximum."
function relative_l2(rec, ref)
    m = ref .> maximum(ref) * 1.0e-4
    return norm(rec[m] .- ref[m]) / norm(ref[m])
end

const u0_dummy = SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

"Forward Monte Carlo: launch velocities drawn from the source VDF."
function prob_func_mc(prob, ctx)
    return remake(prob; u0 = SA[x_source..., rand(ctx.rng, vdf)...])
end

"Forward Liouville: launch velocities uniform in a ball around the drift velocity."
function prob_func_liouville(prob, ctx)
    v = sample_velocity_ball(ctx.rng, vradius; center = SA[V_drift, 0.0, 0.0])
    return remake(prob; u0 = SA[x_source..., v...])
end

"Projections of `f_src` by the midpoint rule on the bin centers, as for reconstructions."
function ana_projections()
    f = [f_src(SA[a, b, c]) for a in v_centers, b in v_centers, c in v_centers]
    return project_vdf(f, dv)
end

@testset "phase space utilities" begin

    @testset "sample_velocity_ball" begin
        R = 100.0e3
        center = SA[10.0e3, 0.0, -5.0e3]
        vs = [sample_velocity_ball(Xoshiro(42 + i), R; center) for i in 1:20000]
        @test all(norm(v .- center) ≤ R for v in vs)
        ## Uniform inside the ball: the mean radius is 3R/4 and the mean is the center.
        @test sum(norm(v .- center) for v in vs) / length(vs) ≈ 0.75R rtol = 0.02
        @test maximum(abs, sum(vs) ./ length(vs) .- center) < 0.05R
        @test sample_velocity_ball(Xoshiro(1), R) == sample_velocity_ball(Xoshiro(1), R)
    end

    @testset "bin_velocity_space" begin
        edges = -100.0e3:20.0e3:100.0e3
        ## The first sample sits at the center of bin (1, 1, 1); the second is outside.
        vs = [SA[-90.0e3, -90.0e3, -90.0e3], SA[1.0e6, 0.0, 0.0]]
        A = bin_velocity_space(vs, [2.0, 5.0], edges)
        @test A[1, 1, 1] == 2.0
        @test sum(A) == 2.0
        ## `average` takes the mean of the weights landing in a bin instead of their sum.
        C = bin_velocity_space([vs[1], vs[1]], [2.0, 4.0], edges; average = true)
        @test C[1, 1, 1] == 3.0
        @test C[1, 1, 2] == 0.0
        ## One range per axis
        D = bin_velocity_space(vs, [2.0, 5.0], (edges, -100.0e3:40.0e3:100.0e3, edges))
        @test size(D) == (10, 5, 10) && D[1, 1, 1] == 2.0
    end

    @testset "project_vdf" begin
        f_xy, f_xz, f_yz = project_vdf(ones(4, 5, 6), (1.0, 2.0, 3.0))
        @test size(f_xy) == (4, 5) && all(==(6 * 3.0), f_xy)
        @test size(f_xz) == (4, 6) && all(==(5 * 2.0), f_xz)
        @test size(f_yz) == (5, 6) && all(==(4 * 1.0), f_yz)
        d = CartesianGridVDF(0:3, 0:2:8, 0:3:15, ones(4, 5, 6))
        @test project_vdf(d) == (f_xy, f_xz, f_yz)
    end

    @testset "refine_vdf_window" begin
        vc = -300.0e3:100.0e3:300.0e3
        dv = 25.0e3
        v0 = first(vc) + dv / 2
        f = zeros(length(vc), length(vc), length(vc))
        f[2, 3, 4] = 1.0
        vx, vy, vz = refine_vdf_window(f, vc, vc, vc, v0, dv; margin = 50.0e3)
        @test step(vx) == dv && step(vy) == dv && step(vz) == dv
        ## Bounds are snapped onto the full grid `v0 + k·dv`.
        @test rem(first(vx) - v0, dv) ≈ 0 atol = 1.0e-6
        @test rem(first(vy) - v0, dv) ≈ 0 atol = 1.0e-6
        @test rem(first(vz) - v0, dv) ≈ 0 atol = 1.0e-6
        ## The populated cell, expanded by the margin, is enclosed.
        @test first(vx) ≤ vc[2] - 50.0e3 && last(vx) ≥ vc[2] + 50.0e3
        @test first(vy) ≤ vc[3] - 50.0e3 && last(vy) ≥ vc[3] + 50.0e3
        @test first(vz) ≤ vc[4] - 50.0e3 && last(vz) ≥ vc[4] + 50.0e3
        ## Without a populated region the coarse grid is returned unchanged.
        @test refine_vdf_window(zeros(size(f)), vc, vc, vc, v0, dv) == (vc, vc, vc)
    end
end

@testset "phase space reconstruction" begin
    ana = ana_projections()

    @testset "forward Monte Carlo" begin
        n = 50000
        prob = TraceProblem(u0_dummy, tspan, param; prob_func = prob_func_mc)
        sols = solve(
            prob, Boris(), EnsembleThreads(); dt,
            trajectories = n, seed = 42
        )
        vxi = [s.u[1][4] for s in sols.u]
        vs, ws_init = get_particle_crossings(sols, detector, vxi)
        ## Each macro-particle stands for `n0 / N` of the source density spread over one
        ## bin volume, and |v_x,src| / |v_x,det| turns its crossing flux into a density.
        w = n0 / (n * dv^3) .* abs.(ws_init) ./ abs.(first.(vs))
        f3d = bin_velocity_space(vs, w, v_edges)
        rec = project_vdf(f3d, dv)
        for i in 1:3
            @test relative_l2(rec[i], ana[i]) < 0.15
        end
        @test sum(f3d) * dv^3 ≈ n0 rtol = 0.1
    end

    @testset "get_particle_crossings serially" begin
        ## Ensembles below the chunking threshold accumulate without threads, and
        ## the multi-weight form has to agree with single-weight calls on them.
        n_small = 120
        prob = TraceProblem(u0_dummy, tspan, param; prob_func = prob_func_liouville)
        sols_small = solve(
            prob, Boris(), EnsembleThreads(); dt,
            trajectories = n_small, seed = 42
        )
        ws0 = [n0 * pdf(vdf, s.u[1][SA[4, 5, 6]]) for s in sols_small.u]
        vxi = [s.u[1][4] for s in sols_small.u]

        vs_t, (ws_t, wvxi_t) = get_particle_crossings(sols_small, detector, (ws0, vxi))
        vs_s, ws_s = get_particle_crossings(sols_small, detector, ws0)
        _, wvxi_s = get_particle_crossings(sols_small, detector, vxi)
        @test !isempty(vs_t)
        @test vs_t == vs_s
        @test ws_t == ws_s
        @test wvxi_t == wvxi_s
        @test length(vs_t) == length(ws_t) == length(wvxi_t)

        ## A scalar weight is shared by every crossing.
        _, ws_unit = get_particle_crossings(sols_small, detector)
        @test ws_unit == fill(1.0, length(vs_t))
    end

    isoutside = (u, p, t) -> u[1] < x_detector - 1.0e5 || u[1] > x_source[1] + 1.0e5
    backward(; kw...) = vdf_backward_trace(
        param, x_detector, source_plane, f_src;
        merge(
            (; v_range = 300.0e3, dv = 50.0e3, dt, tspan = (0.0, -8.0), isoutside),
            (; alg = Boris(), kw...)
        )...
    )
    ## Backward tracing returns `f_src` itself on its grid.
    on_grid(d) = [f_src(SA[a, b, c]) for a in d.vx, b in d.vy, c in d.vz]

    @testset "vdf_backward_trace" begin
        d = backward(; adaptive = false)
        full = -275.0e3:50.0e3:275.0e3
        @test d.vx ≈ full && d.vy ≈ full && d.vz ≈ full
        @test relative_l2(d.f, on_grid(d)) < 0.05

        ## The coarse pass only locates the support: the populated window of the full grid
        ## is retraced, dropping the empty +x half as the source drifts towards -x.
        da = backward()
        @test step(da.vx) == step(da.vz) == 50.0e3
        @test da.vx ⊆ d.vx && da.vy ⊆ d.vy && da.vz ⊆ d.vz
        @test length(da.vx) < length(full)
        @test relative_l2(da.f, on_grid(da)) < 0.05

        ## `bounds` clamps the refinement to a prescribed instead of the probed box.
        yz = (-100.0e3, 100.0e3)
        dc = backward(; bounds = ((first(full), last(full)), yz, yz))
        @test first(dc.vy) ≥ yz[1] && last(dc.vy) ≤ yz[2]
        @test first(dc.vz) ≥ yz[1] && last(dc.vz) ≤ yz[2]

        @test vdf_backward_trace(
            param, x_detector, source_plane, population;
            v_range = 300.0e3, dv = 50.0e3, dt, tspan = (0.0, -8.0), isoutside,
            adaptive = false, alg = Boris()
        ).f == d.f
        @test_throws "no trajectory reached" backward(;
            tspan = (0.0, -0.01), adaptive = false
        )

        ## A grid covering the whole beam recovers its moments.
        dw = backward(; v_range = 800.0e3)
        @test dw.n ≈ n0 rtol = 0.02
        @test mean(dw) ≈ [V_drift, 0, 0] atol = 2.0e3
    end

    @testset "vdf_forward_trace" begin
        d = vdf_forward_trace(
            param, x_source[1], detector, population;
            nparticles = 50000, vradius, tspan, dt, dv,
            center = SA[V_drift, 0.0, 0.0], alg = Boris()
        )
        @test d.vx ≈ v_centers
        for i in 1:3
            @test relative_l2(project_vdf(d)[i], ana[i]) < 0.15
        end
        @test d.n ≈ n0 rtol = 0.1
        @test mean(d) ≈ [V_drift, 0, 0] atol = 10.0e3
    end
end

end # module test_phasespace
