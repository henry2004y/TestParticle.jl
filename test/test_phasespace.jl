# Phase-space (VDF) reconstruction utilities

module test_phasespace

using TestParticle
import TestParticle as TP
using Meshes
using StaticArrays
using Random
using LinearAlgebra: norm
using VelocityDistributionFunctions
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
const p_par = n0 * TP.qᵢ * T_par
const p_perp = n0 * TP.qᵢ * T_perp
const vdf = TP.BiMaxwellian(
    SA[0.0, 0.0, 1.0], SA[V_drift, 0.0, 0.0], p_par, p_perp, n0; m = TP.mᵢ
)
## Source phase-space density [s³/m⁶]; `pdf` returns a normalized VDF.
f_src(v) = n0 * pdf(vdf, v)

const vth_perp = sqrt(2 * p_perp / (n0 * TP.mᵢ))
const vradius = 3 * vth_perp     # velocity-space sampling radius [m/s]

const x_source = SA[300.0e3, 0.0, 0.0]  # source plane [m]
const x_detector = -200.0e3             # detector plane [m]
const tspan = (0.0, 4.0)                # [s]; > transport time 500 km / 400 km/s
const dt = get_gyroperiod(B_mag) / 40   # [s]
const param = prepare(get_E, get_B; species = Proton)
const detector = Meshes.Plane(Meshes.Point(x_detector, 0.0, 0.0), Meshes.Vec(1.0, 0.0, 0.0))
const source_plane = Meshes.Plane(Meshes.Point(x_source...), Meshes.Vec(1.0, 0.0, 0.0))

const dv_km = 50.0                      # velocity bin width [km/s]
const v_edges = -1000:dv_km:1000        # [km/s]
const v_centers = bin_centers(v_edges)  # [km/s]

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

"Analytic 2-D projections on the bin centers shared by the three methods."
ana_projections() = (
    analytic_projection(f_src, 3, v_centers, v_centers, v_centers, dv_km),
    analytic_projection(f_src, 2, v_centers, v_centers, v_centers, dv_km),
    analytic_projection(f_src, 1, v_centers, v_centers, v_centers, dv_km),
)

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

    @testset "bin_centers" begin
        @test bin_centers(0:10:100) == 5:10:95
        @test length(v_centers) == length(v_edges) - 1
    end

    @testset "bin_velocity_space" begin
        edges = -100.0:20.0:100.0  # [km/s]
        ## The first sample sits at the center of bin (1, 1, 1); the second is outside.
        vs = [SA[-90.0e3, -90.0e3, -90.0e3], SA[1.0e6, 0.0, 0.0]]  # [m/s]
        A = bin_velocity_space(vs, [2.0, 5.0], edges)
        @test A[1, 1, 1] == 2.0
        @test sum(A) == 2.0
        ## `vx_source` rescales by |v_x,src| / |v_x,det| = 2.
        B = bin_velocity_space(vs, [2.0, 5.0], edges; vx_source = [-180.0e3, 0.0])
        @test B[1, 1, 1] == 4.0
        ## `average` takes the mean of the weights landing in a bin instead of their sum.
        C = bin_velocity_space([vs[1], vs[1]], [2.0, 4.0], edges; average = true)
        @test C[1, 1, 1] == 3.0
        @test C[1, 1, 2] == 0.0
    end

    @testset "project_vdf" begin
        f_xy, f_xz, f_yz = project_vdf(ones(4, 5, 6), 2.0)
        @test size(f_xy) == (4, 5) && all(==(6 * 2.0), f_xy)
        @test size(f_xz) == (4, 6) && all(==(5 * 2.0), f_xz)
        @test size(f_yz) == (5, 6) && all(==(4 * 2.0), f_yz)
    end

    @testset "analytic_projection" begin
        vth = 100.0e3  # [m/s]
        f_gauss(v) = exp(-(v[1]^2 + v[2]^2 + v[3]^2) / (2 * vth^2))
        g = -400.0:50.0:400.0
        ## `k` selects the integrated axis; the remaining two follow as `i < j`.
        for (k, i, j) in ((3, 1, 2), (2, 1, 3), (1, 2, 3))
            M = analytic_projection(f_gauss, k, g, g, g, 50.0)
            ref = zeros(length(g), length(g))
            for (bi, a) in enumerate(g), (bj, b) in enumerate(g)
                s = 0.0
                for c in g
                    v1 = i == 1 ? a : (j == 1 ? b : c)
                    v2 = i == 2 ? a : (j == 2 ? b : c)
                    v3 = i == 3 ? a : (j == 3 ? b : c)
                    s += f_gauss(SA[v1, v2, v3] * 1.0e3)
                end
                ref[bi, bj] = s * 1.0e18 * 50.0
            end
            @test M ≈ ref
        end
        @test_throws ArgumentError analytic_projection(f_gauss, 4, g, g, g, 50.0)

        ## Integrating a projection of the source VDF over velocity returns `n0` [m⁻³]
        ## and its first moment returns `n0·V_x` [m⁻³·km/s]. The grid has to cover the
        ## bulk velocity on both sides, hence the wide symmetric range.
        g = -1000.0:25.0:1000.0
        M = analytic_projection(f_src, 3, g, g, g, 25.0)
        m = velocity_moments(g, M; dv = 25.0)
        @test m.n ≈ n0 rtol = 1.0e-6
        @test m.nV ≈ n0 * V_drift * 1.0e-3 rtol = 1.0e-6
        @test velocity_moments((g, g, M); dv = 25.0).n == m.n
    end

    @testset "relative_l2" begin
        A = [0.0 1.0; 1.0 2.0]
        @test relative_l2(A, A) == 0.0
        @test relative_l2(2A, A) ≈ 1.0
        ## Cells below `thresh` of the reference maximum are left out, including the
        ## empty cell where the two disagree.
        @test relative_l2([5.0 1.0; 1.0 2.0], A) == 0.0
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
        @test !isempty(vs)
        ## Each macro-particle stands for `n0 / N` of the source density spread over one
        ## bin volume, hence `n0 [km⁻³] / (N · dv³)`.
        w = fill(n0 * 1.0e9 / (n * dv_km^3), length(vs))
        f3d = bin_velocity_space(vs, w, v_edges; vx_source = ws_init)
        rec = project_vdf(f3d, dv_km)
        ## The reconstruction integrates the third axis with the same midpoint rule as
        ## `analytic_projection`, so both share the same quadrature error.
        for i in 1:3
            @test relative_l2(rec[i], ana[i]) < 0.15
        end
        @test velocity_moments(v_centers, rec[1]; dv = dv_km).n ≈ n0 rtol = 0.1
    end

    @testset "forward Liouville" begin
        n = 50000
        prob = TraceProblem(u0_dummy, tspan, param; prob_func = prob_func_liouville)
        sols = solve(
            prob, Boris(), EnsembleThreads(); dt,
            trajectories = n, seed = 42
        )
        ws0 = [n0 * pdf(vdf, s.u[1][SA[4, 5, 6]]) for s in sols.u]
        vxi = [s.u[1][4] for s in sols.u]
        vs, ws = get_particle_crossings(sols, detector, ws0)
        _, ws_vxi = get_particle_crossings(sols, detector, vxi)
        vs_tup, (ws_tup, ws_vxi_tup) = get_particle_crossings(sols, detector, (ws0, vxi))
        @test vs_tup == vs
        @test ws_tup == ws
        @test ws_vxi_tup == ws_vxi
        ## Each sample stands for the source velocity volume `V_ball / N`.
        V_ball = (4 / 3) * π * (vradius * 1.0e-3)^3  # [km³/s³]
        f3d = bin_velocity_space(
            vs_tup, ws_tup .* (V_ball * 1.0e18 / (n * dv_km^3)), v_edges; vx_source = ws_vxi_tup
        )
        rec = project_vdf(f3d, dv_km)
        for i in 1:3
            @test relative_l2(rec[i], ana[i]) < 0.15
        end
        @test velocity_moments(v_centers, rec[1]; dv = dv_km).n ≈ n0 rtol = 0.1
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

    @testset "backward Liouville" begin
        ## Tracing the detector velocity grid back to the source and evaluating the
        ## source VDF there returns `f` on the grid, which here equals `f_src` itself.
        v = -300.0e3:50.0e3:300.0e3
        dims = (length(v), length(v), length(v))
        prob = vdf_grid_problem(v, v, v, SA[x_detector, 0.0, 0.0], param, (0.0, -8.0))
        sols = solve(
            prob, Boris(), EnsembleThreads(); dt = -dt, trajectories = prod(dims),
            isoutside = (u, p, t) -> u[1] < x_detector - 1.0e5 ||
                u[1] > x_source[1] + 1.0e5
        )
        f3d = vdf_backward(sols, source_plane, f_src, dims)
        ref = [f_src(SA[a, b, c]) * 1.0e18 for a in v, b in v, c in v]
        @test relative_l2(f3d, ref) < 0.05
    end

    @testset "vdf_backward_trace adaptive refinement" begin
        dv_bw = 50.0e3                # [m/s]
        v0_bw = -300.0e3 + dv_bw / 2  # first bin center shared by every fine grid
        full = range(v0_bw, -v0_bw; step = dv_bw)
        isoutside_bw = (u, p, t) -> u[1] < x_detector - 1.0e5 ||
            u[1] > x_source[1] + 1.0e5

        f3d_bw, (vx_bw, vy_bw, vz_bw) = vdf_backward_trace(
            param, x_detector, source_plane, f_src;
            v_range = 300.0e3, dv = dv_bw, dt, tspan = (0.0, -8.0),
            isoutside = isoutside_bw
        )
        ## The coarse pass only locates the support: the second pass retraces that
        ## region on the requested spacing, snapped onto the full grid lattice.
        @test step(vx_bw) == dv_bw && step(vy_bw) == dv_bw && step(vz_bw) == dv_bw
        @test rem(first(vy_bw) - v0_bw, dv_bw) ≈ 0 atol = 1.0e-6
        @test rem(first(vz_bw) - v0_bw, dv_bw) ≈ 0 atol = 1.0e-6
        @test first(vx_bw) ≥ first(full) && last(vx_bw) ≤ last(full)
        ## The source drifts towards -x, hence the empty +x half is dropped.
        @test length(vx_bw) < length(full)
        ## Backward tracing returns `f_src` on whichever grid it is sampled.
        ref_bw = [f_src(SA[a, b, c]) * 1.0e18 for a in vx_bw, b in vy_bw, c in vz_bw]
        @test size(f3d_bw) == (length(vx_bw), length(vy_bw), length(vz_bw))
        @test relative_l2(f3d_bw, ref_bw) < 0.05

        ## `bounds` clamps the refinement to a prescribed instead of the probed box.
        yz_bounds = (-100.0e3, 100.0e3)
        _, (vx_cl, vy_cl, vz_cl) = vdf_backward_trace(
            param, x_detector, source_plane, f_src;
            v_range = 300.0e3, dv = dv_bw, dt, tspan = (0.0, -8.0),
            bounds = ((first(full), last(full)), yz_bounds, yz_bounds),
            isoutside = isoutside_bw
        )
        @test step(vz_cl) == dv_bw
        @test first(vy_cl) ≥ yz_bounds[1] && last(vy_cl) ≤ yz_bounds[2]
        @test first(vz_cl) ≥ yz_bounds[1] && last(vz_cl) ≤ yz_bounds[2]
        @test length(vz_cl) < length(full)
    end

    @testset "embed_vdf" begin
        sub_centers = -50.0:20.0:50.0
        full_centers = -100.0:20.0:100.0
        M = ones(length(sub_centers), length(sub_centers))
        full_M = embed_vdf(sub_centers, sub_centers, full_centers, M)
        @test size(full_M) == (length(full_centers), length(full_centers))
        @test sum(full_M) == sum(M)
    end

    @testset "vdf_backward_trace and vdf_forward_trace" begin
        f3d_bw, (vx_bw, vy_bw, vz_bw) = vdf_backward_trace(
            param, x_detector, source_plane, f_src;
            v_range = 300.0e3, dv = 50.0e3, dt, tspan = (0.0, -8.0),
            adaptive = false,
            isoutside = (u, p, t) -> u[1] < x_detector - 1.0e5 ||
                u[1] > x_source[1] + 1.0e5
        )
        ref = [f_src(SA[a, b, c]) * 1.0e18 for a in vx_bw, b in vy_bw, c in vz_bw]
        @test relative_l2(f3d_bw, ref) < 0.05

        f3d_fw, proj_fw = vdf_forward_trace(
            param, x_source[1], detector, vdf, n0;
            nparticles = 50000, vradius, tspan, dt, dv_km,
            center = SA[V_drift, 0.0, 0.0]
        )
        for i in 1:3
            @test relative_l2(proj_fw[i], ana[i]) < 0.15
        end
    end
end

end # module test_phasespace
