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
        sols = TP.solve(
            prob, Boris(), EnsembleThreads(); dt, savestepinterval = 1,
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
        sols = TP.solve(
            prob, Boris(), EnsembleThreads(); dt, savestepinterval = 1,
            trajectories = n, seed = 42
        )
        ws0 = [n0 * pdf(vdf, s.u[1][SA[4, 5, 6]]) for s in sols.u]
        vxi = [s.u[1][4] for s in sols.u]
        vs, ws = get_particle_crossings(sols, detector, ws0)
        _, ws_vxi = get_particle_crossings(sols, detector, vxi)
        ## Each sample stands for the source velocity volume `V_ball / N`.
        V_ball = (4 / 3) * π * (vradius * 1.0e-3)^3  # [km³/s³]
        f3d = bin_velocity_space(
            vs, ws .* (V_ball * 1.0e18 / (n * dv_km^3)), v_edges; vx_source = ws_vxi
        )
        rec = project_vdf(f3d, dv_km)
        for i in 1:3
            @test relative_l2(rec[i], ana[i]) < 0.15
        end
        @test velocity_moments(v_centers, rec[1]; dv = dv_km).n ≈ n0 rtol = 0.1
    end

    @testset "backward Liouville" begin
        ## Tracing the detector velocity grid back to the source and evaluating the
        ## source VDF there returns `f` on the grid, which here equals `f_src` itself.
        v = -300.0e3:50.0e3:300.0e3
        dims = (length(v), length(v), length(v))
        prob = vdf_grid_problem(v, v, v, SA[x_detector, 0.0, 0.0], param, (0.0, -8.0))
        sols = TP.solve(
            prob, Boris(), EnsembleThreads(); dt = -dt, trajectories = prod(dims),
            savestepinterval = 1,
            isoutside = (u, p, t) -> u[1] < x_detector - 1.0e5 ||
                u[1] > x_source[1] + 1.0e5
        )
        f3d = vdf_backward(sols, source_plane, f_src, dims)
        ref = [f_src(SA[a, b, c]) * 1.0e18 for a in v, b in v, c in v]
        @test relative_l2(f3d, ref) < 0.05
    end
end

end # module test_phasespace
