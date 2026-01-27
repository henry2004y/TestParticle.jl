using Test
using TestParticle
using StaticArrays
using LinearAlgebra
using OrdinaryDiffEq
using Random
import TestParticle as TP

@testset "Utility" begin
    @testset "Basic Utility" begin
        # From runtests.jl
        V, B = 440.0e3, 5.0e-9
        r = get_gyroradius(V, B)
        @test r ≈ 919206.1737113602
        @test get_gyroperiod(B) ≈ 13.126233465754511
    end

    @testset "Mean Magnitude" begin
        # Test for get_mean_magnitude
        # 1D
        B1 = [1.0, 2.0, 3.0]
        @test get_mean_magnitude(B1) ≈ 3.7416573867739413

        # 3D vector field
        nx, ny, nz = 2, 2, 2
        B = zeros(3, nx, ny, nz)
        B[1, :, :, :] .= 1.0
        B[2, :, :, :] .= 2.0
        B[3, :, :, :] .= 2.0
        # Magnitude of each vector is sqrt(1^2 + 2^2 + 2^2) = sqrt(9) = 3
        # RMS should be 3.0
        @test get_mean_magnitude(B) ≈ 3.0

        # Verify against manual calculation for a varying field
        B2 = zeros(3, 2)
        B2[:, 1] = [1.0, 0.0, 0.0]
        B2[:, 2] = [0.0, 3.0, 4.0] # mag 5
        # RMS = sqrt((1^2 + 5^2)/2) = sqrt(26/2) = sqrt(13) ≈ 3.60555
        @test get_mean_magnitude(B2) ≈ sqrt(13)
    end

    @testset "Interpolation" begin
        # From runtests.jl
        begin # scalar interpolation
            x = range(-10, 10, length = 4)
            y = range(-10, 10, length = 6)
            z = range(-10, 10, length = 8)
            n = [
                Float32(i + j + k)
                    for i in eachindex(x), j in eachindex(y), k in eachindex(z)
            ]
            nfunc11 = TP.getinterp_scalar(n, x, y, z)
            @test nfunc11(SA[9, 0, 0]) == 11.85
            nfunc12 = TP.getinterp_scalar(n, x, y, z, 1, 2)
            @test nfunc12(SA[20, 0, 0]) == 9.5
            nfunc13 = TP.getinterp_scalar(n, x, y, z, 1, 3)
            @test nfunc13(SA[20, 0, 0]) == 12.0
            nfunc21 = TP.getinterp_scalar(n, x, y, z, 2)
            @test nfunc21(SA[9, 0, 0]) == 11.898528302013874
            nfunc22 = TP.getinterp_scalar(n, x, y, z, 2, 2)
            @test nfunc22(SA[20, 0, 0]) == 9.166666686534882
            nfunc23 = TP.getinterp_scalar(n, x, y, z, 2, 3)
            @test nfunc23(SA[20, 0, 0]) == 12.14705765247345
            nfunc31 = TP.getinterp_scalar(n, x, y, z, 3)
            @test nfunc31(SA[9, 0, 0]) == 11.882999392215163
            nfunc32 = TP.getinterp_scalar(n, x, y, z, 3, 2)
            @test nfunc32(SA[20, 0, 0]) == 9.124999547351358
            nfunc33 = TP.getinterp_scalar(n, x, y, z, 3, 3)
            @test nfunc33(SA[20, 0, 0]) == 12.191176189381315
        end

        begin # spherical interpolation
            r = range(0, 10, length = 11)
            θ = range(0, π, length = 11)
            ϕ = range(0, 2π, length = 11)
            # Vector field
            B = fill(0.0, 3, length(r), length(θ), length(ϕ))
            B[1, :, :, :] .= 1.0
            Bfunc = TP.getinterp(TP.StructuredGrid, B, r, θ, ϕ)
            @test Bfunc(SA[1, 1, 1]) ≈ [0.57735, 0.57735, 0.57735] atol = 1.0e-5
            # Scalar field
            A = ones(length(r), length(θ), length(ϕ))
            Afunc = TP.getinterp_scalar(TP.StructuredGrid, A, r, θ, ϕ)
            @test Afunc(SA[1, 1, 1]) == 1.0
            @test Afunc(SA[0, 0, 0]) == 1.0
        end

        begin # non-uniform spherical interpolation
            r = logrange(1.0, 10.0, length = 11)
            θ = range(0, π, length = 11)
            ϕ = range(0, 2π, length = 11)
            # Vector field
            B = fill(0.0, 3, length(r), length(θ), length(ϕ))
            B[1, :, :, :] .= 1.0
            Bfunc = TP.getinterp(TP.StructuredGrid, B, r, θ, ϕ)
            @test Bfunc(SA[1, 1, 1]) ≈ [0.57735, 0.57735, 0.57735] atol = 1.0e-5
            # Scalar field
            A = ones(length(r), length(θ), length(ϕ))
            Afunc = TP.getinterp_scalar(TP.StructuredGrid, A, r, θ, ϕ)
            @test Afunc(SA[1, 1, 1]) == 1.0
            @test Afunc(SA[0, 0, 0]) == 1.0
        end
    end

    @testset "SVector Interpolation" begin
        # From test_svector_interp.jl
        # Cartesian Grid
        nx, ny, nz = 10, 10, 10
        gridx = range(0.0, 10.0, length = nx)
        gridy = range(0.0, 10.0, length = ny)
        gridz = range(0.0, 10.0, length = nz)

        # Create random field
        Random.seed!(42)
        B_array = rand(3, nx, ny, nz)
        B_svector = reinterpret(reshape, SVector{3, Float64}, B_array)

        # Interpolator using Array{Float64, 4} (Implicitly converted internally)
        itp_implicit = TP.get_interpolator(TP.CartesianGrid, B_array, gridx, gridy, gridz)

        # Interpolator using Array{SVector, 3} (Explicitly passed)
        itp_explicit = TP.getinterp(TP.CartesianGrid, B_svector, gridx, gridy, gridz)

        # Compare
        pt = SA[5.5, 5.5, 5.5]
        @test itp_implicit(pt) ≈ itp_explicit(pt)
        @test itp_explicit(pt) isa SVector{3, Float64}

        # Spherical Grid
        r = range(1.0, 10.0, length = 10)
        theta = range(0.0, π, length = 10)
        phi = range(0.0, 2π, length = 10)

        B_array_sph = rand(3, 10, 10, 10)
        B_svector_sph = reinterpret(reshape, SVector{3, Float64}, B_array_sph)

        itp_sph_implicit = TP.get_interpolator(TP.StructuredGrid, B_array_sph, r, theta, phi)
        itp_sph_explicit = TP.getinterp(TP.StructuredGrid, B_svector_sph, r, theta, phi)

        pt_cart = SA[5.0, 0.0, 0.0] # On x-axis

        @test itp_sph_implicit(pt_cart) ≈ itp_sph_explicit(pt_cart)
        @test itp_sph_explicit(pt_cart) isa SVector{3, Float64}

        # 2D case
        B_array_2d = rand(3, nx, ny)
        B_svector_2d = reinterpret(reshape, SVector{3, Float64}, B_array_2d)

        itp_2d_implicit = TP.getinterp(TP.CartesianGrid, B_array_2d, gridx, gridy)
        itp_2d_explicit = TP.getinterp(TP.CartesianGrid, B_svector_2d, gridx, gridy)

        pt_2d = SA[5.5, 5.5]
        @test itp_2d_implicit(pt_2d) ≈ itp_2d_explicit(pt_2d)
    end

    @testset "CartesianNonUniform" begin
        # From test_svector_interp.jl
        x = range(0.0, 10.0, length = 11)
        y = range(0.0, 10.0, length = 11)
        z = range(0.0, 10.0, length = 11)

        # Vector field
        B = fill(0.0, 3, length(x), length(y), length(z))
        B[1, :, :, :] .= 1.0 # Bx = 1.0

        Bfunc = TP.getinterp(TP.RectilinearGrid, B, x, y, z)
        @test Bfunc(SA[5.0, 5.0, 5.0]) ≈ [1.0, 0.0, 0.0]

        # Scalar field
        A = fill(2.0, length(x), length(y), length(z))
        Afunc = TP.getinterp_scalar(TP.RectilinearGrid, A, x, y, z)
        @test Afunc(SA[5.0, 5.0, 5.0]) ≈ 2.0

        # Check interpolation values
        # Linear gradient
        A_grad = [i + j + k for i in x, j in y, k in z]
        Afunc_grad = TP.getinterp_scalar(TP.RectilinearGrid, A_grad, x, y, z)

        # Center point: 5.0, 5.0, 5.0 -> should match exactly for linear interpolation
        # value = 5.0 + 5.0 + 5.0 = 15.0
        @test Afunc_grad(SA[5.0, 5.0, 5.0]) ≈ 15.0

        # Off-grid point: 5.5, 5.5, 5.5
        # Should be 5.5 + 5.5 + 5.5 = 16.5
        @test Afunc_grad(SA[5.5, 5.5, 5.5]) ≈ 16.5

        # SVector array support
        B_sv = [SA[1.0, 0.0, 0.0] for i in x, j in y, k in z]
        Bfunc_sv = TP.getinterp(TP.RectilinearGrid, B_sv, x, y, z)
        @test Bfunc_sv(SA[5.0, 5.0, 5.0]) ≈ SA[1.0, 0.0, 0.0]

        # Non-uniform grid check
        x_nu = [0.0, 1.0, 4.0, 9.0]
        y_nu = [0.0, 1.0, 4.0, 9.0]
        z_nu = [0.0, 1.0, 4.0, 9.0]
        A_nu = [sqrt(i) + sqrt(j) + sqrt(k) for i in x_nu, j in y_nu, k in z_nu]

        # Check if we can interpolate on this
        Afunc_nu = TP.getinterp_scalar(TP.RectilinearGrid, A_nu, x_nu, y_nu, z_nu)

        # Point (4.0, 4.0, 4.0) -> sqrt(4)+sqrt(4)+sqrt(4) = 2+2+2=6
        @test Afunc_nu(SA[4.0, 4.0, 4.0]) ≈ 6.0

        # Point (1.0, 4.0, 9.0) -> 1+2+3=6
        @test Afunc_nu(SA[1.0, 4.0, 9.0]) ≈ 6.0
    end

    @testset "Time-dependent field interpolation" begin
        # 1. Setup simple time-dependent field E(x, t) = x * t
        # We will use a dummy grid for spatial interpolation just to map x -> x
        # But effectively we want the loader to return a function E_t(x) = x * t_fixed.

        # 3 time steps: t = 0, 1, 2.
        times = [0.0, 1.0, 2.0]

        # Loader function: returns a callable E(x) -> SVector{3}
        # For testing, we verify that it is called lazily.
        load_counts = Dict{Int, Int}()

        function loader(i)
            if !haskey(load_counts, i)
                load_counts[i] = 0
            end
            load_counts[i] += 1

            t_val = times[i]
            # E = scalar * t * x
            return x -> SVector{3}(x[1] * t_val, x[2] * t_val, x[3] * t_val)
        end

        itp = LazyTimeInterpolator(times, loader)

        x = SVector(1.0, 2.0, 3.0)

        # Test t=0.0 (exact match)
        @test itp(x, 0.0) ≈ SVector(0.0, 0.0, 0.0)
        @test load_counts[1] == 1

        # Test t=1.0 (exact match)
        @test itp(x, 1.0) ≈ SVector(1.0, 2.0, 3.0)
        @test load_counts[2] == 1

        # Test t=0.5 (interpolation)
        # E(0.5) should be 0.5 * x
        @test itp(x, 0.5) ≈ SVector(0.5, 1.0, 1.5)
        # Should reuse loaded fields
        @test load_counts[3] == 1

        # Test t=1.5 (interpolation between 1 and 2)
        @test itp(x, 1.5) ≈ SVector(1.5, 3.0, 4.5)
        @test load_counts[3] == 1

        # Test extrapolation/clamping
        @test itp(x, 3.0) ≈ itp(x, 2.0) # Should clamp to end
        @test itp(x, -1.0) ≈ itp(x, 0.0) # Should clamp to start
    end

    @testset "Gyroradius Utility" begin
        # From test_gyroradius.jl
        @testset "Non-relativistic" begin
            # Uniform B field
            B0 = 1.0e-8
            B_func(x, t) = SA[0.0, 0.0, B0]
            E_func(x, t) = SA[0.0, 0.0, 0.0]

            # Particle parameters
            q = TestParticle.qᵢ
            m = TestParticle.mᵢ

            # Initial state: v perpendicular to B
            v_perp = 1.0e5
            x0 = SA[0.0, 0.0, 0.0]
            v0 = SA[v_perp, 0.0, 0.0]
            stateinit = [x0..., v0...]
            tspan = (0.0, 1.0)

            param = prepare(E_func, B_func; species = Ion(1, 1)) # Proton
            prob = ODEProblem(trace!, stateinit, tspan, param)
            sol = solve(prob, Tsit5())

            # Theoretical gyroradius
            # r = m * v_perp / (q * B)
            r_expected = m * v_perp / (q * B0)

            r_calc = get_gyroradius(sol, 0.5)
            @test r_calc ≈ r_expected rtol = 1.0e-5
        end

        @testset "Relativistic" begin
            # Uniform B field
            B0 = 0.01
            B_func(x, t) = SA[0.0, 0.0, B0]
            E_func(x, t) = SA[0.0, 0.0, 0.0]

            # Particle parameters
            q = TestParticle.qₑ
            m = TestParticle.mₑ
            c = TestParticle.c

            # Initial state: relativistic velocity
            # Let v = 0.5c
            # γ = 1/sqrt(1-0.5^2) = 1/sqrt(0.75) ≈ 1.1547
            # u = γv
            v_mag = 0.5 * c
            gamma = 1.0 / sqrt(1 - (v_mag / c)^2)
            u_perp = gamma * v_mag

            x0 = SA[0.0, 0.0, 0.0]
            u0 = SA[u_perp, 0.0, 0.0] # γv
            stateinit = [x0..., u0...]
            tspan = (0.0, 1.0e-9)

            param = prepare(E_func, B_func; species = Electron)
            prob = ODEProblem(trace_relativistic!, stateinit, tspan, param)
            sol = solve(prob, Vern6())

            # Theoretical relativistic gyroradius
            # r = p_perp / (q * B) = m * (γv)_perp / (q * B)
            r_expected = m * u_perp / (abs(q) * B0)

            r_calc = get_gyroradius(sol, 0.5e-9)
            @test r_calc ≈ r_expected rtol = 1.0e-3
        end

        @testset "Zero Field" begin
            # B = 0 -> r = Inf
            B_func(x, t) = SA[0.0, 0.0, 0.0]
            E_func(x, t) = SA[0.0, 0.0, 0.0]

            x0 = SA[0.0, 0.0, 0.0]
            v0 = SA[1.0, 0.0, 0.0]
            stateinit = [x0..., v0...]
            tspan = (0.0, 1.0)

            param = prepare(E_func, B_func; species = Ion(1, 1))
            prob = ODEProblem(trace!, stateinit, tspan, param)
            sol = solve(prob, Tsit5())

            @test get_gyroradius(sol, 0.5) == Inf
        end

        @testset "ExB Drift Correction" begin
            # B = 10 nT, E in y direction
            B0 = 10.0e-9
            E0 = 1.0e-4 # V/m
            B_func(x, t) = SA[0.0, 0.0, B0]
            E_func(x, t) = SA[0.0, E0, 0.0]

            # Drift velocity: v_E = (E x B) / B^2
            # E = (0, E0, 0), B = (0, 0, B0)
            # E x B = (E0*B0, 0, 0)
            # v_E = (E0/B0, 0, 0)
            v_drift = E0 / B0 # 1e-4 / 10e-9 = 10,000 m/s

            # Initial state: particle with exactly the drift velocity
            x0 = SA[0.0, 0.0, 0.0]
            v0 = SA[v_drift, 0.0, 0.0]
            stateinit = [x0..., v0...]
            tspan = (0.0, 1.0)

            param = prepare(E_func, B_func; species = Ion(1, 1))
            prob = ODEProblem(trace!, stateinit, tspan, param)
            sol = solve(prob, Tsit5())

            # With pure drift, gyroradius should be effectively 0
            r_calc = get_gyroradius(sol, 0.5)
            @test r_calc ≈ 0.0 atol = 1.0e-3
        end

        @testset "Vector Dispatch" begin
            q = TestParticle.qᵢ
            m = TestParticle.mᵢ
            B0 = 1.0e-8
            v_perp = 1.0e5

            # Case 1: Vector inputs for V and B, same result as scalar
            V = SA[v_perp, 0.0, 0.0]
            B = SA[0.0, 0.0, B0]
            r_scalar = get_gyroradius(v_perp, B0; q, m)
            r_vector = get_gyroradius(V, B; q, m)
            @test r_vector ≈ r_scalar

            # Case 2: Vector inputs with E field correction
            E0 = 1.0e-4
            E = SA[0.0, E0, 0.0]
            B_vec = SA[0.0, 0.0, 10.0e-9] # same as ExB drift correction test
            v_drift = E0 / 10.0e-9
            V_drift_vec = SA[v_drift, 0.0, 0.0]

            # Pure drift should give ~0 gyroradius
            r_drift = get_gyroradius(V_drift_vec, B_vec, E; q, m)
            @test r_drift ≈ 0.0 atol = 1.0e-3

            # Compare with explicit scalar calculation without E field (which would be wrong)
            # Without E, r would be m*v_drift / qB
            r_wrong = get_gyroradius(v_drift, 10.0e-9; q, m)
            @test r_wrong > 0.1 # Should be large
        end
    end

    @testset "Loop Field" begin
        # From test_loop.jl
        # Parameters
        R = [0.0, 0.0, 0.0]
        a = 1.0
        I = 1.0
        n = [0.0, 0.0, 1.0]

        # Test field at the center
        r_center = [0.0, 0.0, 0.0]
        B_center = TestParticle.getB_loop(r_center, R, a, I, n)
        # B = mu0 * I / (2 * a) in z direction
        B_expected = TestParticle.μ₀ * I / (2 * a) * n
        @test B_center ≈ B_expected atol = 1.0e-15

        # Test field on the axis
        z = 2.0
        r_axis = [0.0, 0.0, z]
        B_axis = TestParticle.getB_loop(r_axis, R, a, I, n)
        # B = mu0 * I * a^2 / (2 * (a^2 + z^2)^(3/2))
        B_mag = TestParticle.μ₀ * I * a^2 / (2 * (a^2 + z^2)^1.5)
        @test B_axis ≈ B_mag * n atol = 1.0e-15

        # Test arbitrary orientation
        n_rot = [1.0, 0.0, 0.0]
        B_center_rot = TestParticle.getB_loop(r_center, R, a, I, n_rot)
        @test B_center_rot ≈ TestParticle.μ₀ * I / (2 * a) * n_rot atol = 1.0e-15

        # Test off-axis field consistency (check B_z at z=0 plane)
        # For z=0, B_z = mu0*I/(2*pi*a) * K(k^2) ? No, let's use the formula.
        # At z=0, B_rho = 0.
        # We can check continuity or just that it runs without error and gives reasonable numbers.

        r_off = [0.5, 0.0, 0.0]
        B_off = TestParticle.getB_loop(r_off, R, a, I, n)
        @test B_off[1] ≈ 0.0 atol = 1.0e-15 # No radial component if z=0? Wait.
        # If z=0, B_rho component:
        # factor * (z / rho) * ... -> 0 because z=0. Correct.
        # So B should be purely in z direction.
        @test abs(B_off[3]) > 0.0

        # Test rotated loop with off-axis point
        # Rotate loop by 90 degrees around Y axis, normal becomes X.
        n_x = [1.0, 0.0, 0.0]
        # Point that was at [0, 0, 2] (z-axis) should now be at [2, 0, 0] (x-axis) relative to loop
        r_x = [2.0, 0.0, 0.0]
        B_x = TestParticle.getB_loop(r_x, R, a, I, n_x)

        # Magnitude should match the previous axis test
        @test norm(B_x) ≈ B_mag atol = 1.0e-15
        @test B_x ≈ B_mag * n_x atol = 1.0e-15

        # Test translated loop
        R_new = [10.0, 10.0, 10.0]
        r_new = R_new + [0.0, 0.0, 2.0] # 2.0 along z from center
        B_trans = TestParticle.getB_loop(r_new, R_new, a, I, n)
        @test B_trans ≈ B_mag * n atol = 1.0e-15
    end

    @testset "Curvature Radius" begin
        # From test_curvature.jl
        function B_circular(x, t)
            r = hypot(x[1], x[3])
            r == 0 && return SA[0.0, 0.0, 0.0]
            # Magnitude = r
            # Vector = (-z, 0, x). Norm = r.
            # B = (-z, 0, x)
            return SA[-x[3], 0.0, x[1]]
        end

        @testset "Zero Field" begin
            B_zero(x, t) = SA[0.0, 0.0, 0.0]
            x = SA[1.0, 0.0, 0.0]
            @test get_curvature_radius(x, 0.0, B_zero) == Inf
        end

        @testset "Uniform Field" begin
            B_uniform(x, t) = SA[0.0, 0.0, 1.0]
            x = SA[1.0, 0.0, 0.0]
            @test get_curvature_radius(x, 0.0, B_uniform) == Inf
        end

        @testset "Circular Field" begin
            # Field lines are circles in x-z plane: B = (-z, 0, x)
            # Radius of curvature at (R, 0, 0) is R.
            R = 5.0
            @test get_curvature_radius(SA[R, 0.0, 0.0], 0.0, B_circular) ≈ R
            @test get_curvature_radius(SA[0.0, 0.0, R], 0.0, B_circular) ≈ R

            theta = π / 4
            @test get_curvature_radius(
                SA[R * cos(theta), 0.0, R * sin(theta)], 0.0, B_circular
            ) ≈ R
        end

        @testset "Adiabaticity" begin
            # Adiabaticity ϵ = ρ / Rc.
            # With R=5, q=1, m=2, μ=0.5:
            # Rc = R = 5.0, |B| = R = 5.0
            # ρ = sqrt(2*μ*m/|B|) / q ≈ 0.632
            # ϵ ≈ 0.126

            R = 5.0
            r = SA[R, 0.0, 0.0]

            q = 1.0
            m = 2.0
            μ = 0.5
            B_circular_local(x, t) = B_circular(x, t) # Explicit capture?

            val = get_adiabaticity(r, B_circular, q, m, μ, 0.0)
            ρ_expected = sqrt(2 * μ * m / R) / q
            Rc_expected = R
            @test val ≈ ρ_expected / Rc_expected

            # Test dispatch
            ts = TestParticle.Species(m, q)
            @test get_adiabaticity(r, B_circular, μ, 0.0; species = ts) ≈ val

            # Test Default Proton
            # Proton: q=qᵢ, m=mᵢ
            # val_p = get_adiabaticity(r, B_circular, 1e-19, 0.0) # μ ~ 1e-19
            # just check it runs
            @test isfinite(get_adiabaticity(r, B_circular, 1.0e-19, 0.0))

            # Test Negative Charge
            # Should result in positive adiabaticity due to abs(q)
            q_neg = -1.0
            val_neg = get_adiabaticity(r, B_circular, q_neg, m, μ, 0.0)
            @test val_neg > 0
            @test val_neg ≈ val

            # Test Zero B Field
            # Should return Inf
            B_zero_test(x, t) = SA[0.0, 0.0, 0.0]
            @test get_adiabaticity(r, B_zero_test, q, m, μ, 0.0) == Inf
        end
    end
end
