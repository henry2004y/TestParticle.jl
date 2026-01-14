# Hybrid tracing with switching logic

import SciMLBase
using DiffEqCallbacks: DiscreteCallback, terminate!

"""
    check_gc_validity(x, v, p, t; epsilon=0.1)

Check if the guiding center approximation is valid.
Condition: r_L / L_B < epsilon.
r_L = v_perp / (q B / m)
L_B = B / |∇B|
"""
function check_gc_validity(x, v, p, t; epsilon = 0.1)
    q2m, m, Efunc, Bfunc, _ = p

    B, Bmag, b, ∇B, JB = get_B_parameters(x, t, Bfunc)

    if Bmag < eps(eltype(x))
        return false
    end

    # Calculate gyroradius r_L
    v_par = (v ⋅ b) .* b
    v_perp = v - v_par
    v_perp_mag = norm(v_perp)

    # If velocity is purely parallel, r_L is 0, GC is valid.
    if v_perp_mag < eps(eltype(v))
        return true
    end

    Ω = abs(q2m * Bmag)
    r_L = v_perp_mag / Ω

    # Calculate scale length L_B
    # L_B can be defined as B / |∇B| or similar.
    # Using 1/|∇B/B| = B / |∇B|
    ∇B_mag = norm(∇B)

    if ∇B_mag < eps(eltype(x))
        return true # Homogeneous field, GC valid
    end

    L_B = Bmag / ∇B_mag

    return r_L / L_B < epsilon
end

# State converters

"""
    gc_to_full(gc_state, p_gc, t)

Convert Guiding Center state (4-element) to Full Particle state (6-element).
Uses a fixed gyrophase of 0.
"""
function gc_to_full(gc_state, p_gc, t; phase = nothing)
    R = gc_state[SA[1:3...]]
    u_par = gc_state[4]

    q, m, μ, Efunc, Bfunc = p_gc

    # Need local B field
    B = Bfunc(R, t)
    Bmag = norm(B)
    b = B / Bmag

    # Calculate v_gyro magnitude from magnetic moment mu = m v_gyro^2 / (2B)
    v_gyro_mag = sqrt(2 * Bmag * μ / m)

    # Construct a perpendicular vector for phase 0.
    # We need a vector perpendicular to b.
    # Find the component with the smallest magnitude
    _, min_idx = findmin(abs.(b))
    arbitrary = SVector{3}(ntuple(i -> i == min_idx ? 1.0 : 0.0, 3))

    perp1 = normalize(arbitrary × b)
    perp2 = b × perp1

    # Random or fixed gyrophase
    θ = isnothing(phase) ? 0.0 : phase
    v_gyro = v_gyro_mag * (cos(θ) * perp1 + sin(θ) * perp2)

    # ExB drift
    E = Efunc(R, t)
    v_E = (E × b) / Bmag

    # Total velocity
    v = u_par * b + v_gyro + v_E

    # Position: x = R + rho.
    # rho = (b × v_gyro) / Omega
    q2m = q / m
    Ω = q2m * Bmag
    rho_vec = (b × v_gyro) / Ω

    x = R + rho_vec

    return vcat(x, v)
end

"""
    full_to_gc(full_state, p_full, t)

Convert Full Particle state (6-element) to Guiding Center state (4-element).
"""
function full_to_gc(full_state, p_full, t)
    x = full_state[SA[1:3...]]
    v = full_state[SA[4:6...]]

    q2m, m, Efunc, Bfunc, Ffunc = p_full

    # Calculate GC position
    R = get_gc(full_state, p_full)

    # Calculate u_par and mu at guiding center
    B_gc = Bfunc(R, t)
    Bmag_gc = norm(B_gc)
    b_gc = B_gc / Bmag_gc

    E_gc = Efunc(R, t)

    # Calculate parallel velocity
    u_par = v ⋅ b_gc

    # Calculate ExB drift at GC
    v_E = (E_gc × b_gc) / Bmag_gc

    # Calculate gyro-velocity (perpendicular velocity in frame of drift)
    v_perp_vec = v - u_par * b_gc - v_E
    v_perp² = v_perp_vec ⋅ v_perp_vec

    # Calculate new magnetic moment to conserve energy (approximately)
    μ_new = m * v_perp² / (2 * Bmag_gc)

    return SVector{4}(R[1], R[2], R[3], u_par), μ_new
end

"""
    get_gc_state_6d(u, mode, p_gc, p_full, t)

Helper to convert state `u` (either 4D GC or 6D Full) to 6D Guiding Center state (Position + GC Velocity).
"""
function get_gc_state_6d(u, mode, p_gc, p_full, t)
    if mode == :GC
        # u is 4D GC state
        v_gc = get_gc_velocity(u, p_gc, t)
        return vcat(u[SA[1:3...]], v_gc)
    else # mode == :Full
        # u is 6D Full state
        # Convert to 4D GC state first
        gc_state_4d, _ = full_to_gc(u, p_full, t)
        v_gc = get_gc_velocity(gc_state_4d, p_gc, t)
        return vcat(gc_state_4d[SA[1:3...]], v_gc)
    end
end

"""
    process_segment!(times, states, sol, mode, p_gc, p_full)

Process solution segment and append converted GC states to result containers.
"""
function process_segment!(times, states, sol, mode, p_gc, p_full)
    for i in eachindex(sol.t)
        t_val = sol.t[i]
        u_val = sol.u[i]
        state_6d = get_gc_state_6d(u_val, mode, p_gc, p_full, t_val)
        push!(times, t_val)
        push!(states, state_6d)
    end
    return
end

# Hybrid solver implementation

"""
    solve_hybrid(prob::ODEProblem, alg; epsilon=0.1, dt, kwargs...)

Solve the particle tracing problem using a hybrid method that switches between
Guiding Center (GC) tracing and Full Particle tracing based on validity conditions.

# Arguments
- `prob`: An `ODEProblem` setup for Guiding Center tracing (initial state is GC).
- `alg`: The solver algorithm (e.g., `Tsit5()`, `Vern9()`).
- `epsilon`: Threshold for GC validity (r_L / L_B < epsilon).
- `dt`: Time step (required for fixed step or initial hint).
- `kwargs`: Additional arguments passed to `solve`.

# Returns
- A `SciMLBase.ODESolution` containing the stitched trajectory (in Guiding Center coordinates).
"""
function solve_hybrid(prob::ODEProblem, alg; epsilon = 0.1, dt = nothing, kwargs...)
    # 1. Unpack GC parameters
    p_gc = prob.p # (q, m, mu, E, B)
    q, m, mu, E, B = p_gc

    # 2. Construct Full parameters
    # Trace param: (q2m, m, E, B, F)
    # Assume F=0 since GC doesn't support it
    q2m = q / m
    F_zero = ZeroField()
    p_full = (q2m, m, E, B, F_zero)

    tspan = prob.tspan
    t_end = tspan[2]

    # Results containers
    # Infer types from problem
    u_type = eltype(prob.u0)
    t_type = eltype(prob.tspan)
    times = t_type[]
    states = SVector{6, u_type}[]

    current_t = tspan[1]
    current_state = prob.u0 # Initial GC state (4-element)
    mode = :GC # Start in GC mode

    # Pre-define checking functions

    # Condition: r_L / L_B > epsilon -> Switch to Full
    function condition_gc_to_full(u, t, integrator)
        # u is GC state. Need to estimate validity.
        # r_L = v_perp / Omega. v_perp = sqrt(2B mu / m)
        R = u[SA[1:3...]]
        B_val = B(R, t)
        B_mag = norm(B_val)

        if B_mag < eps(eltype(u))
            return true
        end

        v_perp = sqrt(2 * B_mag * mu / m)
        Ω = abs(q2m * B_mag)
        r_L = v_perp / Ω

        # L_B
        _, _, _, ∇B, _ = get_B_parameters(R, t, B)
        ∇B_mag = norm(∇B)
        if ∇B_mag < eps(eltype(u))
            return false
        end

        L_B = B_mag / ∇B_mag

        ratio = r_L / L_B
        return ratio > epsilon
    end

    # Condition: r_L / L_B < 0.8 * epsilon -> Switch to GC (Hysteresis)
    function condition_full_to_gc(u, t, integrator)
        # u is Full state.
        x = u[SA[1:3...]]
        v = u[SA[4:6...]]

        # Use the shared validity check
        is_valid_gc = check_gc_validity(x, v, p_full, t; epsilon = 0.8 * epsilon)
        return is_valid_gc # If true, we should switch to GC
    end

    # Loop until t_end
    while current_t < t_end

        remaining_tspan = (current_t, t_end)

        if mode == :GC
            # Setup GC problem
            prob_current = ODEProblem(trace_gc!, current_state, remaining_tspan, p_gc)
            cb = DiscreteCallback(condition_gc_to_full, terminate!)
        else # Full Mode
            # Setup Full problem
            prob_current = ODEProblem(trace!, current_state, remaining_tspan, p_full)
            cb = DiscreteCallback(condition_full_to_gc, terminate!)
        end

        sol = SciMLBase.solve(prob_current, alg; callback = cb, dt = dt, kwargs...)

        # Process results
        process_segment!(times, states, sol, mode, p_gc, p_full)

        # Update state for next step
        current_t = sol.t[end]
        final_u = sol.u[end]

        if sol.retcode == :Terminated
            # Switched!
            if mode == :GC
                mode = :Full
                current_state = gc_to_full(final_u, p_gc, current_t)
            else
                mode = :GC
                current_state, new_mu = full_to_gc(final_u, p_full, current_t)
                # Update p_gc with new mu
                p_gc = (q, m, new_mu, E, B)
            end
        else
            break # Finished
        end
    end

    return SciMLBase.build_solution(prob, alg, times, states)
end
