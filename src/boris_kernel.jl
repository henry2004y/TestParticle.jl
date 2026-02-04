# GPU Boris solver using KernelAbstractions.jl

using KernelAbstractions
const KA = KernelAbstractions
using Adapt
using Interpolations: AbstractInterpolation, AbstractExtrapolation

# Helper functions for GPU interpolation support

"""
    is_interpolation_field(f)

Check if a field function is an interpolation object that can be adapted to GPU.
"""
function is_interpolation_field(f)
    return f isa AbstractInterpolation || f isa AbstractExtrapolation ||
        f isa FieldInterpolator
end

"""
    adapt_field_to_gpu(field::Field, backend::KA.Backend)

Adapt interpolation fields to GPU memory using Adapt.jl.
Analytic functions are returned unchanged.
"""
function adapt_field_to_gpu(field::Field, backend::KA.Backend)
    if backend isa KA.CPU
        return field
    end
    f = field.field_function
    if is_interpolation_field(f)
        # Unwrap FieldInterpolator to get the inner interpolation object
        itp = f isa FieldInterpolator ? f.itp : f

        # Adapt interpolation object to GPU
        adapted_func = Adapt.adapt(backend, itp)

        # Re-wrap in FieldInterpolator to maintain calling convention f(r)
        adapted_wrapper = FieldInterpolator(adapted_func)
        return Field{is_time_dependent(field), typeof(adapted_wrapper)}(adapted_wrapper)
    end
    # Analytic fields don't need adaptation (assuming they are GPU compatible functions)
    return field
end

# Fallback for ZeroField
adapt_field_to_gpu(field::ZeroField, backend::KA.Backend) = field


@inline function get_particle(xv, i)
    return SVector(xv[1, i], xv[2, i], xv[3, i]), SVector(xv[4, i], xv[5, i], xv[6, i])
end

@inline function boris_push_node!(i, xv_in, xv_out, q2m, dt, Efunc, Bfunc, t)
    r_vec, v_vec = get_particle(xv_in, i)

    # Evaluate fields directly
    E_val = Efunc(r_vec, t)
    B_val = Bfunc(r_vec, t)

    qdt_2m = q2m * 0.5 * dt

    v_new = boris_velocity_update(v_vec, E_val, B_val, qdt_2m)
    # Use scalar indexing for GPU compilation
    xv_out[1, i] = r_vec[1] + v_new[1] * dt
    xv_out[2, i] = r_vec[2] + v_new[2] * dt
    xv_out[3, i] = r_vec[3] + v_new[3] * dt
    xv_out[4, i] = v_new[1]
    xv_out[5, i] = v_new[2]
    xv_out[6, i] = v_new[3]

    return
end

@kernel function boris_push_kernel!(
        @Const(xv_in), xv_out, @Const(q2m), @Const(dt),
        Efunc, Bfunc, @Const(t)
    )
    i = @index(Global)
    boris_push_node!(i, xv_in, xv_out, q2m, dt, Efunc, Bfunc, t)
end

@inline function velocity_back_node!(i, xv_out, xv_in, q2m_val, dt_val, Efunc, Bfunc, t)
    r_vec, v_vec = get_particle(xv_in, i)

    # Evaluate fields at current position
    E_val = Efunc(r_vec, t)
    B_val = Bfunc(r_vec, t)

    qdt_2m = q2m_val * 0.5 * dt_val

    v_new = boris_velocity_update(v_vec, E_val, B_val, qdt_2m)
    # Use scalar indexing for GPU compilation
    xv_out[4, i] = v_new[1]
    xv_out[5, i] = v_new[2]
    xv_out[6, i] = v_new[3]

    return
end

@kernel function velocity_back_kernel!(
        xv_out, @Const(xv_in), @Const(q2m_val), @Const(dt_val),
        Efunc, Bfunc, @Const(t)
    )
    i = @index(Global)
    velocity_back_node!(i, xv_out, xv_in, q2m_val, dt_val, Efunc, Bfunc, t)
end

@inline function boris_step!(
        backend::KA.Backend, xv_in, xv_out, q2m, dt, Efunc, Bfunc, t,
        n_particles, workgroup_size
    )
    kernel! = boris_push_kernel!(backend, workgroup_size)
    kernel!(xv_in, xv_out, q2m, dt, Efunc, Bfunc, t; ndrange = n_particles)
    KA.synchronize(backend)
    return
end

@inline function boris_step!(
        ::KA.CPU, xv_in, xv_out, q2m, dt, Efunc, Bfunc, t, n_particles,
        workgroup_size
    )
    @inbounds for i in 1:n_particles
        boris_push_node!(i, xv_in, xv_out, q2m, dt, Efunc, Bfunc, t)
    end
    return
end

@inline function velocity_back_step!(
        backend::KA.Backend, xv_in, xv_out, q2m, dt, Efunc, Bfunc,
        t, n_particles, workgroup_size
    )
    kernel! = velocity_back_kernel!(backend, workgroup_size)
    kernel!(xv_out, xv_in, q2m, dt, Efunc, Bfunc, t; ndrange = n_particles)
    KA.synchronize(backend)
    return
end

@inline function velocity_back_step!(
        ::KA.CPU, xv_in, xv_out, q2m, dt, Efunc, Bfunc, t,
        n_particles, workgroup_size
    )
    @inbounds for i in 1:n_particles
        velocity_back_node!(i, xv_out, xv_in, q2m, dt, Efunc, Bfunc, t)
    end
    return
end


function _leapfrog_to_output(xv, Efunc, Bfunc, t, qdt_2m_half)
    T = eltype(xv)
    # Extract position and velocity (v^{n-1/2})
    r_vec = SVector(xv[1], xv[2], xv[3])
    v_vec = SVector{3}(xv[4], xv[5], xv[6])

    # Evaluate fields at current position and time
    E_val = Efunc(r_vec, t)
    B_val = Bfunc(r_vec, t)

    # Correct velocity to v^n using half-step push
    v_n = boris_velocity_update(v_vec, E_val, B_val, qdt_2m_half)

    return vcat(r_vec, v_n)
end


@inbounds function solve(
        prob::TraceProblem, backend::KA.Backend;
        dt::AbstractFloat, trajectories::Int = 1, savestepinterval::Int = 1,
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        workgroup_size::Int = 256
    )
    (; tspan, p, u0) = prob
    q2m, _, Efunc, Bfunc, _ = p
    T = eltype(u0)

    # Adapt interpolation fields to GPU memory
    Efunc_gpu = adapt_field_to_gpu(Efunc, backend)
    Bfunc_gpu = adapt_field_to_gpu(Bfunc, backend)

    ttotal = tspan[2] - tspan[1]
    nt = round(Int, abs(ttotal / dt))

    nout = 0
    if save_start
        nout += 1
    end
    if save_everystep
        steps = nt ÷ savestepinterval
        last_is_step = (nt > 0) && (nt % savestepinterval == 0)
        nout += steps
        if !save_end && last_is_step
            nout -= 1
        end
        if save_end && !last_is_step
            nout += 1
        end
    elseif save_end
        nout += 1
    end

    n_particles = trajectories
    xv_current = KA.zeros(backend, T, 6, n_particles)
    xv_next = KA.zeros(backend, T, 6, n_particles)

    # Optimization for CPU backend: alias buffers to avoid allocations
    is_cpu_accessible = xv_current isa Array

    if is_cpu_accessible
        xv_init = xv_current
    else
        xv_init = zeros(T, 6, n_particles)
    end

    for i in 1:n_particles
        new_prob = prob.prob_func(prob, i, false)
        u0_i = new_prob.u0
        xv_init[:, i] .= u0_i
    end

    if !is_cpu_accessible
        copyto!(xv_current, xv_init)
    end

    # Buffer for particle positions on CPU (used for saving data)
    if is_cpu_accessible
        xv_cpu_buffer = xv_current
    else
        xv_cpu_buffer = zeros(T, 6, n_particles)
    end

    sols = Vector{
        typeof(build_solution(prob, :boris, [tspan[1]], [SVector{6, T}(u0)])),
    }(undef, trajectories)

    saved_data = [Vector{SVector{6, T}}(undef, nout) for _ in 1:trajectories]
    saved_times = [Vector{typeof(tspan[1] + dt)}(undef, nout) for _ in 1:trajectories]
    iout_counters = zeros(Int, trajectories)

    if save_start
        if !is_cpu_accessible
            copyto!(xv_cpu_buffer, xv_current)
        end
        # If is_cpu_accessible, xv_cpu_buffer aliases xv_current, so it's already up to date
        for i in 1:n_particles
            iout_counters[i] += 1
            saved_data[i][iout_counters[i]] = SVector{6, T}(xv_cpu_buffer[:, i])
            saved_times[i][iout_counters[i]] = tspan[1]
        end
    end

    # Initial backward half-step
    velocity_back_step!(
        backend, xv_current, xv_current, q2m, -0.5 * dt,
        Efunc_gpu, Bfunc_gpu, tspan[1], n_particles, workgroup_size
    )

    for it in 1:nt
        t = tspan[1] + (it - 0.5) * dt

        boris_step!(
            backend, xv_current, xv_next, q2m, dt,
            Efunc_gpu, Bfunc_gpu, t, n_particles, workgroup_size
        )

        xv_current, xv_next = xv_next, xv_current

        if save_everystep && it % savestepinterval == 0
            if !is_cpu_accessible
                copyto!(xv_cpu_buffer, xv_current)
            end

            t_current = tspan[1] + it * dt
            qdt_2m_half = q2m * 0.5 * (0.5 * dt)

            for i in 1:n_particles
                if iout_counters[i] < nout
                    iout_counters[i] += 1
                    saved_data[i][iout_counters[i]] = _leapfrog_to_output(
                        @view(xv_cpu_buffer[:, i]), Efunc, Bfunc, t_current, qdt_2m_half
                    )
                    saved_times[i][iout_counters[i]] = t_current
                end
            end
        end
    end

    if save_end
        if !is_cpu_accessible
            copyto!(xv_cpu_buffer, xv_current)
        end
        t_current = tspan[2]
        qdt_2m_half = q2m * 0.5 * (0.5 * dt)

        for i in 1:n_particles
            if iout_counters[i] < nout
                iout_counters[i] += 1
                saved_data[i][iout_counters[i]] = _leapfrog_to_output(
                    @view(xv_cpu_buffer[:, i]), Efunc, Bfunc, t_current, qdt_2m_half
                )
                saved_times[i][iout_counters[i]] = t_current
            end
        end
    end

    for i in 1:trajectories
        actual_len = iout_counters[i]
        if actual_len < nout
            resize!(saved_data[i], actual_len)
            resize!(saved_times[i], actual_len)
        end

        interp = LinearInterpolation(saved_times[i], saved_data[i])
        sols[i] = build_solution(
            prob, :boris, saved_times[i], saved_data[i];
            interp = interp, retcode = ReturnCode.Default, stats = nothing
        )
    end

    return sols
end
