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
    return f isa AbstractInterpolation || f isa AbstractExtrapolation
end

"""
    adapt_field_to_gpu(field::Field, backend::KA.Backend)

Adapt interpolation fields to GPU memory using Adapt.jl.
Analytic functions are returned unchanged.
"""
function adapt_field_to_gpu(field::Field, backend::KA.Backend)
    f = field.field_function
    if is_interpolation_field(f)
        # Adapt interpolation object to GPU
        adapted_func = Adapt.adapt(backend, f)
        return Field{is_time_dependent(field), typeof(adapted_func)}(adapted_func)
    end
    # Analytic fields don't need adaptation
    return field
end

# Fallback for ZeroField
adapt_field_to_gpu(field::ZeroField, backend::KA.Backend) = field


@kernel function boris_push_kernel!(
        @Const(xv_in), xv_out, @Const(q2m), @Const(dt),
        @Const(Ex_arr), @Const(Ey_arr), @Const(Ez_arr),
        @Const(Bx_arr), @Const(By_arr), @Const(Bz_arr)
    )
    i = @index(Global)

    x = xv_in[1, i]
    y = xv_in[2, i]
    z = xv_in[3, i]
    vx = xv_in[4, i]
    vy = xv_in[5, i]
    vz = xv_in[6, i]

    Ex = Ex_arr[i]
    Ey = Ey_arr[i]
    Ez = Ez_arr[i]
    Bx = Bx_arr[i]
    By = By_arr[i]
    Bz = Bz_arr[i]

    qdt_2m = q2m * 0.5 * dt

    vx_minus = vx + qdt_2m * Ex
    vy_minus = vy + qdt_2m * Ey
    vz_minus = vz + qdt_2m * Ez

    tx = qdt_2m * Bx
    ty = qdt_2m * By
    tz = qdt_2m * Bz

    t_mag2 = tx * tx + ty * ty + tz * tz
    factor = 2 / (1 + t_mag2)
    sx = factor * tx
    sy = factor * ty
    sz = factor * tz

    vpx = vx_minus + (vy_minus * tz - vz_minus * ty)
    vpy = vy_minus + (vz_minus * tx - vx_minus * tz)
    vpz = vz_minus + (vx_minus * ty - vy_minus * tx)

    vx_plus = vx_minus + (vpy * sz - vpz * sy)
    vy_plus = vy_minus + (vpz * sx - vpx * sz)
    vz_plus = vz_minus + (vpx * sy - vpy * sx)

    vx_new = vx_plus + qdt_2m * Ex
    vy_new = vy_plus + qdt_2m * Ey
    vz_new = vz_plus + qdt_2m * Ez

    xv_out[1, i] = x + vx_new * dt
    xv_out[2, i] = y + vy_new * dt
    xv_out[3, i] = z + vz_new * dt
    xv_out[4, i] = vx_new
    xv_out[5, i] = vy_new
    xv_out[6, i] = vz_new
end

"""
GPU kernel to evaluate interpolated fields directly on GPU.
This eliminates CPU-GPU data transfers for numerical fields.
"""
@kernel function evaluate_interp_fields_kernel!(
        Ex_arr, Ey_arr, Ez_arr, Bx_arr, By_arr, Bz_arr,
        @Const(xv), E_interp, B_interp, t
    )
    i = @index(Global)

    x = xv[1, i]
    y = xv[2, i]
    z = xv[3, i]

    # Evaluate interpolation directly on GPU
    # For 3D interpolation, call with (x, y, z)
    E_val = E_interp(x, y, z)
    B_val = B_interp(x, y, z)

    Ex_arr[i] = E_val[1]
    Ey_arr[i] = E_val[2]
    Ez_arr[i] = E_val[3]
    Bx_arr[i] = B_val[1]
    By_arr[i] = B_val[2]
    Bz_arr[i] = B_val[3]
end

function evaluate_fields_on_particles!(
        Ex_cpu, Ey_cpu, Ez_cpu, Bx_cpu, By_cpu, Bz_cpu, xv, Efunc, Bfunc, t
    )
    n_particles = size(xv, 2)
    xv_cpu = Array(xv)

    for i in 1:n_particles
        r = SVector(xv_cpu[1, i], xv_cpu[2, i], xv_cpu[3, i])
        E_val = Efunc(r, t)
        B_val = Bfunc(r, t)

        Ex_cpu[i] = E_val[1]
        Ey_cpu[i] = E_val[2]
        Ez_cpu[i] = E_val[3]
        Bx_cpu[i] = B_val[1]
        By_cpu[i] = B_val[2]
        Bz_cpu[i] = B_val[3]
    end

    return
end

function solve(
        prob::TraceProblem, backend::KA.Backend;
        dt::AbstractFloat, trajectories::Int = 1, savestepinterval::Int = 1,
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true
    )
    (; tspan, p, u0) = prob
    q2m, m, Efunc, Bfunc, _ = p
    T = eltype(u0)

    # Check if fields are interpolation objects and adapt to GPU
    use_gpu_interp = is_interpolation_field(Efunc.field_function) ||
        is_interpolation_field(Bfunc.field_function)

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

    Ex_arr = KA.zeros(backend, T, n_particles)
    Ey_arr = KA.zeros(backend, T, n_particles)
    Ez_arr = KA.zeros(backend, T, n_particles)
    Bx_arr = KA.zeros(backend, T, n_particles)
    By_arr = KA.zeros(backend, T, n_particles)
    Bz_arr = KA.zeros(backend, T, n_particles)

    xv_init = zeros(T, 6, n_particles)
    for i in 1:n_particles
        new_prob = prob.prob_func(prob, i, false)
        u0_i = new_prob.u0
        xv_init[:, i] .= u0_i
    end
    copyto!(xv_current, xv_init)

    Ex_cpu = zeros(T, n_particles)
    Ey_cpu = zeros(T, n_particles)
    Ez_cpu = zeros(T, n_particles)
    Bx_cpu = zeros(T, n_particles)
    By_cpu = zeros(T, n_particles)
    Bz_cpu = zeros(T, n_particles)

    # Initial field evaluation
    if use_gpu_interp
        # Use GPU kernel for interpolation
        eval_kernel! = evaluate_interp_fields_kernel!(backend, 256)
        eval_kernel!(
            Ex_arr, Ey_arr, Ez_arr, Bx_arr, By_arr, Bz_arr,
            xv_current, Efunc_gpu.field_function, Bfunc_gpu.field_function, tspan[1];
            ndrange = n_particles
        )
        KA.synchronize(backend)
    else
        # CPU evaluation for analytic fields
        evaluate_fields_on_particles!(
            Ex_cpu, Ey_cpu, Ez_cpu, Bx_cpu, By_cpu, Bz_cpu,
            xv_current, Efunc, Bfunc, tspan[1]
        )
        copyto!(Ex_arr, Ex_cpu)
        copyto!(Ey_arr, Ey_cpu)
        copyto!(Ez_arr, Ez_cpu)
        copyto!(Bx_arr, Bx_cpu)
        copyto!(By_arr, By_cpu)
        copyto!(Bz_arr, Bz_cpu)
    end

    kernel! = boris_push_kernel!(backend, 256)

    @kernel function velocity_back_kernel!(
            xv_out, @Const(xv_in), @Const(q2m_val), @Const(dt_val),
            @Const(Ex), @Const(Ey), @Const(Ez), @Const(Bx), @Const(By), @Const(Bz)
        )
        i = @index(Global)
        vx = xv_in[4, i]
        vy = xv_in[5, i]
        vz = xv_in[6, i]

        qdt_2m = q2m_val * 0.5 * dt_val
        vx_minus = vx + qdt_2m * Ex[i]
        vy_minus = vy + qdt_2m * Ey[i]
        vz_minus = vz + qdt_2m * Ez[i]

        tx = qdt_2m * Bx[i]
        ty = qdt_2m * By[i]
        tz = qdt_2m * Bz[i]
        t_mag2 = tx * tx + ty * ty + tz * tz
        factor = 2 / (1 + t_mag2)
        sx = factor * tx
        sy = factor * ty
        sz = factor * tz

        vpx = vx_minus + (vy_minus * tz - vz_minus * ty)
        vpy = vy_minus + (vz_minus * tx - vx_minus * tz)
        vpz = vz_minus + (vx_minus * ty - vy_minus * tx)

        vx_plus = vx_minus + (vpy * sz - vpz * sy)
        vy_plus = vy_minus + (vpz * sx - vpx * sz)
        vz_plus = vz_minus + (vpx * sy - vpy * sx)

        xv_out[4, i] = vx_plus + qdt_2m * Ex[i]
        xv_out[5, i] = vy_plus + qdt_2m * Ey[i]
        xv_out[6, i] = vz_plus + qdt_2m * Ez[i]
    end

    vback_kernel! = velocity_back_kernel!(backend, 256)
    vback_kernel!(
        xv_current, xv_current, q2m, -0.5 * dt,
        Ex_arr, Ey_arr, Ez_arr, Bx_arr, By_arr, Bz_arr; ndrange = n_particles
    )
    KA.synchronize(backend)

    sols = Vector{
        typeof(build_solution(prob, :boris, [tspan[1]], [SVector{6, T}(u0)])),
    }(undef, trajectories)

    saved_data = [Vector{SVector{6, T}}(undef, nout) for _ in 1:trajectories]
    saved_times = [Vector{typeof(tspan[1] + dt)}(undef, nout) for _ in 1:trajectories]
    iout_counters = zeros(Int, trajectories)

    if save_start
        xv_cpu = Array(xv_current)
        for i in 1:n_particles
            iout_counters[i] += 1
            saved_data[i][iout_counters[i]] = SVector{6, T}(xv_cpu[:, i])
            saved_times[i][iout_counters[i]] = tspan[1]
        end
    end

    for it in 1:nt
        t = tspan[1] + (it - 0.5) * dt

        if use_gpu_interp
            # Use GPU kernel for interpolation - no CPU-GPU transfers needed
            eval_kernel! = evaluate_interp_fields_kernel!(backend, 256)
            eval_kernel!(
                Ex_arr, Ey_arr, Ez_arr, Bx_arr, By_arr, Bz_arr,
                xv_current, Efunc_gpu.field_function, Bfunc_gpu.field_function, t;
                ndrange = n_particles
            )
            KA.synchronize(backend)
        else
            # Fallback to CPU evaluation for analytic fields
            evaluate_fields_on_particles!(
                Ex_cpu, Ey_cpu, Ez_cpu, Bx_cpu, By_cpu, Bz_cpu,
                xv_current, Efunc, Bfunc, t
            )
            copyto!(Ex_arr, Ex_cpu)
            copyto!(Ey_arr, Ey_cpu)
            copyto!(Ez_arr, Ez_cpu)
            copyto!(Bx_arr, Bx_cpu)
            copyto!(By_arr, By_cpu)
            copyto!(Bz_arr, Bz_cpu)
        end

        kernel!(
            xv_current, xv_next, q2m, dt,
            Ex_arr, Ey_arr, Ez_arr, Bx_arr, By_arr, Bz_arr;
            ndrange = n_particles
        )
        KA.synchronize(backend)

        xv_current, xv_next = xv_next, xv_current

        if save_everystep && it % savestepinterval == 0
            xv_cpu = Array(xv_current)
            for i in 1:n_particles
                if iout_counters[i] < nout
                    iout_counters[i] += 1
                    saved_data[i][iout_counters[i]] = SVector{6, T}(xv_cpu[:, i])
                    saved_times[i][iout_counters[i]] = tspan[1] + it * dt
                end
            end
        end
    end

    if save_end
        xv_cpu = Array(xv_current)
        for i in 1:n_particles
            if iout_counters[i] < nout
                iout_counters[i] += 1
                saved_data[i][iout_counters[i]] = SVector{6, T}(xv_cpu[:, i])
                saved_times[i][iout_counters[i]] = tspan[2]
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
