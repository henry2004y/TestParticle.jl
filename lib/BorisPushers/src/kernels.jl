# GPU Boris kernels using KernelAbstractions.jl

"""
    adapt_field_to_gpu(field, backend::Backend)

Adapt interpolation fields or parameters to GPU memory using Adapt.jl.
Analytic functions and CPU backend return the input unchanged.
"""
adapt_field_to_gpu(field, ::CPU) = field
adapt_field_to_gpu(field, backend::Backend) = Adapt.adapt(backend, field)

"""
    adapt_params(p, backend::Backend)

Adapt parameter container `p` for device execution.
Default implementation adapts electric and magnetic fields through `get_EField`
and `get_BField`.
"""
adapt_params(p, ::CPU) = p
function adapt_params(p, backend::Backend)
    E = adapt_field_to_gpu(get_EField(p), backend)
    B = adapt_field_to_gpu(get_BField(p), backend)
    q2m = get_q2m(p)
    return (q2m, nothing, E, B)
end

@inline function boris_update_xv!(i, xv_in, xv_out, p, dt, t, alg)
    @inbounds begin
        r = SVector(xv_in[i, 1], xv_in[i, 2], xv_in[i, 3])
        v_half = SVector(xv_in[i, 4], xv_in[i, 5], xv_in[i, 6])
    end

    r_new, v_half_new = advance_boris(v_half, r, dt, t, p, alg)

    @inbounds begin
        xv_out[i, 1] = r_new[1]
        xv_out[i, 2] = r_new[2]
        xv_out[i, 3] = r_new[3]
        xv_out[i, 4] = v_half_new[1]
        xv_out[i, 5] = v_half_new[2]
        xv_out[i, 6] = v_half_new[3]
    end

    return
end

@inline function boris_retard_v!(i, xv_in, xv_out, p, dt, t, alg)
    @inbounds begin
        r = SVector(xv_in[i, 1], xv_in[i, 2], xv_in[i, 3])
        v = SVector(xv_in[i, 4], xv_in[i, 5], xv_in[i, 6])
    end

    v_half = update_velocity_half(v, r, dt, t, p, alg)

    @inbounds begin
        xv_out[i, 4] = v_half[1]
        xv_out[i, 5] = v_half[2]
        xv_out[i, 6] = v_half[3]
    end

    return
end

@kernel function boris_velocity_kernel!(
        xv_out, @Const(xv_in), p, @Const(dt), @Const(t), @Const(offset), @Const(alg)
    )
    idx = @index(Global)
    i = idx + offset
    boris_retard_v!(i, xv_in, xv_out, p, dt, t, alg)
end

@kernel function boris_update_kernel!(
        @Const(xv_in), xv_out, p, @Const(dt), @Const(t), @Const(offset), @Const(alg)
    )
    idx = @index(Global)
    i = idx + offset
    boris_update_xv!(i, xv_in, xv_out, p, dt, t, alg)
end

@inline function boris_step!(
        backend::Backend, xv_in, xv_out, p, dt, t, irange, workgroup_size, alg
    )
    offset = irange.start - 1
    n_particles = length(irange)
    kernel! = boris_update_kernel!(backend, workgroup_size)
    kernel!(xv_in, xv_out, p, dt, t, offset, alg; ndrange = n_particles)
    synchronize(backend)
    return
end

@inline function boris_step!(::CPU, xv_in, xv_out, p, dt, t, irange, workgroup_size, alg)
    @inbounds for i in irange
        boris_update_xv!(i, xv_in, xv_out, p, dt, t, alg)
    end
    return
end

@inline function boris_velocity_step!(
        backend::Backend, xv_in, xv_out, p, dt, t, irange, workgroup_size, alg
    )
    offset = irange.start - 1
    n_particles = length(irange)
    kernel! = boris_velocity_kernel!(backend, workgroup_size)
    kernel!(xv_out, xv_in, p, dt, t, offset, alg; ndrange = n_particles)
    synchronize(backend)
    return
end

@inline function boris_velocity_step!(
        ::CPU, xv_in, xv_out, p, dt, t, irange, workgroup_size, alg
    )
    @inbounds for i in irange
        boris_retard_v!(i, xv_in, xv_out, p, dt, t, alg)
    end
    return
end

@kernel function boris_trajectory_kernel!(
        xv, p, @Const(dt), @Const(t0), @Const(nt), @Const(offset), @Const(alg)
    )
    idx = @index(Global)
    i = idx + offset
    @inbounds begin
        r = SVector(xv[i, 1], xv[i, 2], xv[i, 3])
        v = SVector(xv[i, 4], xv[i, 5], xv[i, 6])
    end

    v_half = update_velocity_half(v, r, dt, t0, p, alg)

    for it in 1:nt
        t = t0 + (it - 1) * dt
        r, v_half = advance_boris(v_half, r, dt, t, p, alg)
    end

    t_end = t0 + nt * dt
    v_end = update_velocity_node(v_half, r, dt, t_end, p, alg)

    @inbounds begin
        xv[i, 1] = r[1]
        xv[i, 2] = r[2]
        xv[i, 3] = r[3]
        xv[i, 4] = v_end[1]
        xv[i, 5] = v_end[2]
        xv[i, 6] = v_end[3]
    end
end

@kernel function boris_saveat_kernel!(
        saved_data, @Const(xv), @Const(plan_times), p, @Const(dt),
        @Const(t0), @Const(nt), @Const(nsave), @Const(dir),
        @Const(save_start), @Const(save_end), @Const(offset), @Const(alg)
    )
    idx = @index(Global)
    i = idx + offset
    @inbounds begin
        r = SVector(xv[i, 1], xv[i, 2], xv[i, 3])
        v = SVector(xv[i, 4], xv[i, 5], xv[i, 6])
    end

    iout = 0
    if save_start
        iout += 1
        @inbounds begin
            saved_data[i, 1, iout] = r[1]
            saved_data[i, 2, iout] = r[2]
            saved_data[i, 3, iout] = r[3]
            saved_data[i, 4, iout] = v[1]
            saved_data[i, 5, iout] = v[2]
            saved_data[i, 6, iout] = v[3]
        end
    end

    v_half = update_velocity_half(v, r, dt, t0, p, alg)
    isave = 1
    r_prev = r
    v_half_prev = v_half
    t_prev = t0

    for it in 1:nt
        t_current = t0 + it * dt
        r, v_half = advance_boris(v_half_prev, r_prev, dt, t_prev, p, alg)

        if isave <= nsave && @inbounds saveat_reached(plan_times[isave], t_current, dir)
            v_prev = update_velocity_node(v_half_prev, r_prev, dt, t_prev, p, alg)
            v_cur = update_velocity_node(v_half, r, dt, t_current, p, alg)

            y_prev = SVector(
                r_prev[1], r_prev[2], r_prev[3],
                v_prev[1], v_prev[2], v_prev[3]
            )
            y_cur = SVector(
                r[1], r[2], r[3],
                v_cur[1], v_cur[2], v_cur[3]
            )

            while isave <= nsave &&
                    @inbounds saveat_reached(plan_times[isave], t_current, dir)
                t_target = @inbounds plan_times[isave]
                y_target = saveat_interpolate(t_prev, y_prev, t_current, y_cur, t_target)
                iout += 1
                @inbounds begin
                    saved_data[i, 1, iout] = y_target[1]
                    saved_data[i, 2, iout] = y_target[2]
                    saved_data[i, 3, iout] = y_target[3]
                    saved_data[i, 4, iout] = y_target[4]
                    saved_data[i, 5, iout] = y_target[5]
                    saved_data[i, 6, iout] = y_target[6]
                end
                isave += 1
            end
        end

        r_prev = r
        v_half_prev = v_half
        t_prev = t_current
    end

    if save_end
        t_end = t0 + nt * dt
        v_end = update_velocity_node(v_half, r, dt, t_end, p, alg)
        iout += 1
        @inbounds begin
            saved_data[i, 1, iout] = r[1]
            saved_data[i, 2, iout] = r[2]
            saved_data[i, 3, iout] = r[3]
            saved_data[i, 4, iout] = v_end[1]
            saved_data[i, 5, iout] = v_end[2]
            saved_data[i, 6, iout] = v_end[3]
        end
    end
end

@kernel function boris_everystep_kernel!(
        saved_data, @Const(xv), p, @Const(dt),
        @Const(t0), @Const(nt), @Const(save_start), @Const(save_end),
        @Const(offset), @Const(alg)
    )
    idx = @index(Global)
    i = idx + offset
    @inbounds begin
        r = SVector(xv[i, 1], xv[i, 2], xv[i, 3])
        v = SVector(xv[i, 4], xv[i, 5], xv[i, 6])
    end

    iout = 0
    if save_start
        iout += 1
        @inbounds begin
            saved_data[i, 1, iout] = r[1]
            saved_data[i, 2, iout] = r[2]
            saved_data[i, 3, iout] = r[3]
            saved_data[i, 4, iout] = v[1]
            saved_data[i, 5, iout] = v[2]
            saved_data[i, 6, iout] = v[3]
        end
    end

    v_half = update_velocity_half(v, r, dt, t0, p, alg)

    for it in 1:nt
        t_prev = t0 + (it - 1) * dt
        r, v_half = advance_boris(v_half, r, dt, t_prev, p, alg)
        t_current = t0 + it * dt

        if it < nt || save_end
            v_node = update_velocity_node(v_half, r, dt, t_current, p, alg)
            iout += 1
            @inbounds begin
                saved_data[i, 1, iout] = r[1]
                saved_data[i, 2, iout] = r[2]
                saved_data[i, 3, iout] = r[3]
                saved_data[i, 4, iout] = v_node[1]
                saved_data[i, 5, iout] = v_node[2]
                saved_data[i, 6, iout] = v_node[3]
            end
        end
    end
end
