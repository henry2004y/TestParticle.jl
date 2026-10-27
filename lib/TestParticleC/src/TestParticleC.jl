module TestParticleC

using BorisPushers: update_velocity_boris, update_velocity_multistep
using StaticArrays: SVector
using PrecompileTools: @setup_workload, @compile_workload

# --- Generic Internal Solver Kernels ---

@inline function _update_velocity_multistep(
        v::SVector{3, T}, E::SVector{3, T}, B::SVector{3, T},
        q2m::T, dt::T, n::Int, ::Val{N}
    ) where {T, N}
    param = (q2m = q2m, E = (r, t) -> E, B = (r, t) -> B)
    return update_velocity_multistep(
        v, SVector{3, T}(0, 0, 0), dt, zero(T), n, Val{N}(), param
    )
end

@inline function _multistep_dispatch(v, E, B, q2m, dt, n::Int, order_val)
    if order_val isa Val{4} || order_val == 4
        return _update_velocity_multistep(v, E, B, q2m, dt, n, Val{4}())
    elseif order_val isa Val{6} || order_val == 6
        return _update_velocity_multistep(v, E, B, q2m, dt, n, Val{6}())
    else
        return _update_velocity_multistep(v, E, B, q2m, dt, n, Val{2}())
    end
end

@inline function _update_velocity_boris!(
        v::Ptr{T}, E::Ptr{T}, B::Ptr{T}, qdt2m::T
    ) where {T}
    v_in = SVector{3, T}(unsafe_load(v, 1), unsafe_load(v, 2), unsafe_load(v, 3))
    E_vec = SVector{3, T}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, T}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))
    v_out = update_velocity_boris(v_in, E_vec, B_vec, qdt2m)
    unsafe_store!(v, v_out[1], 1)
    unsafe_store!(v, v_out[2], 2)
    unsafe_store!(v, v_out[3], 3)
    return nothing
end

@inline function _boris_step_staggered!(
        x::Ptr{T}, v_half::Ptr{T}, q2m::T, dt::T, E::Ptr{T}, B::Ptr{T}
    ) where {T}
    r = SVector{3, T}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vh = SVector{3, T}(
        unsafe_load(v_half, 1), unsafe_load(v_half, 2), unsafe_load(v_half, 3)
    )
    E_vec = SVector{3, T}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, T}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))

    vh_next = update_velocity_boris(vh, E_vec, B_vec, q2m * dt / 2)
    r_next = r + vh_next * dt

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)
    unsafe_store!(v_half, vh_next[1], 1)
    unsafe_store!(v_half, vh_next[2], 2)
    unsafe_store!(v_half, vh_next[3], 3)
    return nothing
end

@inline function _boris_push!(
        x::Ptr{T}, v::Ptr{T}, q2m::T, dt::T, E::Ptr{T}, B::Ptr{T}
    ) where {T}
    r = SVector{3, T}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vn = SVector{3, T}(unsafe_load(v, 1), unsafe_load(v, 2), unsafe_load(v, 3))
    E_vec = SVector{3, T}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, T}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))

    vh = update_velocity_boris(vn, E_vec, B_vec, -q2m * dt / 4)
    vh_next = update_velocity_boris(vh, E_vec, B_vec, q2m * dt / 2)
    r_next = r + vh_next * dt
    vn_next = update_velocity_boris(vh_next, E_vec, B_vec, q2m * dt / 4)

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)
    unsafe_store!(v, vn_next[1], 1)
    unsafe_store!(v, vn_next[2], 2)
    unsafe_store!(v, vn_next[3], 3)
    return nothing
end

@inline function _multistep_boris_push!(
        x::Ptr{T}, v::Ptr{T}, q2m::T, dt::T, E::Ptr{T}, B::Ptr{T},
        n::Integer, order_val
    ) where {T}
    r = SVector{3, T}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vn = SVector{3, T}(unsafe_load(v, 1), unsafe_load(v, 2), unsafe_load(v, 3))
    E_vec = SVector{3, T}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, T}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))
    n_sub = Int(max(1, n))

    vh = _multistep_dispatch(vn, E_vec, B_vec, q2m, -dt / 2, n_sub, order_val)
    vh_next = _multistep_dispatch(vh, E_vec, B_vec, q2m, dt, n_sub, order_val)
    r_next = r + vh_next * dt
    vn_next = _multistep_dispatch(vh_next, E_vec, B_vec, q2m, dt / 2, n_sub, order_val)

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)
    unsafe_store!(v, vn_next[1], 1)
    unsafe_store!(v, vn_next[2], 2)
    unsafe_store!(v, vn_next[3], 3)
    return nothing
end

@inline function _multistep_boris_step_staggered!(
        x::Ptr{T}, v_half::Ptr{T}, q2m::T, dt::T, E::Ptr{T}, B::Ptr{T},
        n::Integer, order_val
    ) where {T}
    r = SVector{3, T}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vh = SVector{3, T}(
        unsafe_load(v_half, 1), unsafe_load(v_half, 2), unsafe_load(v_half, 3)
    )
    E_vec = SVector{3, T}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, T}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))
    n_sub = Int(max(1, n))

    vh_next = _multistep_dispatch(vh, E_vec, B_vec, q2m, dt, n_sub, order_val)
    r_next = r + vh_next * dt

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)
    unsafe_store!(v_half, vh_next[1], 1)
    unsafe_store!(v_half, vh_next[2], 2)
    unsafe_store!(v_half, vh_next[3], 3)
    return nothing
end

@inline function _boris_propose_dt(q2m::T, B::Ptr{T}, safety::T) where {T}
    B_vec = SVector{3, T}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))
    Bmag = sqrt(sum(abs2, B_vec))
    (iszero(Bmag) || iszero(q2m)) && return zero(T)
    return (2 * T(π) * safety) / (abs(q2m) * Bmag)
end

@inline function _boris_resync_velocity_half!(
        v_half::Ptr{T}, q2m::T, dt_prev::T, dt_next::T, E::Ptr{T}, B::Ptr{T}
    ) where {T}
    vh = SVector{3, T}(
        unsafe_load(v_half, 1), unsafe_load(v_half, 2), unsafe_load(v_half, 3)
    )
    E_vec = SVector{3, T}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, T}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))

    v_node = update_velocity_boris(vh, E_vec, B_vec, q2m * dt_prev / 4)
    vh_new = update_velocity_boris(v_node, E_vec, B_vec, -q2m * dt_next / 4)

    unsafe_store!(v_half, vh_new[1], 1)
    unsafe_store!(v_half, vh_new[2], 2)
    unsafe_store!(v_half, vh_new[3], 3)
    return nothing
end

@inline function _adaptive_boris_step!(
        x::Ptr{T}, v_half::Ptr{T}, dt_current::T, dt_next::Ptr{T},
        q2m::T, E::Ptr{T}, B::Ptr{T}, safety::T
    ) where {T}
    r = SVector{3, T}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vh = SVector{3, T}(
        unsafe_load(v_half, 1), unsafe_load(v_half, 2), unsafe_load(v_half, 3)
    )
    E_vec = SVector{3, T}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, T}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))

    vh_next = update_velocity_boris(vh, E_vec, B_vec, q2m * dt_current / 2)
    r_next = r + vh_next * dt_current

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)

    Bmag = sqrt(sum(abs2, B_vec))
    dt_prop = (!iszero(Bmag) && !iszero(q2m)) ?
        (2 * T(π) * safety) / (abs(q2m) * Bmag) : dt_current

    if dt_prop != dt_current && dt_prop > zero(T)
        v_node = update_velocity_boris(vh_next, E_vec, B_vec, q2m * dt_current / 4)
        vh_resynced = update_velocity_boris(v_node, E_vec, B_vec, -q2m * dt_prop / 4)
        unsafe_store!(v_half, vh_resynced[1], 1)
        unsafe_store!(v_half, vh_resynced[2], 2)
        unsafe_store!(v_half, vh_resynced[3], 3)
    else
        unsafe_store!(v_half, vh_next[1], 1)
        unsafe_store!(v_half, vh_next[2], 2)
        unsafe_store!(v_half, vh_next[3], 3)
    end

    unsafe_store!(dt_next, dt_prop)
    return nothing
end

@inline function _boris_push_batch_soa!(
        n::Csize_t,
        x::Ptr{T}, y::Ptr{T}, z::Ptr{T},
        vx::Ptr{T}, vy::Ptr{T}, vz::Ptr{T},
        q2m_arr::Ptr{T}, q2m_scalar::T, dt::T,
        Ex::Ptr{T}, Ey::Ptr{T}, Ez::Ptr{T},
        Bx::Ptr{T}, By::Ptr{T}, Bz::Ptr{T}
    ) where {T}
    has_arr = q2m_arr != C_NULL
    for i in 1:Int(n)
        q2m = has_arr ? unsafe_load(q2m_arr, i) : q2m_scalar
        r = SVector{3, T}(unsafe_load(x, i), unsafe_load(y, i), unsafe_load(z, i))
        vn = SVector{3, T}(unsafe_load(vx, i), unsafe_load(vy, i), unsafe_load(vz, i))
        E_vec = SVector{3, T}(unsafe_load(Ex, i), unsafe_load(Ey, i), unsafe_load(Ez, i))
        B_vec = SVector{3, T}(unsafe_load(Bx, i), unsafe_load(By, i), unsafe_load(Bz, i))

        vh = update_velocity_boris(vn, E_vec, B_vec, -q2m * dt / 4)
        vh_next = update_velocity_boris(vh, E_vec, B_vec, q2m * dt / 2)
        r_next = r + vh_next * dt
        vn_next = update_velocity_boris(vh_next, E_vec, B_vec, q2m * dt / 4)

        unsafe_store!(x, r_next[1], i)
        unsafe_store!(y, r_next[2], i)
        unsafe_store!(z, r_next[3], i)
        unsafe_store!(vx, vn_next[1], i)
        unsafe_store!(vy, vn_next[2], i)
        unsafe_store!(vz, vn_next[3], i)
    end
    return nothing
end

@inline function _boris_step_staggered_batch_soa!(
        n::Csize_t,
        x::Ptr{T}, y::Ptr{T}, z::Ptr{T},
        vhx::Ptr{T}, vhy::Ptr{T}, vhz::Ptr{T},
        q2m_arr::Ptr{T}, q2m_scalar::T, dt::T,
        Ex::Ptr{T}, Ey::Ptr{T}, Ez::Ptr{T},
        Bx::Ptr{T}, By::Ptr{T}, Bz::Ptr{T}
    ) where {T}
    has_arr = q2m_arr != C_NULL
    for i in 1:Int(n)
        q2m = has_arr ? unsafe_load(q2m_arr, i) : q2m_scalar
        r = SVector{3, T}(unsafe_load(x, i), unsafe_load(y, i), unsafe_load(z, i))
        vh = SVector{3, T}(unsafe_load(vhx, i), unsafe_load(vhy, i), unsafe_load(vhz, i))
        E_vec = SVector{3, T}(unsafe_load(Ex, i), unsafe_load(Ey, i), unsafe_load(Ez, i))
        B_vec = SVector{3, T}(unsafe_load(Bx, i), unsafe_load(By, i), unsafe_load(Bz, i))

        vh_next = update_velocity_boris(vh, E_vec, B_vec, q2m * dt / 2)
        r_next = r + vh_next * dt

        unsafe_store!(x, r_next[1], i)
        unsafe_store!(y, r_next[2], i)
        unsafe_store!(z, r_next[3], i)
        unsafe_store!(vhx, vh_next[1], i)
        unsafe_store!(vhy, vh_next[2], i)
        unsafe_store!(vhz, vh_next[3], i)
    end
    return nothing
end

@inline function _eval_c_field!(
        field_fn::Ptr{Cvoid}, t::Float64, r_buf::Vector{Float64},
        E_buf::Vector{Float64}, B_buf::Vector{Float64}, user_data::Ptr{Cvoid}
    )
    ccall(
        field_fn, Cvoid,
        (Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid}),
        t, pointer(r_buf), pointer(E_buf), pointer(B_buf), user_data
    )
end

@inline function _eval_c_field!(
        field_fn::Ptr{Cvoid}, t::Float32, r_buf::Vector{Float32},
        E_buf::Vector{Float32}, B_buf::Vector{Float32}, user_data::Ptr{Cvoid}
    )
    ccall(
        field_fn, Cvoid,
        (Cfloat, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cfloat}, Ptr{Cvoid}),
        t, pointer(r_buf), pointer(E_buf), pointer(B_buf), user_data
    )
end

@inline function _trace_boris!(
        x::Ptr{T}, v::Ptr{T}, t_start::T, t_end::T, dt::T, q2m::T,
        field_fn::Ptr{Cvoid}, user_data::Ptr{Cvoid},
        max_steps::Csize_t, n_steps_taken::Ptr{Csize_t}, trajectory_out::Ptr{T}
    ) where {T}
    (field_fn == C_NULL || dt <= zero(T) || max_steps == 0) && return Cint(-1)

    r = SVector{3, T}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vn = SVector{3, T}(unsafe_load(v, 1), unsafe_load(v, 2), unsafe_load(v, 3))

    r_buf = zeros(T, 3)
    E_buf = zeros(T, 3)
    B_buf = zeros(T, 3)

    t = t_start
    steps = 0
    save_traj = trajectory_out != C_NULL

    while t < t_end && steps < Int(max_steps)
        r_buf[1] = r[1]
        r_buf[2] = r[2]
        r_buf[3] = r[3]

        _eval_c_field!(field_fn, t, r_buf, E_buf, B_buf, user_data)

        E_vec = SVector{3, T}(E_buf[1], E_buf[2], E_buf[3])
        B_vec = SVector{3, T}(B_buf[1], B_buf[2], B_buf[3])

        if save_traj
            base = 7 * steps
            unsafe_store!(trajectory_out, t, base + 1)
            unsafe_store!(trajectory_out, r[1], base + 2)
            unsafe_store!(trajectory_out, r[2], base + 3)
            unsafe_store!(trajectory_out, r[3], base + 4)
            unsafe_store!(trajectory_out, vn[1], base + 5)
            unsafe_store!(trajectory_out, vn[2], base + 6)
            unsafe_store!(trajectory_out, vn[3], base + 7)
        end

        dt_step = min(dt, t_end - t)
        vh = update_velocity_boris(vn, E_vec, B_vec, -q2m * dt_step / 4)
        vh_next = update_velocity_boris(vh, E_vec, B_vec, q2m * dt_step / 2)
        r = r + vh_next * dt_step
        vn = update_velocity_boris(vh_next, E_vec, B_vec, q2m * dt_step / 4)

        t += dt_step
        steps += 1
    end

    unsafe_store!(x, r[1], 1)
    unsafe_store!(x, r[2], 2)
    unsafe_store!(x, r[3], 3)
    unsafe_store!(v, vn[1], 1)
    unsafe_store!(v, vn[2], 2)
    unsafe_store!(v, vn[3], 3)

    if n_steps_taken != C_NULL
        unsafe_store!(n_steps_taken, Csize_t(steps))
    end

    return Cint(0)
end

# --- C-Callable Bindings (Float64 & Float32) ---

Base.@ccallable function testparticle_update_velocity_boris_f64(
        v::Ptr{Cdouble}, E::Ptr{Cdouble}, B::Ptr{Cdouble}, qdt2m::Cdouble
    )::Cvoid
    _update_velocity_boris!(v, E, B, qdt2m)
end

Base.@ccallable function testparticle_update_velocity_boris_f32(
        v::Ptr{Cfloat}, E::Ptr{Cfloat}, B::Ptr{Cfloat}, qdt2m::Cfloat
    )::Cvoid
    _update_velocity_boris!(v, E, B, qdt2m)
end

Base.@ccallable function testparticle_boris_step_staggered_f64(
        x::Ptr{Cdouble}, v_half::Ptr{Cdouble}, q2m::Cdouble, dt::Cdouble,
        E::Ptr{Cdouble}, B::Ptr{Cdouble}
    )::Cvoid
    _boris_step_staggered!(x, v_half, q2m, dt, E, B)
end

Base.@ccallable function testparticle_boris_step_staggered_f32(
        x::Ptr{Cfloat}, v_half::Ptr{Cfloat}, q2m::Cfloat, dt::Cfloat,
        E::Ptr{Cfloat}, B::Ptr{Cfloat}
    )::Cvoid
    _boris_step_staggered!(x, v_half, q2m, dt, E, B)
end

Base.@ccallable function testparticle_boris_push_f64(
        x::Ptr{Cdouble}, v::Ptr{Cdouble}, q2m::Cdouble, dt::Cdouble,
        E::Ptr{Cdouble}, B::Ptr{Cdouble}
    )::Cvoid
    _boris_push!(x, v, q2m, dt, E, B)
end

Base.@ccallable function testparticle_boris_push_f32(
        x::Ptr{Cfloat}, v::Ptr{Cfloat}, q2m::Cfloat, dt::Cfloat,
        E::Ptr{Cfloat}, B::Ptr{Cfloat}
    )::Cvoid
    _boris_push!(x, v, q2m, dt, E, B)
end

Base.@ccallable function testparticle_boris_push_batch_f64(
        n::Csize_t, x::Ptr{Cdouble}, v::Ptr{Cdouble}, q2m::Cdouble, dt::Cdouble,
        E::Ptr{Cdouble}, B::Ptr{Cdouble}
    )::Cvoid
    iszero(n) && return nothing
    for i in 0:(n - 1)
        offset = 3 * i
        testparticle_boris_push_f64(
            x + offset * sizeof(Cdouble), v + offset * sizeof(Cdouble),
            q2m, dt, E + offset * sizeof(Cdouble), B + offset * sizeof(Cdouble)
        )
    end
    return nothing
end

testparticle_boris_push_batch_f64(n::Integer, x, v, q2m, dt, E, B) =
    testparticle_boris_push_batch_f64(Csize_t(n), x, v, q2m, dt, E, B)

Base.@ccallable function testparticle_boris_push_batch_f32(
        n::Csize_t, x::Ptr{Cfloat}, v::Ptr{Cfloat}, q2m::Cfloat, dt::Cfloat,
        E::Ptr{Cfloat}, B::Ptr{Cfloat}
    )::Cvoid
    iszero(n) && return nothing
    for i in 0:(n - 1)
        offset = 3 * i
        testparticle_boris_push_f32(
            x + offset * sizeof(Cfloat), v + offset * sizeof(Cfloat),
            q2m, dt, E + offset * sizeof(Cfloat), B + offset * sizeof(Cfloat)
        )
    end
    return nothing
end

testparticle_boris_push_batch_f32(n::Integer, x, v, q2m, dt, E, B) =
    testparticle_boris_push_batch_f32(Csize_t(n), x, v, q2m, dt, E, B)

Base.@ccallable function testparticle_multistep_boris_push_f64(
        x::Ptr{Cdouble}, v::Ptr{Cdouble}, q2m::Cdouble, dt::Cdouble,
        E::Ptr{Cdouble}, B::Ptr{Cdouble}, n::Cint, order::Cint
    )::Cvoid
    _multistep_boris_push!(x, v, q2m, dt, E, B, n, order)
end

Base.@ccallable function testparticle_multistep_boris_push_f32(
        x::Ptr{Cfloat}, v::Ptr{Cfloat}, q2m::Cfloat, dt::Cfloat,
        E::Ptr{Cfloat}, B::Ptr{Cfloat}, n::Cint, order::Cint
    )::Cvoid
    _multistep_boris_push!(x, v, q2m, dt, E, B, n, order)
end

Base.@ccallable function testparticle_multistep_boris_step_staggered_f64(
        x::Ptr{Cdouble}, v_half::Ptr{Cdouble}, q2m::Cdouble, dt::Cdouble,
        E::Ptr{Cdouble}, B::Ptr{Cdouble}, n::Cint, order::Cint
    )::Cvoid
    _multistep_boris_step_staggered!(x, v_half, q2m, dt, E, B, n, order)
end

Base.@ccallable function testparticle_multistep_boris_step_staggered_f32(
        x::Ptr{Cfloat}, v_half::Ptr{Cfloat}, q2m::Cfloat, dt::Cfloat,
        E::Ptr{Cfloat}, B::Ptr{Cfloat}, n::Cint, order::Cint
    )::Cvoid
    _multistep_boris_step_staggered!(x, v_half, q2m, dt, E, B, n, order)
end

Base.@ccallable function testparticle_hyper_boris4_push_f64(
        x::Ptr{Cdouble}, v::Ptr{Cdouble}, q2m::Cdouble, dt::Cdouble,
        E::Ptr{Cdouble}, B::Ptr{Cdouble}, n::Cint
    )::Cvoid
    _multistep_boris_push!(x, v, q2m, dt, E, B, n, Val{4}())
end

Base.@ccallable function testparticle_hyper_boris4_push_f32(
        x::Ptr{Cfloat}, v::Ptr{Cfloat}, q2m::Cfloat, dt::Cfloat,
        E::Ptr{Cfloat}, B::Ptr{Cfloat}, n::Cint
    )::Cvoid
    _multistep_boris_push!(x, v, q2m, dt, E, B, n, Val{4}())
end

Base.@ccallable function testparticle_hyper_boris6_push_f64(
        x::Ptr{Cdouble}, v::Ptr{Cdouble}, q2m::Cdouble, dt::Cdouble,
        E::Ptr{Cdouble}, B::Ptr{Cdouble}, n::Cint
    )::Cvoid
    _multistep_boris_push!(x, v, q2m, dt, E, B, n, Val{6}())
end

Base.@ccallable function testparticle_hyper_boris6_push_f32(
        x::Ptr{Cfloat}, v::Ptr{Cfloat}, q2m::Cfloat, dt::Cfloat,
        E::Ptr{Cfloat}, B::Ptr{Cfloat}, n::Cint
    )::Cvoid
    _multistep_boris_push!(x, v, q2m, dt, E, B, n, Val{6}())
end

Base.@ccallable function testparticle_multistep_boris_push_batch_f64(
        n_particles::Csize_t, x::Ptr{Cdouble}, v::Ptr{Cdouble}, q2m::Cdouble,
        dt::Cdouble, E::Ptr{Cdouble}, B::Ptr{Cdouble}, n_sub::Cint, order::Cint
    )::Cvoid
    iszero(n_particles) && return nothing
    for i in 0:(n_particles - 1)
        offset = 3 * i
        testparticle_multistep_boris_push_f64(
            x + offset * sizeof(Cdouble), v + offset * sizeof(Cdouble),
            q2m, dt, E + offset * sizeof(Cdouble), B + offset * sizeof(Cdouble),
            n_sub, order
        )
    end
    return nothing
end

testparticle_multistep_boris_push_batch_f64(
    n::Integer, x, v, q2m, dt, E, B, n_sub, order
) = testparticle_multistep_boris_push_batch_f64(
    Csize_t(n), x, v, q2m, dt, E, B, Cint(n_sub), Cint(order)
)

Base.@ccallable function testparticle_multistep_boris_push_batch_f32(
        n_particles::Csize_t, x::Ptr{Cfloat}, v::Ptr{Cfloat}, q2m::Cfloat,
        dt::Cfloat, E::Ptr{Cfloat}, B::Ptr{Cfloat}, n_sub::Cint, order::Cint
    )::Cvoid
    iszero(n_particles) && return nothing
    for i in 0:(n_particles - 1)
        offset = 3 * i
        testparticle_multistep_boris_push_f32(
            x + offset * sizeof(Cfloat), v + offset * sizeof(Cfloat),
            q2m, dt, E + offset * sizeof(Cfloat), B + offset * sizeof(Cfloat),
            n_sub, order
        )
    end
    return nothing
end

testparticle_multistep_boris_push_batch_f32(
    n::Integer, x, v, q2m, dt, E, B, n_sub, order
) = testparticle_multistep_boris_push_batch_f32(
    Csize_t(n), x, v, q2m, dt, E, B, Cint(n_sub), Cint(order)
)

Base.@ccallable function testparticle_boris_propose_dt_f64(
        q2m::Cdouble, B::Ptr{Cdouble}, safety::Cdouble
    )::Cdouble
    _boris_propose_dt(q2m, B, safety)
end

Base.@ccallable function testparticle_boris_propose_dt_f32(
        q2m::Cfloat, B::Ptr{Cfloat}, safety::Cfloat
    )::Cfloat
    _boris_propose_dt(q2m, B, safety)
end

Base.@ccallable function testparticle_boris_resync_velocity_half_f64(
        v_half::Ptr{Cdouble}, q2m::Cdouble, dt_prev::Cdouble, dt_next::Cdouble,
        E::Ptr{Cdouble}, B::Ptr{Cdouble}
    )::Cvoid
    _boris_resync_velocity_half!(v_half, q2m, dt_prev, dt_next, E, B)
end

Base.@ccallable function testparticle_boris_resync_velocity_half_f32(
        v_half::Ptr{Cfloat}, q2m::Cfloat, dt_prev::Cfloat, dt_next::Cfloat,
        E::Ptr{Cfloat}, B::Ptr{Cfloat}
    )::Cvoid
    _boris_resync_velocity_half!(v_half, q2m, dt_prev, dt_next, E, B)
end

Base.@ccallable function testparticle_adaptive_boris_step_f64(
        x::Ptr{Cdouble}, v_half::Ptr{Cdouble}, dt_current::Cdouble,
        dt_next::Ptr{Cdouble}, q2m::Cdouble, E::Ptr{Cdouble}, B::Ptr{Cdouble},
        safety::Cdouble
    )::Cvoid
    _adaptive_boris_step!(x, v_half, dt_current, dt_next, q2m, E, B, safety)
end

Base.@ccallable function testparticle_adaptive_boris_step_f32(
        x::Ptr{Cfloat}, v_half::Ptr{Cfloat}, dt_current::Cfloat,
        dt_next::Ptr{Cfloat}, q2m::Cfloat, E::Ptr{Cfloat}, B::Ptr{Cfloat},
        safety::Cfloat
    )::Cvoid
    _adaptive_boris_step!(x, v_half, dt_current, dt_next, q2m, E, B, safety)
end

Base.@ccallable function testparticle_boris_push_batch_soa_f64(
        n::Csize_t,
        x::Ptr{Cdouble}, y::Ptr{Cdouble}, z::Ptr{Cdouble},
        vx::Ptr{Cdouble}, vy::Ptr{Cdouble}, vz::Ptr{Cdouble},
        q2m_arr::Ptr{Cdouble}, q2m_scalar::Cdouble, dt::Cdouble,
        Ex::Ptr{Cdouble}, Ey::Ptr{Cdouble}, Ez::Ptr{Cdouble},
        Bx::Ptr{Cdouble}, By::Ptr{Cdouble}, Bz::Ptr{Cdouble}
    )::Cvoid
    _boris_push_batch_soa!(
        n, x, y, z, vx, vy, vz, q2m_arr, q2m_scalar, dt, Ex, Ey, Ez, Bx, By, Bz
    )
end

Base.@ccallable function testparticle_boris_push_batch_soa_f32(
        n::Csize_t,
        x::Ptr{Cfloat}, y::Ptr{Cfloat}, z::Ptr{Cfloat},
        vx::Ptr{Cfloat}, vy::Ptr{Cfloat}, vz::Ptr{Cfloat},
        q2m_arr::Ptr{Cfloat}, q2m_scalar::Cfloat, dt::Cfloat,
        Ex::Ptr{Cfloat}, Ey::Ptr{Cfloat}, Ez::Ptr{Cfloat},
        Bx::Ptr{Cfloat}, By::Ptr{Cfloat}, Bz::Ptr{Cfloat}
    )::Cvoid
    _boris_push_batch_soa!(
        n, x, y, z, vx, vy, vz, q2m_arr, q2m_scalar, dt, Ex, Ey, Ez, Bx, By, Bz
    )
end

Base.@ccallable function testparticle_boris_step_staggered_batch_soa_f64(
        n::Csize_t,
        x::Ptr{Cdouble}, y::Ptr{Cdouble}, z::Ptr{Cdouble},
        vhx::Ptr{Cdouble}, vhy::Ptr{Cdouble}, vhz::Ptr{Cdouble},
        q2m_arr::Ptr{Cdouble}, q2m_scalar::Cdouble, dt::Cdouble,
        Ex::Ptr{Cdouble}, Ey::Ptr{Cdouble}, Ez::Ptr{Cdouble},
        Bx::Ptr{Cdouble}, By::Ptr{Cdouble}, Bz::Ptr{Cdouble}
    )::Cvoid
    _boris_step_staggered_batch_soa!(
        n, x, y, z, vhx, vhy, vhz, q2m_arr, q2m_scalar, dt, Ex, Ey, Ez, Bx, By, Bz
    )
end

Base.@ccallable function testparticle_boris_push_batch_multispecies_f64(
        n::Csize_t, x::Ptr{Cdouble}, v::Ptr{Cdouble}, q2m_arr::Ptr{Cdouble},
        dt::Cdouble, E::Ptr{Cdouble}, B::Ptr{Cdouble}
    )::Cvoid
    iszero(n) && return nothing
    for i in 0:(n - 1)
        offset = 3 * i
        q2m = unsafe_load(q2m_arr, i + 1)
        testparticle_boris_push_f64(
            x + offset * sizeof(Cdouble), v + offset * sizeof(Cdouble),
            q2m, dt, E + offset * sizeof(Cdouble), B + offset * sizeof(Cdouble)
        )
    end
    return nothing
end

Base.@ccallable function testparticle_trace_boris_f64(
        x::Ptr{Cdouble}, v::Ptr{Cdouble}, t_start::Cdouble, t_end::Cdouble,
        dt::Cdouble, q2m::Cdouble, field_fn::Ptr{Cvoid}, user_data::Ptr{Cvoid},
        max_steps::Csize_t, n_steps_taken::Ptr{Csize_t}, trajectory_out::Ptr{Cdouble}
    )::Cint
    _trace_boris!(
        x, v, t_start, t_end, dt, q2m, field_fn, user_data,
        max_steps, n_steps_taken, trajectory_out
    )
end

testparticle_trace_boris_f64(
    x, v, t_start, t_end, dt, q2m, field_fn::Base.CFunction,
    user_data, max_steps, n_steps_taken, trajectory_out
) = testparticle_trace_boris_f64(
    x, v, t_start, t_end, dt, q2m, Base.unsafe_convert(Ptr{Cvoid}, field_fn),
    user_data, max_steps, n_steps_taken, trajectory_out
)

Base.@ccallable function testparticle_trace_boris_f32(
        x::Ptr{Cfloat}, v::Ptr{Cfloat}, t_start::Cfloat, t_end::Cfloat,
        dt::Cfloat, q2m::Cfloat, field_fn::Ptr{Cvoid}, user_data::Ptr{Cvoid},
        max_steps::Csize_t, n_steps_taken::Ptr{Csize_t}, trajectory_out::Ptr{Cfloat}
    )::Cint
    _trace_boris!(
        x, v, t_start, t_end, dt, q2m, field_fn, user_data,
        max_steps, n_steps_taken, trajectory_out
    )
end

testparticle_trace_boris_f32(
    x, v, t_start, t_end, dt, q2m, field_fn::Base.CFunction,
    user_data, max_steps, n_steps_taken, trajectory_out
) = testparticle_trace_boris_f32(
    x, v, t_start, t_end, dt, q2m, Base.unsafe_convert(Ptr{Cvoid}, field_fn),
    user_data, max_steps, n_steps_taken, trajectory_out
)

# --- Workload Precompilation ---

@setup_workload begin
    @compile_workload begin
        x64 = [0.0, 0.0, 0.0]
        v64 = [1.0e5, 0.0, 0.0]
        E64 = [0.0, 1.0e-10, 0.0]
        B64 = [0.0, 0.0, 1.0e-9]
        testparticle_boris_push_f64(
            pointer(x64), pointer(v64), 1.0e8, 1.0e-9, pointer(E64), pointer(B64)
        )
        testparticle_boris_step_staggered_f64(
            pointer(x64), pointer(v64), 1.0e8, 1.0e-9, pointer(E64), pointer(B64)
        )
        testparticle_update_velocity_boris_f64(
            pointer(v64), pointer(E64), pointer(B64), 1.0e8 * 1.0e-9 / 2
        )
        testparticle_boris_push_batch_f64(
            Csize_t(1), pointer(x64), pointer(v64), 1.0e8, 1.0e-9,
            pointer(E64), pointer(B64)
        )
        testparticle_multistep_boris_push_f64(
            pointer(x64), pointer(v64), 1.0e8, 1.0e-9, pointer(E64), pointer(B64),
            Cint(2), Cint(4)
        )
        testparticle_multistep_boris_step_staggered_f64(
            pointer(x64), pointer(v64), 1.0e8, 1.0e-9, pointer(E64), pointer(B64),
            Cint(2), Cint(4)
        )
        testparticle_hyper_boris4_push_f64(
            pointer(x64), pointer(v64), 1.0e8, 1.0e-9, pointer(E64), pointer(B64), Cint(2)
        )
        testparticle_hyper_boris6_push_f64(
            pointer(x64), pointer(v64), 1.0e8, 1.0e-9, pointer(E64), pointer(B64), Cint(2)
        )
        testparticle_multistep_boris_push_batch_f64(
            Csize_t(1), pointer(x64), pointer(v64), 1.0e8, 1.0e-9,
            pointer(E64), pointer(B64), Cint(2), Cint(4)
        )
        testparticle_boris_propose_dt_f64(1.0e8, pointer(B64), 0.1)
        testparticle_boris_resync_velocity_half_f64(
            pointer(v64), 1.0e8, 1.0e-9, 2.0e-9, pointer(E64), pointer(B64)
        )
        dt_out64 = [0.0]
        testparticle_adaptive_boris_step_f64(
            pointer(x64), pointer(v64), 1.0e-9, pointer(dt_out64), 1.0e8,
            pointer(E64), pointer(B64), 0.1
        )

        x_c64 = [0.0]
        y_c64 = [0.0]
        z_c64 = [0.0]
        vx_c64 = [1.0e5]
        vy_c64 = [0.0]
        vz_c64 = [0.0]
        Ex_c64 = [0.0]
        Ey_c64 = [0.0]
        Ez_c64 = [0.0]
        Bx_c64 = [0.0]
        By_c64 = [0.0]
        Bz_c64 = [1.0]
        testparticle_boris_push_batch_soa_f64(
            Csize_t(1), pointer(x_c64), pointer(y_c64), pointer(z_c64),
            pointer(vx_c64), pointer(vy_c64), pointer(vz_c64),
            Ptr{Cdouble}(C_NULL), 1.0e8, 1.0e-9,
            pointer(Ex_c64), pointer(Ey_c64), pointer(Ez_c64),
            pointer(Bx_c64), pointer(By_c64), pointer(Bz_c64)
        )
        testparticle_boris_step_staggered_batch_soa_f64(
            Csize_t(1), pointer(x_c64), pointer(y_c64), pointer(z_c64),
            pointer(vx_c64), pointer(vy_c64), pointer(vz_c64),
            Ptr{Cdouble}(C_NULL), 1.0e8, 1.0e-9,
            pointer(Ex_c64), pointer(Ey_c64), pointer(Ez_c64),
            pointer(Bx_c64), pointer(By_c64), pointer(Bz_c64)
        )
        q2m_arr64 = [1.0e8]
        testparticle_boris_push_batch_multispecies_f64(
            Csize_t(1), pointer(x64), pointer(v64), pointer(q2m_arr64), 1.0e-9,
            pointer(E64), pointer(B64)
        )

        x32 = Float32[0.0, 0.0, 0.0]
        v32 = Float32[1.0e5, 0.0, 0.0]
        E32 = Float32[0.0, 1.0e-10, 0.0]
        B32 = Float32[0.0, 0.0, 1.0e-9]
        testparticle_boris_push_f32(
            pointer(x32), pointer(v32), 1.0f8, 1.0f-9, pointer(E32), pointer(B32)
        )
        testparticle_boris_step_staggered_f32(
            pointer(x32), pointer(v32), 1.0f8, 1.0f-9, pointer(E32), pointer(B32)
        )
        testparticle_update_velocity_boris_f32(
            pointer(v32), pointer(E32), pointer(B32), 1.0f8 * 1.0f-9 / 2.0f0
        )
        testparticle_boris_push_batch_f32(
            Csize_t(1), pointer(x32), pointer(v32), 1.0f8, 1.0f-9,
            pointer(E32), pointer(B32)
        )
        testparticle_multistep_boris_push_f32(
            pointer(x32), pointer(v32), 1.0f8, 1.0f-9, pointer(E32), pointer(B32),
            Cint(2), Cint(4)
        )
        testparticle_multistep_boris_step_staggered_f32(
            pointer(x32), pointer(v32), 1.0f8, 1.0f-9, pointer(E32), pointer(B32),
            Cint(2), Cint(4)
        )
        testparticle_hyper_boris4_push_f32(
            pointer(x32), pointer(v32), 1.0f8, 1.0f-9, pointer(E32), pointer(B32), Cint(2)
        )
        testparticle_hyper_boris6_push_f32(
            pointer(x32), pointer(v32), 1.0f8, 1.0f-9, pointer(E32), pointer(B32), Cint(2)
        )
        testparticle_multistep_boris_push_batch_f32(
            Csize_t(1), pointer(x32), pointer(v32), 1.0f8, 1.0f-9,
            pointer(E32), pointer(B32), Cint(2), Cint(4)
        )
        testparticle_boris_propose_dt_f32(1.0f8, pointer(B32), 0.1f0)
        testparticle_boris_resync_velocity_half_f32(
            pointer(v32), 1.0f8, 1.0f-9, 2.0f-9, pointer(E32), pointer(B32)
        )
        dt_out32 = Float32[0.0]
        testparticle_adaptive_boris_step_f32(
            pointer(x32), pointer(v32), 1.0f-9, pointer(dt_out32), 1.0f8,
            pointer(E32), pointer(B32), 0.1f0
        )
    end
end

end # module TestParticleC
