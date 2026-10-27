module TestParticleC

using BorisPushers: update_velocity_boris, update_velocity_multistep
using StaticArrays: SVector
using PrecompileTools: @setup_workload, @compile_workload

@inline function _update_velocity_multistep(
        v::SVector{3, T}, E::SVector{3, T}, B::SVector{3, T},
        q2m::T, dt::T, n::Int, ::Val{N}
    ) where {T, N}
    param = (q2m = q2m, E = (r, t) -> E, B = (r, t) -> B)
    return update_velocity_multistep(
        v, SVector{3, T}(0, 0, 0), dt, zero(T), n, Val{N}(), param
    )
end

@inline function _multistep_dispatch(v, E, B, q2m, dt, n::Int, order::Integer)
    if order == 4
        return _update_velocity_multistep(v, E, B, q2m, dt, n, Val{4}())
    elseif order == 6
        return _update_velocity_multistep(v, E, B, q2m, dt, n, Val{6}())
    else
        return _update_velocity_multistep(v, E, B, q2m, dt, n, Val{2}())
    end
end

Base.@ccallable function testparticle_update_velocity_boris_f64(
        v::Ptr{Cdouble},
        E::Ptr{Cdouble},
        B::Ptr{Cdouble},
        qdt2m::Cdouble
    )::Cvoid
    v_in = SVector{3, Float64}(unsafe_load(v, 1), unsafe_load(v, 2), unsafe_load(v, 3))
    E_vec = SVector{3, Float64}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, Float64}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))

    v_out = update_velocity_boris(v_in, E_vec, B_vec, qdt2m)
    unsafe_store!(v, v_out[1], 1)
    unsafe_store!(v, v_out[2], 2)
    unsafe_store!(v, v_out[3], 3)
    return nothing
end

Base.@ccallable function testparticle_update_velocity_boris_f32(
        v::Ptr{Cfloat},
        E::Ptr{Cfloat},
        B::Ptr{Cfloat},
        qdt2m::Cfloat
    )::Cvoid
    v_in = SVector{3, Float32}(unsafe_load(v, 1), unsafe_load(v, 2), unsafe_load(v, 3))
    E_vec = SVector{3, Float32}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, Float32}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))

    v_out = update_velocity_boris(v_in, E_vec, B_vec, qdt2m)
    unsafe_store!(v, v_out[1], 1)
    unsafe_store!(v, v_out[2], 2)
    unsafe_store!(v, v_out[3], 3)
    return nothing
end

Base.@ccallable function testparticle_boris_step_staggered_f64(
        x::Ptr{Cdouble},
        v_half::Ptr{Cdouble},
        q2m::Cdouble,
        dt::Cdouble,
        E::Ptr{Cdouble},
        B::Ptr{Cdouble}
    )::Cvoid
    r = SVector{3, Float64}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vh = SVector{3, Float64}(
        unsafe_load(v_half, 1), unsafe_load(v_half, 2), unsafe_load(v_half, 3)
    )
    E_vec = SVector{3, Float64}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, Float64}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))

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

Base.@ccallable function testparticle_boris_step_staggered_f32(
        x::Ptr{Cfloat},
        v_half::Ptr{Cfloat},
        q2m::Cfloat,
        dt::Cfloat,
        E::Ptr{Cfloat},
        B::Ptr{Cfloat}
    )::Cvoid
    r = SVector{3, Float32}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vh = SVector{3, Float32}(
        unsafe_load(v_half, 1), unsafe_load(v_half, 2), unsafe_load(v_half, 3)
    )
    E_vec = SVector{3, Float32}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, Float32}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))

    vh_next = update_velocity_boris(vh, E_vec, B_vec, q2m * dt / 2.0f0)
    r_next = r + vh_next * dt

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)
    unsafe_store!(v_half, vh_next[1], 1)
    unsafe_store!(v_half, vh_next[2], 2)
    unsafe_store!(v_half, vh_next[3], 3)
    return nothing
end

Base.@ccallable function testparticle_boris_push_f64(
        x::Ptr{Cdouble},
        v::Ptr{Cdouble},
        q2m::Cdouble,
        dt::Cdouble,
        E::Ptr{Cdouble},
        B::Ptr{Cdouble}
    )::Cvoid
    r = SVector{3, Float64}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vn = SVector{3, Float64}(unsafe_load(v, 1), unsafe_load(v, 2), unsafe_load(v, 3))
    E_vec = SVector{3, Float64}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, Float64}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))

    v_half = update_velocity_boris(vn, E_vec, B_vec, -q2m * dt / 4)
    v_half_next = update_velocity_boris(v_half, E_vec, B_vec, q2m * dt / 2)
    r_next = r + v_half_next * dt
    vn_next = update_velocity_boris(v_half_next, E_vec, B_vec, q2m * dt / 4)

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)
    unsafe_store!(v, vn_next[1], 1)
    unsafe_store!(v, vn_next[2], 2)
    unsafe_store!(v, vn_next[3], 3)
    return nothing
end

Base.@ccallable function testparticle_boris_push_f32(
        x::Ptr{Cfloat},
        v::Ptr{Cfloat},
        q2m::Cfloat,
        dt::Cfloat,
        E::Ptr{Cfloat},
        B::Ptr{Cfloat}
    )::Cvoid
    r = SVector{3, Float32}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vn = SVector{3, Float32}(unsafe_load(v, 1), unsafe_load(v, 2), unsafe_load(v, 3))
    E_vec = SVector{3, Float32}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, Float32}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))

    v_half = update_velocity_boris(vn, E_vec, B_vec, -q2m * dt / 4.0f0)
    v_half_next = update_velocity_boris(v_half, E_vec, B_vec, q2m * dt / 2.0f0)
    r_next = r + v_half_next * dt
    vn_next = update_velocity_boris(v_half_next, E_vec, B_vec, q2m * dt / 4.0f0)

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)
    unsafe_store!(v, vn_next[1], 1)
    unsafe_store!(v, vn_next[2], 2)
    unsafe_store!(v, vn_next[3], 3)
    return nothing
end

Base.@ccallable function testparticle_boris_push_batch_f64(
        n::Csize_t,
        x::Ptr{Cdouble},
        v::Ptr{Cdouble},
        q2m::Cdouble,
        dt::Cdouble,
        E::Ptr{Cdouble},
        B::Ptr{Cdouble}
    )::Cvoid
    for i in 0:(n - 1)
        offset = 3 * i
        testparticle_boris_push_f64(
            x + offset * sizeof(Cdouble),
            v + offset * sizeof(Cdouble),
            q2m,
            dt,
            E + offset * sizeof(Cdouble),
            B + offset * sizeof(Cdouble)
        )
    end
    return nothing
end

testparticle_boris_push_batch_f64(n::Integer, x, v, q2m, dt, E, B) =
    testparticle_boris_push_batch_f64(Csize_t(n), x, v, q2m, dt, E, B)

Base.@ccallable function testparticle_boris_push_batch_f32(
        n::Csize_t,
        x::Ptr{Cfloat},
        v::Ptr{Cfloat},
        q2m::Cfloat,
        dt::Cfloat,
        E::Ptr{Cfloat},
        B::Ptr{Cfloat}
    )::Cvoid
    for i in 0:(n - 1)
        offset = 3 * i
        testparticle_boris_push_f32(
            x + offset * sizeof(Cfloat),
            v + offset * sizeof(Cfloat),
            q2m,
            dt,
            E + offset * sizeof(Cfloat),
            B + offset * sizeof(Cfloat)
        )
    end
    return nothing
end

testparticle_boris_push_batch_f32(n::Integer, x, v, q2m, dt, E, B) =
    testparticle_boris_push_batch_f32(Csize_t(n), x, v, q2m, dt, E, B)

Base.@ccallable function testparticle_multistep_boris_push_f64(
        x::Ptr{Cdouble},
        v::Ptr{Cdouble},
        q2m::Cdouble,
        dt::Cdouble,
        E::Ptr{Cdouble},
        B::Ptr{Cdouble},
        n::Cint,
        order::Cint
    )::Cvoid
    r = SVector{3, Float64}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vn = SVector{3, Float64}(unsafe_load(v, 1), unsafe_load(v, 2), unsafe_load(v, 3))
    E_vec = SVector{3, Float64}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, Float64}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))
    n_sub = Int(max(1, n))

    vh = _multistep_dispatch(vn, E_vec, B_vec, q2m, -dt / 2, n_sub, order)
    vh_next = _multistep_dispatch(vh, E_vec, B_vec, q2m, dt, n_sub, order)
    r_next = r + vh_next * dt
    vn_next = _multistep_dispatch(vh_next, E_vec, B_vec, q2m, dt / 2, n_sub, order)

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)
    unsafe_store!(v, vn_next[1], 1)
    unsafe_store!(v, vn_next[2], 2)
    unsafe_store!(v, vn_next[3], 3)
    return nothing
end

Base.@ccallable function testparticle_multistep_boris_push_f32(
        x::Ptr{Cfloat},
        v::Ptr{Cfloat},
        q2m::Cfloat,
        dt::Cfloat,
        E::Ptr{Cfloat},
        B::Ptr{Cfloat},
        n::Cint,
        order::Cint
    )::Cvoid
    r = SVector{3, Float32}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vn = SVector{3, Float32}(unsafe_load(v, 1), unsafe_load(v, 2), unsafe_load(v, 3))
    E_vec = SVector{3, Float32}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, Float32}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))
    n_sub = Int(max(1, n))

    vh = _multistep_dispatch(vn, E_vec, B_vec, q2m, -dt / 2.0f0, n_sub, order)
    vh_next = _multistep_dispatch(vh, E_vec, B_vec, q2m, dt, n_sub, order)
    r_next = r + vh_next * dt
    vn_next = _multistep_dispatch(vh_next, E_vec, B_vec, q2m, dt / 2.0f0, n_sub, order)

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)
    unsafe_store!(v, vn_next[1], 1)
    unsafe_store!(v, vn_next[2], 2)
    unsafe_store!(v, vn_next[3], 3)
    return nothing
end

Base.@ccallable function testparticle_multistep_boris_step_staggered_f64(
        x::Ptr{Cdouble},
        v_half::Ptr{Cdouble},
        q2m::Cdouble,
        dt::Cdouble,
        E::Ptr{Cdouble},
        B::Ptr{Cdouble},
        n::Cint,
        order::Cint
    )::Cvoid
    r = SVector{3, Float64}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vh = SVector{3, Float64}(
        unsafe_load(v_half, 1), unsafe_load(v_half, 2), unsafe_load(v_half, 3)
    )
    E_vec = SVector{3, Float64}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, Float64}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))
    n_sub = Int(max(1, n))

    vh_next = _multistep_dispatch(vh, E_vec, B_vec, q2m, dt, n_sub, order)
    r_next = r + vh_next * dt

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)
    unsafe_store!(v_half, vh_next[1], 1)
    unsafe_store!(v_half, vh_next[2], 2)
    unsafe_store!(v_half, vh_next[3], 3)
    return nothing
end

Base.@ccallable function testparticle_multistep_boris_step_staggered_f32(
        x::Ptr{Cfloat},
        v_half::Ptr{Cfloat},
        q2m::Cfloat,
        dt::Cfloat,
        E::Ptr{Cfloat},
        B::Ptr{Cfloat},
        n::Cint,
        order::Cint
    )::Cvoid
    r = SVector{3, Float32}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vh = SVector{3, Float32}(
        unsafe_load(v_half, 1), unsafe_load(v_half, 2), unsafe_load(v_half, 3)
    )
    E_vec = SVector{3, Float32}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, Float32}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))
    n_sub = Int(max(1, n))

    vh_next = _multistep_dispatch(vh, E_vec, B_vec, q2m, dt, n_sub, order)
    r_next = r + vh_next * dt

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)
    unsafe_store!(v_half, vh_next[1], 1)
    unsafe_store!(v_half, vh_next[2], 2)
    unsafe_store!(v_half, vh_next[3], 3)
    return nothing
end

Base.@ccallable function testparticle_hyper_boris4_push_f64(
        x::Ptr{Cdouble},
        v::Ptr{Cdouble},
        q2m::Cdouble,
        dt::Cdouble,
        E::Ptr{Cdouble},
        B::Ptr{Cdouble},
        n::Cint
    )::Cvoid
    r = SVector{3, Float64}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vn = SVector{3, Float64}(unsafe_load(v, 1), unsafe_load(v, 2), unsafe_load(v, 3))
    E_vec = SVector{3, Float64}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, Float64}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))
    n_sub = Int(max(1, n))

    vh = _update_velocity_multistep(vn, E_vec, B_vec, q2m, -dt / 2, n_sub, Val{4}())
    vh_next = _update_velocity_multistep(vh, E_vec, B_vec, q2m, dt, n_sub, Val{4}())
    r_next = r + vh_next * dt
    vn_next = _update_velocity_multistep(vh_next, E_vec, B_vec, q2m, dt / 2, n_sub, Val{4}())

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)
    unsafe_store!(v, vn_next[1], 1)
    unsafe_store!(v, vn_next[2], 2)
    unsafe_store!(v, vn_next[3], 3)
    return nothing
end

Base.@ccallable function testparticle_hyper_boris4_push_f32(
        x::Ptr{Cfloat},
        v::Ptr{Cfloat},
        q2m::Cfloat,
        dt::Cfloat,
        E::Ptr{Cfloat},
        B::Ptr{Cfloat},
        n::Cint
    )::Cvoid
    r = SVector{3, Float32}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vn = SVector{3, Float32}(unsafe_load(v, 1), unsafe_load(v, 2), unsafe_load(v, 3))
    E_vec = SVector{3, Float32}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, Float32}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))
    n_sub = Int(max(1, n))

    vh = _update_velocity_multistep(vn, E_vec, B_vec, q2m, -dt / 2.0f0, n_sub, Val{4}())
    vh_next = _update_velocity_multistep(vh, E_vec, B_vec, q2m, dt, n_sub, Val{4}())
    r_next = r + vh_next * dt
    vn_next = _update_velocity_multistep(
        vh_next, E_vec, B_vec, q2m, dt / 2.0f0, n_sub, Val{4}()
    )

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)
    unsafe_store!(v, vn_next[1], 1)
    unsafe_store!(v, vn_next[2], 2)
    unsafe_store!(v, vn_next[3], 3)
    return nothing
end

Base.@ccallable function testparticle_hyper_boris6_push_f64(
        x::Ptr{Cdouble},
        v::Ptr{Cdouble},
        q2m::Cdouble,
        dt::Cdouble,
        E::Ptr{Cdouble},
        B::Ptr{Cdouble},
        n::Cint
    )::Cvoid
    r = SVector{3, Float64}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vn = SVector{3, Float64}(unsafe_load(v, 1), unsafe_load(v, 2), unsafe_load(v, 3))
    E_vec = SVector{3, Float64}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, Float64}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))
    n_sub = Int(max(1, n))

    vh = _update_velocity_multistep(vn, E_vec, B_vec, q2m, -dt / 2, n_sub, Val{6}())
    vh_next = _update_velocity_multistep(vh, E_vec, B_vec, q2m, dt, n_sub, Val{6}())
    r_next = r + vh_next * dt
    vn_next = _update_velocity_multistep(vh_next, E_vec, B_vec, q2m, dt / 2, n_sub, Val{6}())

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)
    unsafe_store!(v, vn_next[1], 1)
    unsafe_store!(v, vn_next[2], 2)
    unsafe_store!(v, vn_next[3], 3)
    return nothing
end

Base.@ccallable function testparticle_hyper_boris6_push_f32(
        x::Ptr{Cfloat},
        v::Ptr{Cfloat},
        q2m::Cfloat,
        dt::Cfloat,
        E::Ptr{Cfloat},
        B::Ptr{Cfloat},
        n::Cint
    )::Cvoid
    r = SVector{3, Float32}(unsafe_load(x, 1), unsafe_load(x, 2), unsafe_load(x, 3))
    vn = SVector{3, Float32}(unsafe_load(v, 1), unsafe_load(v, 2), unsafe_load(v, 3))
    E_vec = SVector{3, Float32}(unsafe_load(E, 1), unsafe_load(E, 2), unsafe_load(E, 3))
    B_vec = SVector{3, Float32}(unsafe_load(B, 1), unsafe_load(B, 2), unsafe_load(B, 3))
    n_sub = Int(max(1, n))

    vh = _update_velocity_multistep(vn, E_vec, B_vec, q2m, -dt / 2.0f0, n_sub, Val{6}())
    vh_next = _update_velocity_multistep(vh, E_vec, B_vec, q2m, dt, n_sub, Val{6}())
    r_next = r + vh_next * dt
    vn_next = _update_velocity_multistep(
        vh_next, E_vec, B_vec, q2m, dt / 2.0f0, n_sub, Val{6}()
    )

    unsafe_store!(x, r_next[1], 1)
    unsafe_store!(x, r_next[2], 2)
    unsafe_store!(x, r_next[3], 3)
    unsafe_store!(v, vn_next[1], 1)
    unsafe_store!(v, vn_next[2], 2)
    unsafe_store!(v, vn_next[3], 3)
    return nothing
end

Base.@ccallable function testparticle_multistep_boris_push_batch_f64(
        n_particles::Csize_t,
        x::Ptr{Cdouble},
        v::Ptr{Cdouble},
        q2m::Cdouble,
        dt::Cdouble,
        E::Ptr{Cdouble},
        B::Ptr{Cdouble},
        n_sub::Cint,
        order::Cint
    )::Cvoid
    for i in 0:(n_particles - 1)
        offset = 3 * i
        testparticle_multistep_boris_push_f64(
            x + offset * sizeof(Cdouble),
            v + offset * sizeof(Cdouble),
            q2m,
            dt,
            E + offset * sizeof(Cdouble),
            B + offset * sizeof(Cdouble),
            n_sub,
            order
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
        n_particles::Csize_t,
        x::Ptr{Cfloat},
        v::Ptr{Cfloat},
        q2m::Cfloat,
        dt::Cfloat,
        E::Ptr{Cfloat},
        B::Ptr{Cfloat},
        n_sub::Cint,
        order::Cint
    )::Cvoid
    for i in 0:(n_particles - 1)
        offset = 3 * i
        testparticle_multistep_boris_push_f32(
            x + offset * sizeof(Cfloat),
            v + offset * sizeof(Cfloat),
            q2m,
            dt,
            E + offset * sizeof(Cfloat),
            B + offset * sizeof(Cfloat),
            n_sub,
            order
        )
    end
    return nothing
end

testparticle_multistep_boris_push_batch_f32(
    n::Integer, x, v, q2m, dt, E, B, n_sub, order
) = testparticle_multistep_boris_push_batch_f32(
    Csize_t(n), x, v, q2m, dt, E, B, Cint(n_sub), Cint(order)
)

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
    end
end

end # module TestParticleC
