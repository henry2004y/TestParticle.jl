# Boris particle pusher

@inline ODE_DEFAULT_ISOUTOFDOMAIN(u, p, t) = false

# The Boris velocity update itself lives in BorisPushers, which owns the
# physics; `update_velocity_boris` and `update_velocity` are imported from there.

@inline function _prepare_saved_data(xv, p, t, ::Val{SaveFields}, ::Val{SaveWork}) where {SaveFields, SaveWork}
    data = xv

    # Pre-declare variables to share between blocks if both are true
    local E_field, magnetic_props

    if SaveFields
        r = get_x(xv)
        T = eltype(xv)
        E_field = SVector{3, T}(p.E(r, t))

        # We need magnetic properties for work, so if SaveWork is also true, compute them now
        if SaveWork
            magnetic_props = get_magnetic_properties(r, t, p.B)
            B_vec = SVector{3, T}(magnetic_props[1])
            data = vcat(data, E_field, B_vec)
        else
            B_vec = SVector{3, T}(p.B(r, t))
            data = vcat(data, E_field, B_vec)
        end
    end

    if SaveWork
        # If SaveFields was true, we already computed E_field and magnetic_props
        if SaveFields
            work = get_work_rates(xv, p, t, magnetic_props, E_field)
        else
            work = get_work_rates(xv, p, t)
        end
        data = vcat(data, work)
    end

    return data
end

"""
    get_fields(sol::AbstractODESolution)

Return the electric and magnetic fields from the solution `sol`.
"""
function get_fields(sol::AbstractODESolution)
    p = sol.prob.p
    Efunc, Bfunc = p.E, p.B

    E = map((u, t) -> Efunc(get_x(u), t), sol.u, sol.t)
    B = map((u, t) -> Bfunc(get_x(u), t), sol.u, sol.t)

    return E, B
end

"""
    get_work(sol::AbstractODESolution)

Return the work done by the electric field from the solution `sol`.
"""
function get_work(sol::AbstractODESolution)
    p = sol.prob.p
    return map((u, t) -> get_work_rates(u, p, t), sol.u, sol.t)
end
