# Kept for compatibility until the next breaking release, then removed.

export get_work, TraceHybridProblem

TraceHybridProblem(u0, tspan, p; kw...) = TraceProblem(u0, tspan, p; kw...)

PotentialField(phi, A; kw...) = PotentialField(; phi, A, kw...)
PotentialField(A; phi = nothing, kw...) = PotentialField(; phi, A, kw...)

prepare_canonical(phi, A; kw...) = prepare_canonical(; phi, A, kw...)
prepare_canonical(A; phi = nothing, kw...) = prepare_canonical(; phi, A, kw...)

get_work(sol::AbstractODESolution) = get_work_rates(sol)

get_work_rates(xv, p, t, magnetic_properties, E = nothing) =
    get_work_rates(xv, p, t; magnetic_properties, E)

function _relativistic(::Type{M}) where {M}
    N = M === FullOrbit ? RelativisticOrbit : M === Canonical ? RelativisticCanonical :
        throw(ArgumentError("$M has no relativistic equations."))
    Base.depwarn("`relativistic = true` is deprecated, use `TraceProblem{$(nameof(N))}`.", :TraceProblem)
    return N
end

function get_B_parameters(x, t, Bfunc)
    (; B, Bmag, b̂, ∇B, JB) = get_magnetic_properties(x, t, Bfunc)
    return B, Bmag, b̂, ∇B, JB
end
