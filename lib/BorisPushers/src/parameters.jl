"""
    get_q2m(p)

Return the charge-to-mass ratio `q/m` from the parameter container `p`.

The Boris solvers place no constraint on the type of `p`; they only require
`get_q2m`, `get_EField` and `get_BField` to be defined for it. The default
methods below assume the layout `(q2m, m, E, B, ...)`, so any other package
carrying the fields in a different type can support the solvers by adding
methods to these three functions.
"""
get_q2m(p) = p[1]

"""
    get_EField(p)

Return the electric field function `E(x, t)` from the parameter container `p`.
"""
get_EField(p) = p[3]

"""
    get_BField(p)

Return the magnetic field function `B(x, t)` from the parameter container `p`.
"""
get_BField(p) = p[4]

"""
    CachedFields(p, E, B, r, t)

A parameter container that answers `get_q2m`, `get_EField` and `get_BField` with
the values already taken at `(r, t)`, delegating everything else to the container
`p` it wraps.

A step needs the fields at one point only, but a solver asks for that point
twice: once to synchronise the state at the end of a step, and once to advance
the next one from the same node. Carrying the values across the step boundary
saves the second evaluation, and it does so without the step functions having to
know where the fields came from, since they ask for them in the same way.
"""
struct CachedFields{P, E, B, R, T}
    p::P
    E::E
    B::B
    r::R
    t::T
end

get_q2m(f::CachedFields) = get_q2m(f.p)
get_EField(f::CachedFields) = Base.Returns(f.E)
get_BField(f::CachedFields) = Base.Returns(f.B)

@inline function _node_fields(p, r, t)
    T = eltype(r)
    E = SVector{3, T}(get_EField(p)(r, t))
    B = SVector{3, T}(get_BField(p)(r, t))
    return CachedFields(p, E, B, r, t)
end

"""
    _fields_at(fields, p, r, t)

The fields at `(r, t)`, taking the cached ones if they were taken there and
evaluating them otherwise. The state can move between steps, which a callback
does, and fields carried across such a move would belong to the wrong place.
"""
@inline function _fields_at(fields::CachedFields, p, r, t)
    return fields.r == r && fields.t == t ? fields : _node_fields(p, r, t)
end
