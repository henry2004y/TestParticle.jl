"""
    CachedFields(q2m, E, B, r, t)

A parameter container whose field functions return the values already taken at
`(r, t)`.

A step needs the fields at one point only, but a solver asks for that point
twice: once to synchronise the state at the end of a step, and once to advance
the next one from the same node. Carrying the values across the step boundary
saves the second evaluation, and it does so without the step functions having to
know where the fields came from, since they read them in the same way.
"""
struct CachedFields{Q, E, B, R, T}
    q2m::Q
    E::Returns{E}
    B::Returns{B}
    r::R
    t::T
end

@inline function _node_fields(p, r, t)
    T = eltype(r)
    E = SVector{3, T}(p.E(r, t))
    B = SVector{3, T}(p.B(r, t))
    return CachedFields(p.q2m, Returns(E), Returns(B), r, t)
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
