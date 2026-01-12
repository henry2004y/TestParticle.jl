"""
    trace_fieldline(u0, p, tspan::Tuple{<:Real, <:Real}; mode::Symbol=:both, kw...)

Helper function to create `ODEProblem`s for tracing magnetic field lines.

# Arguments

  - `u0`: Initial position.
  - `p`: Parameter tuple containing the magnetic field. It can also be the magnetic field array (numerical) or function (analytical/interpolator).
  - `tspan`: Time span (arc length) for tracing. typically `(0.0, L)`.
  - `mode`: Tracing mode, one of `:forward`, `:backward`, `:both`.
  - `kw`: Keyword arguments passed to `prepare` when `p` is an array.

# Returns

  - If `mode` is `:forward` or `:backward`, returns a single `ODEProblem`.
  - If `mode` is `:both`, returns a vector of two `ODEProblem`s (forward and backward).
"""
function trace_fieldline(u0, p, tspan::Tuple; mode::Symbol = :both, kw...)
    if p isa AbstractArray
        order = get(kw, :order, 1)
        bc = get(kw, :bc, 1)
        if haskey(kw, :x) && haskey(kw, :y) && haskey(kw, :z)
            itp = getinterp(p, kw[:x], kw[:y], kw[:z], order, bc)
        elseif haskey(kw, :grid)
            itp = getinterp(p, makegrid(kw[:grid])..., order, bc)
        else
            throw(ArgumentError("Missing grid information (x, y, z) or grid object."))
        end
        field = Field(itp)
    elseif p isa AbstractField
        field = p
    elseif p isa Function
        field = Field(p)
    elseif p isa Tuple
        field = get_BField(p)
    else
        field = p
    end

    if mode == :forward
        return ODEProblem(trace_fieldline!, u0, tspan, field)
    elseif mode == :backward
        return ODEProblem(trace_fieldline!, u0, (tspan[1], -tspan[2]), field)
    elseif mode == :both
        prob_fwd = ODEProblem(trace_fieldline!, u0, tspan, field)
        prob_bwd = ODEProblem(trace_fieldline!, u0, (tspan[1], -tspan[2]), field)
        return [prob_fwd, prob_bwd]
    else
        throw(ArgumentError("Unsupported mode: $mode. Expect :forward, :backward, or :both."))
    end
end
