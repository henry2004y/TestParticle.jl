"""
    TraceFieldlineProblem(u0, p, tspan::Tuple{<:Real, <:Real}; mode::Symbol=:both)

Helper function to create `ODEProblem`s for tracing magnetic field lines.

# Arguments

  - `u0`: Initial position.
  - `p`: Parameter tuple containing the magnetic field, or a function (analytical/interpolator).
  - `tspan`: Time span (arc length) for tracing. typically `(0.0, L)`.
  - `mode`: Tracing mode, one of `:forward`, `:backward`, `:both`.

# Returns

  - If `mode` is `:forward` or `:backward`, returns a single `ODEProblem`.
  - If `mode` is `:both`, returns a `Tuple` of two `ODEProblem`s (forward and backward).
"""
function TraceFieldlineProblem(u0, p, tspan::Tuple; mode::Symbol = :both)
    field =
    if p isa Function
        Field(p)
    elseif p isa Tuple
        get_BField(p)
    else
        p
    end

    if mode == :forward
        return ODEProblem(trace_fieldline!, u0, tspan, field)
    elseif mode == :backward
        return ODEProblem(trace_fieldline!, u0, (tspan[1], -tspan[2]), field)
    elseif mode == :both
        prob_fwd = ODEProblem(trace_fieldline!, u0, tspan, field)
        prob_bwd = ODEProblem(trace_fieldline!, u0, (tspan[1], -tspan[2]), field)
        return (prob_fwd, prob_bwd)
    else
        throw(ArgumentError("Unsupported mode: $mode. Expect :forward, :backward, or :both."))
    end
end
