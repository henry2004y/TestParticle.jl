"""
    trace_fieldline(u0, p, tspan::Tuple{<:Real, <:Real}; mode::Symbol=:both)

Helper function to create `ODEProblem`s for tracing magnetic field lines.

# Arguments
- `u0`: Initial position.
- `p`: Parameter tuple containing the magnetic field.
- `tspan`: Time span (arc length) for tracing. typically `(0.0, L)`.
- `mode`: Tracing mode, one of `:forward`, `:backward`, `:both`.

# Returns
- If `mode` is `:forward` or `:backward`, returns a single `ODEProblem`.
- If `mode` is `:both`, returns a vector of two `ODEProblem`s (forward and backward).
"""
function trace_fieldline(u0, p, tspan::Tuple{<:Real, <:Real}; mode::Symbol=:both)
   if mode == :forward
      return ODEProblem(trace_fieldline!, u0, tspan, p)
   elseif mode == :backward
      return ODEProblem(trace_fieldline!, u0, (tspan[1], -tspan[2]), p)
   elseif mode == :both
      prob_fwd = ODEProblem(trace_fieldline!, u0, tspan, p)
      prob_bwd = ODEProblem(trace_fieldline!, u0, (tspan[1], -tspan[2]), p)
      return [prob_fwd, prob_bwd]
   else
      throw(ArgumentError("Unsupported mode: $mode. Expect :forward, :backward, or :both."))
   end
end
