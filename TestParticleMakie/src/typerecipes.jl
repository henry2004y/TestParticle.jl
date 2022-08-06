# Type conversion from ODESolution generated by `TestParticle` to `Makie`.

# restrict the plot type to Lines
Makie.plottype(sol::AbstractODESolution, arg...; kw...) = Lines

# prescribe the keyword arguments for the plot function
Makie.used_attributes(::Makie.PlotFunc, ::AbstractODESolution, arg...; kw...) = (:vars, :tspan, :to_3d)

# conversion from ODESolution to Point set
function Makie.convert_arguments(P::PointBased, sol::AbstractODESolution; vars=nothing, tspan=nothing, to_3d::Bool=false)
    # calculate the time span
    t_min, t_max = get_tspan(sol, tspan)

    # the type of vars must be Tuple, Number, Function or AbstractArray
    # but it will not be handled properly when its type belong to AbstractArray
    if vars === nothing
        # default figure is the orbit
        var = (1, 2, 3)
    elseif isa(vars, Integer) | isa(vars, Function)
        var = (0, vars)
    elseif isa(vars, AbstractArray)
        # Makie cannot handle multipule lines properly.
        @warn "This function can only plot one line at a time. Only the first tuple or number in vars will be plotted."
        if isa(vars[1], Tuple)
            var = vars[1]
        elseif isa(vars[1], Integer) | isa(vars[1], Function)
            var = (0, vars[1])
        else
            throw(ArgumentError("Invalid variable type."))
        end
    elseif isa(vars, Tuple)
        var = vars
    else
        throw(ArgumentError("Invalid variable type."))
    end

    # maybe the step can be decided by user
    t_step = (t_max - t_min) / 1e4
    t = collect(t_min:t_step:t_max)
    u = sol.(t)
    # the structure of subarray: [x, y, z, vx, vy, vz, t]
    append!.(u, t)

    dim = length(var)
    points = []

    for x in var
        if isa(x, Function)
            push!(points, x.(u))
        elseif x == 0
            push!(points, t)
        elseif 1<=x<=6
            # the variable maybe a phase space coordinate
            push!(points, getindex.(u, x))
        else
            throw(ArgumentError("The dimension is out of range."))
        end
    end

    # when to_3d=true, convert 2d points to 3d, and rearrange the order
    if to_3d & (dim == 2)
        diff = setdiff([1, 2, 3], var)
        if length(diff) > 1
            @warn "There are some coordinates that do not belong to position coordinates!"
        else
            push!(points, zeros(length(t)))
            # vector permutation
            order = zeros(Int, 3)
            order[[var..., diff[1]]] .= [1, 2, 3]
            points = points[order]
            dim = 3
        end
    end

    if dim == 2
        return (Point2f.(points...),)
    elseif dim == 3
        return (Point3f.(points...),)
    end
end


# get the minimum and maximum values of the time span
function get_tspan(sol::AbstractODESolution, tspan)
    if tspan === nothing
        t_min = sol.t[1]
        t_max = sol.t[end]
    else
        @assert isa(tspan, Tuple)
        t_min = tspan[1] < sol.t[1] ? sol.t[1] : tspan[1]
        t_max = tspan[2] > sol.t[end] ? sol.t[end] : tspan[2]
    end
    return t_min, t_max
end

const var_names = ("t", "x", "y", "z", "vx", "vy", "vz")
get_label(value) = string(value)
get_label(value::Integer) = value > 6 ? string(value) : var_names[value+1]
get_label(value::Tuple) = join(get_label.(value), ", ")