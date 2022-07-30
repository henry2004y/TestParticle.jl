# Plot recipes based on Makie.jl

using MakieCore
using MakieCore: PointBased 
using GeometryBasics: Point2f, Point3f
using SciMLBase: AbstractODESolution

# type recipes
MakieCore.plottype(sol::AbstractODESolution, vars=nothing, tspan=nothing) = Lines

function MakieCore.convert_arguments(P::PointBased, sol::AbstractODESolution, vars=nothing, tspan=nothing)
    if tspan === nothing
        t_min = sol.t[1]
        t_max = sol.t[end]
    else
        t_min = tspan[1]
        t_max = tspan[2]
    end

    if vars === nothing
        vars = [(1, 2, 3)]
    elseif isa(vars, Number)
        vars = [(0, vars)]
    elseif isa(vars, AbstractArray)
        var_list = []
        for var in vars
            if isa(var, tuple)
                push!(var_list, var)
            elseif isa(var, Number)
                push!(var_list, (0, var))
            else
                throw(ArgumentError("Invalid variable type."))
            end
        end
        vars = var_list
    elseif isa(vars, Tuple)
        if isa(vars[1], AbstractArray)
            if isa(vars[2], AbstractArray)
                vars = collect(zip(vars[1], vars[2]))
            elseif isa(vars[2], Number)
                vars = [(var, vars[2]) for var in vars[1]]
            else
                throw(ArgumentError("Invalid variable type."))
            end
        else
            if isa(vars[2], AbstractArray)
                vars = [(vars[1], var) for var in vars[2]]
            elseif isa(vars[2], Number)
                vars = [vars]
            else
                throw(ArgumentError("Invalid variable type."))
            end
        end
    end

    t_step = (t_max - t_min) / 1e4
    t = collect(t_min:t_step:t_max)
    u = sol.(t)
    prepend!.(u, t)

    plot_list = []

    for var in vars
        if length(var) == 2
            points = [Point2f(r[var[1]+1], r[var[2]+1]) for r in u]
            push!(plot_list, points)
        elseif length(var) == 3
            points = [Point3f(r[var[1]+1], r[var[2]+1], r[var[3]+1]) for r in u]
            push!(plot_list, points)
        else
            throw(ArgumentError("Invalid variable number."))
        end
    end

    return Tuple(plot_list)
end