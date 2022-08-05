# interactive plots

function orbit(sol::AbstractODESolution; vars=nothing, tspan=nothing, to_3d::Bool=false, interactive::Bool=true)
    if vars === nothing 
        vars = [(1, 2, 3)]
    elseif isa(vars, Tuple)
        if isa(vars[1], AbstractArray)
            if isa(vars[2], AbstractArray)
                vars = collect(zip(vars[1], vars[2]))
            elseif isa(vars[2], Integer) | isa(vars[2], Function)
                vars = [(var, vars[2]) for var in vars[1]]
            else
                throw(ArgumentError("Invalid variable type."))
            end
        else
            if isa(vars[2], AbstractArray)
                vars = [(vars[1], var) for var in vars[2]]
            elseif isa(vars[2], Integer) | isa(vars[2], Function)
                vars = [vars]
            else
                throw(ArgumentError("Invalid variable type."))
            end
        end
    elseif isa(vars, AbstractArray)
        nothing
    else
        vars = [vars]
    end

    fig = Figure()
    if isempty([var for var in vars if isa(var, Tuple) && (length(var) == 3)])
        ax = Axis(fig[1, 1:2])
        dim = 2
    else
        ax = LScene(fig[1, 1])
        dim = 3
        # ax = Axis(fig[1, 1:2], aspect=:data)
    end

    if interactive
        t_min, t_max = get_tspan(sol, tspan)
        step = (t_max - t_min) / 1e3
        time = SliderGrid(fig[2, 1], (label="Time", range=range(t_min+step, stop=t_max, step=step), format = "{:.3f}s", startvalue = t_max))
        tspan = lift(time.sliders[1].value) do t
            (t_min, t)
        end
        if dim == 2
            reset_button = Button(fig[2, 2], label="reset")
            on(reset_button.clicks) do n
                reset_limits!(ax)
            end
        end
    end

    for var in vars
        label = "(" * get_label(var) * ")"
        lines!(ax, sol, vars=var, tspan=tspan, to_3d=to_3d, label=label)
    end
    axislegend(ax)
    return fig
end


function monitor(sol::AbstractODESolution; vars=nothing, tspan=nothing)
    if vars === nothing 
        vars = [4, 5, 6]
    end
    @assert isa(vars, AbstractArray)
    @assert length(vars) == 3


    fig = Figure()
    ax1 = Axis3(fig[1:3, 1:3], aspect=:data)
    axs = [Axis(fig[i, 4:6], ylabel=get_label(x), xlabel=get_label(0)) for (i, x) in enumerate(vars)]

    t_min, t_max = get_tspan(sol, tspan)
    step = (t_max - t_min) / 1e3

    time = SliderGrid(fig[5, 1:6], (label="Time", range=range(t_min+step, stop=t_max, step=step), format = "{:.3f}s", startvalue = t_max))
    frame = SliderGrid(fig[4, 1:4], (label="Speed", range=range(1, stop=500, step=1), format = "{:d} frames/s", startvalue = 100))
    tspan = lift(time.sliders[1].value) do t
        autolimits!.(axs)
        # autolimits!(ax1)
        (t_min, t)
    end

    reset_button = Button(fig[4, 5], label="reset")
    run_button = Button(fig[4, 6], label="run")
    isrunning = Observable(false)
    first_time = true
    on(reset_button.clicks) do n
        reset_limits!(ax1)
        reset_limits!.(axs)
    end

    on(run_button.clicks) do n; isrunning[] = !isrunning[]; end
    on(run_button.clicks) do n
        t = time.sliders[1]
        if first_time
            t.value[] = t_min+step
            first_time = !first_time
        end
        @async while isrunning[] && (t.value[] < t_max)
            t.value[] += step
            sleep(1/frame.sliders[1].value[])
            reset_limits!(ax1)
            isopen(ax1.scene) || break
        end
    end

    lines!(ax1, sol, vars=(1, 2, 3), tspan=tspan)
    for (i, var) in enumerate(vars)
        @assert typeof(var) <: Union{Integer, Function}
        lines!(axs[i], sol, vars=var, tspan=tspan)
    end

    fig
end