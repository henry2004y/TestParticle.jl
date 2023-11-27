# interactive plots

"""
    orbit(sol::AbstractODESolution; vars=nothing, tspan=nothing, to_3d::Bool=false, interactive::Bool=true)

A plot recipe for plotting orbit or other figures related to six phase space coordinates and time.

# Arguments
- `sol::AbstractODESolution`: the solution returned by the ODE solver.

# Keywords
- `vars`: the argument used to choose variables to be plotted. Default value is [(1, 2, 3)].
- `tspan::Tuple`: the span of time to be plotted. For example, tspan = (0, 1).
- `to_3d::Bool`: whether to force the points to be plotted in 3d. If `true`, the order of
coordinates will be rearranged for plotting as 3d points.
- `interactive::Bool`: whether to show the figure in interactive mode.

## vars
`vars` can be an Integer, Function, Tuple and Array. The basic form for `vars` is Tuple.
For example, if you want to plot the variable `x` as a function of time or the orbit of a
particle, you can use it like this:
```julia
orbit(sol, vars=(0, 1))
orbit(sol, vars=(1, 2, 3))
```
If the length of the tuple is 3, it will be converted to 3d.
To plot a function of time, position or velocity, you can first define the function. The
arguments of the function must be a 7-dimensional array, the elements of which correspond to
the phase space coordinates and time, and the return value must be a Number. For example,
```julia
Eₖ(xu) = mₑ*(xu[4]^2 + xu[5]^2 + xu[6]^2)/2
orbit(sol, vars=(0, Eₖ))
```
This command will plot the kinetic energy as a function of x. The above form is equivalent to
```julia
orbit(sol, vars=Eₖ)
orbit(sol, vars=1)
```
Plotting multiple lines at one time are supported by providing a tuple of indexes. For example,
```julia
orbit(sol, vars=[1, (1, 2), Eₖ])
```

Advanced usages: if a tuple contains vectors, they will be expanded automatically.
For example,
```julia
orbit(sol, vars=(0, [1, 2, 3]))
orbit(sol, vars=([1, 2, 3], [4, 5, 6]))
```
are equivalent to
```julia
orbit(sol, vars=[(0, 1), (0, 2), (0, 3)])
orbit(sol, vars=[(1, 4), (2, 5), (3, 6)])
```
"""
function orbit(sol::AbstractODESolution; vars=nothing, tspan=nothing, to_3d::Bool=false, interactive::Bool=true)
    if vars === nothing 
        vars = [(1, 2, 3)]
    elseif isa(vars, Tuple)
        if isa(vars[1], AbstractArray)
            if isa(vars[2], AbstractArray)
                # for example, ([1, 2, 3], [4, 5, 6]) -> [(1, 4), (2, 5), (3, 6)]
                vars = collect(zip(vars[1], vars[2]))
            elseif isa(vars[2], Integer) | isa(vars[2], Function)
                # for example, ([1, 3], 2) -> [(1, 2), (3, 2)]
                vars = [(var, vars[2]) for var in vars[1]]
            else
                throw(ArgumentError("Invalid variable type."))
            end
        else
            if isa(vars[2], AbstractArray)
                # for example, (1, [4, 5, 6]) -> [(1, 4), (1, 5), (1, 6)]
                vars = [(vars[1], var) for var in vars[2]]
            elseif isa(vars[2], Integer) | isa(vars[2], Function)
                # for example, (1, 2) -> [(1, 2)]
                vars = [vars]
            else
                throw(ArgumentError("Invalid variable type."))
            end
        end
    elseif isa(vars, AbstractArray)
        nothing
    else
        # The validity of vars will be checked in the `convert_arguments` function.
        vars = [vars]
    end

    fig = Figure()
    if isempty([var for var in vars if isa(var, Tuple) && (length(var) == 3)])
        ax = Axis(fig[1, 1:2])
        dim = 2
    else
        # LScene can be zoomed, but Axis3 cannot do the same thing.
        ax = LScene(fig[1, 1])
        dim = 3
        # ax = Axis3(fig[1, 1:2], aspect=:data)
    end

    if interactive
        t1, t2 = get_tspan(sol, tspan)
        step = (t2 - t1) / 1e3
        # Makie does not support the flag 'g' or 'G'.
        # https://docs.julialang.org/en/v1/stdlib/Printf/#man-printf
        time_format = t2 < 2e-3 ? "{:.3e} s" : "{:.3f} s"
        # the minimum of the range cannot be t1.
        time = SliderGrid(fig[2, 1],
            (label="Time",
             range=range(t1+step, stop=t2, step=step),
             format=time_format,
             startvalue=t2))
        # convert tspan to an Observable
        tspan = lift(time.sliders[1].value) do t
            (t1, t)
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
        lines!(ax, sol; vars=var, tspan, to_3d, label)
    end
    # Axis3 cannot displays legend correctly.
    axislegend(ax)
    return fig
end


"""
    monitor(sol::AbstractODESolution; vars=nothing, tspan=nothing)

A plot recipe for monitor the orbit of a particle and other physics quantities.

# Arguments
- `sol::AbstractODESolution`: the solution returned by the ODE solver.

# Keywords
- `vars`: the argument used to choose the variables to be plotted. Default value is [4, 5, 6].
- `tspan::Tuple`: the span of time to be plotted. For example, tspan = (0, 1).

## vars
The usage of `vars` for this function is different from those in [`orbit`](@ref). It can
only be a list and its length is equal to 3. The type of its elements must be Integer or
Function. The form of a function is the same as those prescribed by `orbit`. For example,
```julia
Eₖ(xu) = mₑ*(xu[4]^2 + xu[5]^2 + xu[6]^2)/2
monitor(sol, vars=[1, 2, Eₖ])
```
"""
function monitor(sol::AbstractODESolution; vars=nothing, tspan=nothing)
    if isnothing(vars)
        # default subplots are vx, vy and vz
        vars = [4, 5, 6]
    end
    @assert isa(vars, AbstractArray)
    @assert length(vars) == 3

    fig = Figure()
    ax1 = Axis3(fig[1:3, 1:3], aspect=:data)
    axs = [Axis(fig[i, 4:6], ylabel=get_label(x), xlabel=get_label(0))
        for (i, x) in enumerate(vars)]

    t1, t2 = get_tspan(sol, tspan)
    step = (t2 - t1) / 1e3

    # Makie does not support the flag 'g' or 'G'.
    # https://docs.julialang.org/en/v1/stdlib/Printf/#man-printf
    time_format = t2 < 2e-3 ? "{:.3e} s" : "{:.3f} s"
    time = SliderGrid(fig[5, 1:6],
        (label="Time",
         range=range(t1+step, stop=t2, step=step),
         format=time_format,
         startvalue=t2))
    # a shorter time for one frame seems useless
    frame = SliderGrid(fig[4, 1:4],
        (label="Speed",
         range=range(1, stop=200, step=1),
         format="{:d} frames/s",
         startvalue=100))
    tspan = lift(time.sliders[1].value) do t
        # According to the doc of Makie, Axis3 have the method `autolimits!` but not.
        autolimits!.(axs)
        # autolimits!(ax1)
        (t1, t)
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
        # first run should start from the beginning
        if first_time
            t.value[] = t1+step
            first_time = !first_time
        end
        tend = t1 < t2 ? t2 : t1
        @async while isrunning[] && (t.value[] < tend)
            t.value[] += step
            sleep(1/frame.sliders[1].value[])
            reset_limits!(ax1)
            # ensure the window is still alive
            isopen(ax1.scene) || break
        end
    end

    lines!(ax1, sol, vars=(1, 2, 3), tspan=tspan)
    for (i, var) in enumerate(vars)
        @assert typeof(var) <: Union{Integer, Function}
        lines!(axs[i], sol, vars=var, tspan=tspan)
    end

    return fig
end
