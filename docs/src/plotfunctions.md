# Plot Functions
Before using the plot recipes of `Testparticle.jl`, you should install the package `TestParticleMakie.jl` firstly. All plot recipes are based on `Makie.jl`.

## Basic usage
The simplest usage is directly calling the `plot` or `lines` function provided by `Makie`. For example, 
```julia
plot(sol)
```
Without other arguments, it will plot the 3d orbit of this particle. And other figures are supported, you can use an argument named `vars` to choose the variables will be plotted. The basic form for `vars` is Tuple. For example, if you want to plot the variable `x` as a function of time or the orbit, you can use it like this:
```julia
plot(sol, vars=(0, 1))  # 0 represents time, 1 represents x
plot(sol, vars=(1, 2, 3))
```
If you want to plot a function of time, position or velocity, you can define a function firstly. The arguments of the function must be a 7-dimensional array, the elements of which is corresponding to the phase space coordinates and time at a certain time respectively, and the return value must be a Number. For example,
```julia
Eₖ(xu) = mₑ*(xu[4]^2 + xu[5]^2 + xu[6]^2)/2
lines(sol, vars=(0, Eₖ))
```
This command will plot the kinetic energy as a function of x. For simplicity, the above two forms are equivalent to
```julia
lines(sol, vars=Eₖ)
lines(sol, vars=1)
```
When the independent variable is not explicit, it will be set as time automatically.

Because of the restriction of `Makie.jl`, `plot` and `lines` function cannot plot multiple lines at one call. If you want to plot multiple lines in a figure, you can use `plot!`, `lines!` or the another plot recipe `orbit`.

If you want to choose the plotting time span, you can use the argument `tspan`. For example,
```julia
lines(sol, tspan=(0, 1))
```

## Advanced usage
There are two plot recipes, which are more powerful than just using `plot` and `lines`. The `orbit` is designed to observe the orbit of a particle or the variation of other variables. It can vary time interactively. The `monitor` is designed to monitor the variation of some physics quantity when the evolution of orbit.

### `orbit`
The basic usages are same as above and just need to change the function name to `orbit`. And multiple lines can be plotted at one time by this function. For example,
```julia
orbit(sol, vars=[1, (1, 2), Eₖ])
```
If a tuple contains lists, they will be expanded automatically. For example,
```julia
orbit(sol, vars=(0, [1, 2, 3]))
orbit(sol, vars=([1, 2, 3], [4, 5, 6]))
```
are equivalent to
```julia
orbit(sol, vars=[(0, 1), (0, 2), (0, 3)])
orbit(sol, vars=[(1, 4), (2, 5), (3, 6)])
```

#### Interactive components
![](../figures/orbit_example.png)

The slider can control the time span, and the maximum of time span will be display in the right of the slider.

The `reset` button can reset the scale of lines when the axis limits are changed. Or when you drag the slider, clicking `reset` button will reset the axis limits to fit data.

### `monitor`
The usage of `vars` for this function is different from those in `orbit`. It can only be a list and its length is equal to 3. The type of its elements must be Integer or Function. The form of function is same as those prescribed by `orbit`. For example,
```julia
Eₖ(xu) = mₑ*(xu[4]^2 + xu[5]^2 + xu[6]^2)/2
monitor(sol, vars=[1, 2, Eₖ])
```

#### Interactive components
![](../figures/monitor_example.png)

After first click of the `run` button, the evolution of orbit will be displayed from the beginning. For other time, it will start from the time set by the time slider. The functionality of `reset` button is same as above.

The time slider can control the time span. The speed slider can control the speed of the animation.