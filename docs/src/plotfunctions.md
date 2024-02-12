# Plot Functions

TestParticle.jl wraps types of `AbstractODESolution` from DifferentialEquations.jl. The plotting recipes for Plots.jl and Makie.jl work automatically for the particle tracing solutions.

Before using the plot recipes of `Testparticle.jl`, you need to import `Makie` package via `using GLMakie` or `using CairoMakie`, which depends on your choice for the backend. All plot recipes depend on `Makie.jl`.

## Convention

By convention, we use integers to represent the 7 dimensions in the input argument `vars` for all plotting methods:

| Component | Meaning |
|-----------|---------|
| 0         | time    |
| 1         | x       |
| 2         | y       |
| 3         | z       |
| 4         | vx      |
| 5         | vy      |
| 6         | vz      |

## Basic recipe

The simplest usage is directly calling the `plot` or `lines` function provided by `Makie`. For example,

```julia
plot(sol)
```

Without other arguments, it will plot the timeseries of x, y, z, vx, vy, vz. The keyword argument `idxs` can be used to select the variables to be plotted. Please refer to [Choose Variables](https://docs.sciml.ai/DiffEqDocs/stable/basics/plot/#plot_vars) for details.

If you want to plot a function of time, position or velocity, you can first define the function. For example,

```julia
Eₖ(t, vx, vy, vz) = (t, mₑ*(vx^2 + vy^2 + vz^2)/2)
lines(sol, idxs=(Eₖ, 0, 4, 5, 6))
```

will plot the kinetic energy as a function of x.

You can choose the plotting time span via the keyword argument `tspan`. For example,

```julia
lines(sol, tspan=(0, 1))
```

The plots can be customized further:

```julia
f = Figure(size=(1200, 800), fontsize=18)
ax = Axis3(f[1,1],
    title = "Particle trajectory",
    xlabel = "X",
    ylabel = "Y",
    zlabel = "Z",
)

plot!(sol)
```

Multiple particle trajectories saved as the type `EnsembleSolution` is also supported by the Makie recipe.

## Advanced recipe

!!! note
    We currently rely on the Makie recipe implemented in DiffEqBase.jl. This recipe depends on an experimental API introduced in [Makie v0.20](https://blog.makie.org/blogposts/v0.20/), which may be unstable and contains [bugs](https://github.com/MakieOrg/Makie.jl/issues/3623). Once it becomes more stabilized with bug fixes, we will recover the interactive widgets support, e.g. `orbit` and `monitor`.

### `orbit`

![](figures/orbit_example.png)

The slider can control the time span, and the maximum of time span will be displayed on the right of the slider.

The `reset` button can reset the scale of lines when the axis limits change. When you drag the slider, clicking `reset` button will reset the axis limits to fit the data.

### `monitor`

```@raw html
<video width="75%" height="auto" controls loop>
<source src="https://raw.githubusercontent.com/TCLiuu/TestParticleResource/master/videos/monitor.mp4?raw=true" type="video/mp4">
</video>
```

![](figures/monitor_example.png)

After first click of the `run` button, the evolution of the orbit will be displayed from the beginning. For other times, it will start from the time set by the time slider. The functionality of `reset` button is the same as above.

The time slider controls the time span. The speed slider controls the speed of the animation.
