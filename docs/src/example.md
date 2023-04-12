# Examples

Multiple [demonstrations](https://github.com/henry2004y/TestParticle.jl/tree/master/examples) are provided.
For all the tracing methods, we provide both an inplace version (with `!` at the end of the function name) and a non-inplace version using StaticArrays. The non-inplace version requires the initial conditions to be static a static vector. Use them at your convenience.

## Choice of numerical algorithm

By default DifferentialEquations.jl applies `Tsit5` to an ODE problem.
However, it is not always guaranteed to work. For example, the demo case of electron tracing in the magnetic bottle with strong magnetic field is tested to work only with fixed timestep algorithms like `Euler` and the Adams-Bashforth family.
Take you some time to figure out which algorithm works for your problem!

## Multiple particles tracing

There are two ways to trace multiple particles simultaneously:

1. Extracting the solution in a loop with varying initial conditions. See an example [here](https://github.com/henry2004y/TestParticle.jl/tree/master/examples/demo_multiple.jl).
2. Constructing the [Ensemble Simulations](https://diffeq.sciml.ai/stable/features/ensemble/). One example can be found [here](https://github.com/henry2004y/TestParticle.jl/tree/master/examples/demo_ensemble.jl). However, note that currently the ensemble type replicates the parameters for each solution, which is very memory inefficient for tracing in a numeric field.

## Gallery

- Tracing proton in a uniform EM field
![](figures/ion_uniformEM.png)

The electric field is parallel to the magnetic field in z direction, so the motion consists of a cyclotron gyration and an acceleration along z.

- Tracing electrons in a uniform EM field
![](figures/electrons_uniformEM.png)

- Tracing electron and proton in the same uniform EM field (real physical parameters), [source](https://github.com/henry2004y/TestParticle.jl/tree/master/examples/demo_electron_proton.jl)
![](figures/electron_ion_uniformEM.png)

Due to the fact that ``m_p / m_e \doteq 1836``, the proton gyro-radius is 1800 times larger than the electron, if they start with the same velocity as in this case. In more common cases we would compare electrons and protons with the same energy, and their gyro-radii differ by a factor of ``\sqrt{m_p/m_e} \sim 40``.

- Tracing in a magnetic dipole field corresponds to a Harris current sheet with strength 20 nT and width 0.4 Re for protons at 50 keV, [source](https://github.com/henry2004y/TestParticle.jl/tree/master/examples/demo_currentsheet.jl)
![](figures/ion_trajectory_current_sheet.png)

- Tracing in a Earth-like magnetic dipole field for protons at 50 keV, [source](https://github.com/henry2004y/TestParticle.jl/tree/master/examples/demo_proton_dipole.jl)
![](figures/ion_trajectory_dipole.png)

This is a combination of grad-B drift, curvature drift, and the bounce motion between mirror points.

- Tracing electrons close to the speed of light in a [magnetic bottle](https://en.wikipedia.org/wiki/Magnetic_mirror#Magnetic_bottles). The left and right subplots show the solutions with relativistic/non-relativistic Lorentz forces, respectively. [source](https://github.com/henry2004y/TestParticle.jl/tree/master/examples/demo_magneticbottle.jl)
![](figures/electron_magnetic_bottle.png)

- Tracing relativistic proton in a [Tokamak](https://en.wikipedia.org/wiki/Tokamak)[^1]
![](figures/ion_tokamak.png), [source](https://github.com/henry2004y/TestParticle.jl/tree/master/examples/demo_tokamak.jl)

As you can see from the trajectory, this proton will escape after a certain time.

[^1]: A excellent introduction video to Tokamak can be found [here](https://www.youtube.com/watch?v=0JqBfYwQcqg) in Mandarin.

- Tracing protons in the simulated EM field near Ganymede from MHD-EPIC
![](../figures/proton_ganymede_mhdepic.png)
