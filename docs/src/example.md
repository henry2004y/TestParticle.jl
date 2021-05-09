# Examples

Multiple [demonstrations](https://github.com/henry2004y/TestParticle.jl/tree/master/examples) are provided.
For all the tracing methods, we provide both an inplace version (with `!` at the end of the function name) and a non-inplace version using StaticArrays. The non-inplace version requires the initial conditions to be static a static vector. Use them at your convenience.

## Choice of numerical algorithm

By default DifferentialEquations.jl applies `Tsit5` to a general ODE problem.
However, it is not always guaranteed to work. For example, the demo case of electron tracing in the magnetic bottle with strong magnetic field is tested to work only with fixed timestep algorithms like `Euler` and the Adams-Bashforth family.
It would take you some time to figure out which algorithm works for your problem!

## Gallery

- Tracing in a Earth-like magnetic dipole field for protons at 50 keV
![](../figures/ion_trajectory_dipole.png)

- Tracing proton in a uniform EM field in z
![](../figures/ion_uniformEM.png)

- Tracing in a magnetic dipole field corresponds to a Harris current sheet with strength 20 nT and width 0.4 Re for protons at 50 keV
![](../figures/ion_trajectory_current_sheet.png)

- Tracing electron and proton in a uniform EM field (real physical parameters)
![](../figures/electron_ion_uniformEM.png)

- Tracing electrons in a uniform EM field
![](../figures/electrons_uniformEM.png)

- Tracing relativistic electron close to the speed of light in a [magnetic bottle](https://en.wikipedia.org/wiki/Magnetic_mirror#Magnetic_bottles)
![](../figures/electron_magnetic_bottle.png)

- Tracing relativistic proton in a [Tokamak](https://en.wikipedia.org/wiki/Tokamak)[^1]
![](../figures/ion_tokamak.png)
As you can see from the trajectory, this proton will escape after a certain time.


[^1]: A excellent introduction video to Tokamak can be found [here](https://www.youtube.com/watch?v=0JqBfYwQcqg) in Mandarin.