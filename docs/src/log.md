# Development Log

I want to take advantage of existing powerful packages as much as possible.
In order to do this, the functions provides by this package aim at working seamlessly with other solid packages.

* Less is more. Functions that are not directly relevant to tracing should not be placed in the main module.

* Wrapperless: there is no need to build wrappers around well-documented functions from other packages; you just need to provide your specific methods. Simplicity is key.

* Plot.jl is really bad at 3D related plotting. I would hope Makie.jl performs better.

* The current API for including a external force term like gravity and pressure gradient is ok, but can be made more elegant.

* Future issues will be mentioned in [GitHub Issues](https://github.com/henry2004y/TestParticle.jl/issues).