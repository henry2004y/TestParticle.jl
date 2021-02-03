# Development Log

I want to take advantage of existing powerful packages as much as possible.
In order to do this, the functions provides by this package aim at working seamlessly with other solid packages.

* Handling the tracing outside the given mesh is a problem. Right now I set the field values outside the simulation domain to NaN.

* Less is more. Functions that are not directly relevant to tracing should not be placed in the main module.

* For some unknown reasons, the results from analytic dipole field are not exactly the same under Linux, Mac and Windows.