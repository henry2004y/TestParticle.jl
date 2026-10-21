# Internal

`using TestParticle` imports the tracing setup and execution APIs, end-to-end VDF
tracers, main result accessors, and plotting entry points. Diagnostics, conversions,
sampling and distribution APIs, interpolation policies, and storage types are public
but not exported. Access them through the module or an explicit import:

```julia
using TestParticle
using TestParticle: get_energy, Maxwellian, ClampExtrap

# Alternatively, use TestParticle.get_energy, TestParticle.Maxwellian, etc.
```

Both exported and public-only names are supported API. If an unqualified name from
an earlier release is no longer imported, qualify it with `TestParticle.` or add an
explicit import; its implementation remains available.

## Public APIs

```@autodocs
Modules = [TestParticle]
Private = false
```

## Private types and methods

```@autodocs
Modules = [TestParticle]
Public = false
```
