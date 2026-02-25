"""Python wrapper for TestParticle.jl."""

from juliacall import Main as jl

# Load the Julia package
jl.seval("using TestParticle")

# Re-export Julia functions and types
prepare = jl.TestParticle.prepare
prepare_gc = jl.TestParticle.prepare_gc
get_gc = jl.TestParticle.get_gc
get_gc_func = jl.TestParticle.get_gc_func

# Functions with ! in Julia are accessed with _b suffix in Python via juliacall
trace_b = jl.TestParticle.trace_b
trace_relativistic_b = jl.TestParticle.trace_relativistic_b
trace_normalized_b = jl.TestParticle.trace_normalized_b
trace_relativistic_normalized_b = jl.TestParticle.trace_relativistic_normalized_b

trace = jl.TestParticle.trace
trace_relativistic = jl.TestParticle.trace_relativistic
trace_normalized = jl.TestParticle.trace_normalized
trace_relativistic_normalized = jl.TestParticle.trace_relativistic_normalized

trace_gc_b = jl.TestParticle.trace_gc_b
trace_gc_drifts_b = jl.TestParticle.trace_gc_drifts_b
trace_gc_flr_b = jl.TestParticle.trace_gc_flr_b
trace_gc_exb_b = jl.TestParticle.trace_gc_exb_b
trace_fieldline_b = jl.TestParticle.trace_fieldline_b
trace_fieldline = jl.TestParticle.trace_fieldline

get_gc_velocity = jl.TestParticle.get_gc_velocity
full_to_gc = jl.TestParticle.full_to_gc
gc_to_full = jl.TestParticle.gc_to_full

Proton = jl.TestParticle.Proton
Electron = jl.TestParticle.Electron
Ion = jl.TestParticle.Ion

Maxwellian = jl.TestParticle.Maxwellian
BiMaxwellian = jl.TestParticle.BiMaxwellian
Kappa = jl.TestParticle.Kappa
BiKappa = jl.TestParticle.BiKappa

AdaptiveBoris = jl.TestParticle.AdaptiveBoris
AdaptiveHybrid = jl.TestParticle.AdaptiveHybrid

CurrentLoop = jl.TestParticle.CurrentLoop
getB_loop = jl.TestParticle.getB_loop

get_gyrofrequency = jl.TestParticle.get_gyrofrequency
get_gyroperiod = jl.TestParticle.get_gyroperiod
get_gyroradius = jl.TestParticle.get_gyroradius
get_velocity = jl.TestParticle.get_velocity
get_energy = jl.TestParticle.get_energy
get_mean_magnitude = jl.TestParticle.get_mean_magnitude
energy2velocity = jl.TestParticle.energy2velocity
get_curvature_radius = jl.TestParticle.get_curvature_radius
get_adiabaticity = jl.TestParticle.get_adiabaticity

sample_unit_sphere = jl.TestParticle.sample_unit_sphere
get_number_density_flux = jl.TestParticle.get_number_density_flux

getB_zpinch = jl.TestParticle.getB_zpinch
getB_bottle = jl.TestParticle.getB_bottle
getB_mirror = jl.TestParticle.getB_mirror
getB_tokamak_coil = jl.TestParticle.getB_tokamak_coil

orbit = jl.TestParticle.orbit
monitor = jl.TestParticle.monitor

get_fields = jl.TestParticle.get_fields
get_work = jl.TestParticle.get_work

LazyTimeInterpolator = jl.TestParticle.LazyTimeInterpolator

TraceProblem = jl.TestParticle.TraceProblem
TraceGCProblem = jl.TestParticle.TraceGCProblem
TraceHybridProblem = jl.TestParticle.TraceHybridProblem

CartesianGrid = jl.TestParticle.CartesianGrid
RectilinearGrid = jl.TestParticle.RectilinearGrid
StructuredGrid = jl.TestParticle.StructuredGrid


__version__ = "0.1.0"

__all__ = [
    "prepare", "prepare_gc", "get_gc", "get_gc_func",
    "trace_b", "trace_relativistic_b", "trace_normalized_b", "trace_relativistic_normalized_b",
    "trace", "trace_relativistic", "trace_normalized", "trace_relativistic_normalized",
    "trace_gc_b",
    "trace_gc_drifts_b", "trace_gc_flr_b", "trace_gc_exb_b", "trace_fieldline_b", "trace_fieldline",
    "get_gc_velocity", "full_to_gc", "gc_to_full",
    "Proton", "Electron", "Ion",
    "Maxwellian", "BiMaxwellian", "Kappa", "BiKappa",
    "AdaptiveBoris", "AdaptiveHybrid",
    "CurrentLoop", "getB_loop",
    "get_gyrofrequency",
    "get_gyroperiod", "get_gyroradius", "get_velocity", "get_energy", "get_mean_magnitude",
    "energy2velocity", "get_curvature_radius", "get_adiabaticity",
    "sample_unit_sphere", "get_number_density_flux",
    "getB_zpinch", "getB_bottle", "getB_mirror", "getB_tokamak_coil",
    "orbit", "monitor",
    "get_fields", "get_work",
    "LazyTimeInterpolator",
    "TraceProblem", "TraceGCProblem", "TraceHybridProblem",
    "CartesianGrid", "RectilinearGrid", "StructuredGrid",
    "__version__", "jl"
]
