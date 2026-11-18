# # Velocity Sampling
#
# This example samples initial velocities from velocity distribution functions (VDFs)
# and compares the samples with the analytic marginal PDFs. The models come from
# [PlasmaPopulations.jl](https://github.com/JuliaPlasma/PlasmaPopulations.jl): they are
# parameterized by temperatures [eV] or thermal speeds `vth = √(2T/m)`, live in the
# field-aligned frame (``\hat{z} \parallel \mathbf{B}``), and are placed in another frame
# with `Oriented(vdf, B; u)`.

using TestParticle
using PlasmaPopulations: BiMaxwellian, BiKappa, Oriented, vth_perp_para
using Gamma: gamma
using Random
using StaticArrays
using CairoMakie
CairoMakie.activate!(type = "png") #hide
import DisplayAs #hide

rng = Xoshiro(1234);

# Number of samples
N = 100_000;

# Plotting configuration
theme = Theme(fontsize = 18)
set_theme!(theme)

# ## 1. Maxwellian Distribution
# 1D projection:
# ```math
# f(v_x) \propto \exp\left(-\frac{v_x^2}{v_{th}^2}\right)
# ```

## Isotropic proton Maxwellian at 10 eV
vdf = BiMaxwellian(Proton; T = 10.0)
println(vdf)

## Sample
vs = rand(rng, vdf, N)
vx = [v[1] for v in vs]

## Theoretical PDF
vth = vdf.vth_perp
v_grid = range(-4vth, 4vth, length = 100)
pdf_maxwellian(v, vth) = exp(-v^2 / vth^2) / (sqrt(π) * vth)

## Visualization
f = Figure()
ax = Axis(f[1, 1], title = "Maxwellian", xlabel = "vx", ylabel = "PDF")
hist!(ax, vx, normalization = :pdf, bins = 50, label = "Sampled", color = (:blue, 0.5))
lines!(
    ax, v_grid, pdf_maxwellian.(v_grid, vth), label = "Theory", color = :red, linewidth = 2
)
axislegend(ax)

f = DisplayAs.PNG(f) #hide

# ## 2. Bi-Maxwellian Distribution
# ```math
# f(v_\parallel) \propto \exp\left(-\frac{v_\parallel^2}{v_{th,\parallel}^2}\right)
# ```
# 2D magnitude:
# ```math
# f(v_\perp) \propto \frac{v_\perp}{v_{th,\perp}^2}
# \exp\left(-\frac{v_\perp^2}{v_{th,\perp}^2}\right)
# ```

## Anisotropic, placed in the frame of a field along x
B = SA[1.0, 0.0, 0.0]
vdf_bi = Oriented(BiMaxwellian(Proton; T_perp = 5.0, T_para = 20.0), B)
println(vdf_bi)

## Sample
vs = rand(rng, vdf_bi, N)
vpar = [v[1] for v in vs] # B is along x
vperp = [hypot(v[2], v[3]) for v in vs]

## Theoretical PDFs
vthperp, vthpar = vth_perp_para(vdf_bi)
pdf_bi_perp(v, vth) = (2v / vth^2) * exp(-v^2 / vth^2)

## Visualization
f = Figure(size = (800, 400))
ax1 = Axis(f[1, 1], title = "Bi-Maxwellian Parallel", xlabel = "v_par", ylabel = "PDF")
hist!(ax1, vpar, normalization = :pdf, bins = 50, label = "Sampled", color = (:blue, 0.5))
v_grid_par = range(-4vthpar, 4vthpar, length = 100)
lines!(
    ax1, v_grid_par, pdf_maxwellian.(v_grid_par, vthpar),
    label = "Theory", color = :red, linewidth = 2
)

ax2 = Axis(
    f[1, 2], title = "Bi-Maxwellian Perpendicular", xlabel = "v_perp", ylabel = "PDF"
)
hist!(ax2, vperp, normalization = :pdf, bins = 50, label = "Sampled", color = (:blue, 0.5))
v_grid_perp = range(0, 4vthperp, length = 100)
lines!(
    ax2, v_grid_perp, pdf_bi_perp.(v_grid_perp, vthperp),
    label = "Theory", color = :red, linewidth = 2
)
axislegend(ax2)

f = DisplayAs.PNG(f) #hide

# ## 3. Kappa Distribution
#
# ```math
# f(\mathbf{v}) \propto \left(1 + \frac{v^2}{A}\right)^{-\kappa-1},
# \quad A = (\kappa - 3/2)\, v_{th}^2,
# ```
# so the temperature ``T = m v_{th}^2 / 2`` does not depend on ``\kappa``.
# Its 1D projection is
# ```math
# f(v_x) = \frac{\Gamma(\kappa)}{\sqrt{\pi A}\,\Gamma(\kappa - 1/2)}
# \left(1 + \frac{v_x^2}{A}\right)^{-\kappa}.
# ```

kappa = 3.0
vdf_kappa = BiKappa(Proton; T = 10.0, kappa)
println(vdf_kappa)

## Sample
vs = rand(rng, vdf_kappa, N)
vx = [v[1] for v in vs]

## Theoretical PDF (1D projection)
function pdf_kappa(v, vth, k)
    A = (k - 3 / 2) * vth^2
    return gamma(k) / (sqrt(π * A) * gamma(k - 1 / 2)) * (1 + v^2 / A)^(-k)
end

## Visualization
f = Figure()
ax = Axis(f[1, 1], title = "Kappa (kappa=3)", xlabel = "vx", ylabel = "PDF", yscale = log10)
hist!(ax, vx, normalization = :pdf, bins = 30, label = "Sampled", color = (:blue, 0.5))
lines!(
    ax, v_grid, pdf_kappa.(v_grid, vdf_kappa.vth_perp, kappa),
    label = "Theory", color = :red, linewidth = 2
)
axislegend(ax)

f = DisplayAs.PNG(f) #hide

# ## 4. Bi-Kappa Distribution
#
# Bi-Maxwellian anisotropy with power-law tails; each 1D projection has the kappa form
# above with its own `vth`.

kappa = 4.0
vdf_bikappa = BiKappa(Proton; T_perp = 5.0, T_para = 20.0, kappa)
println(vdf_bikappa)

## Sample
vs = rand(rng, vdf_bikappa, N)
vpar = [v[3] for v in vs]

## Visualization
f = Figure()
ax = Axis(
    f[1, 1], title = "Bi-Kappa Parallel (Z)",
    xlabel = "v_par", ylabel = "PDF", yscale = log10
)
hist!(ax, vpar, normalization = :pdf, bins = 30, label = "Sampled", color = (:blue, 0.5))
lines!(
    ax, v_grid_par, pdf_kappa.(v_grid_par, vdf_bikappa.vth_para, kappa),
    label = "Theory", color = :red, linewidth = 2
)
axislegend(ax)

f = DisplayAs.PNG(f) #hide
