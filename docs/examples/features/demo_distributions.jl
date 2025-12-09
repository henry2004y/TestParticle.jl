# # Sampling from VDFs
#
# This example demonstrates how to sample from various velocity distribution functions (VDFs)
# including Maxwellian, Bi-Maxwellian, Kappa, Bi-Kappa, and Self-Similar distributions.
# We compare the sampled distributions with their theoretical probability density functions (PDFs).

using TestParticle
using Random
using StaticArrays
using Statistics
using SpecialFunctions: gamma
using CairoMakie
CairoMakie.activate!(type = "png") #hide
import DisplayAs #hide

Random.seed!(1234)

# Number of samples
N = 100_000

# Plotting configuration
theme = Theme(fontsize = 18)
set_theme!(theme)

# ## 1. Maxwellian Distribution
# 1D projection:
# ```math
# v \propto \exp\left(-\frac{v^2}{v_{th}^2}\right)
# ```

## Parameters
u0 = [0.0, 0.0, 0.0]
p = 1.0 # Thermal pressure
n = 1.0 # Number density
m = 1.0 # Mass

## Construct distribution
vdf = Maxwellian(u0, p, n; m)
println(vdf)

## Sample
vs = [sample(vdf) for _ in 1:N]
vx = [v[1] for v in vs]

## Theoretical PDF
vth = vdf.vth
v_grid = range(-4vth, 4vth, length = 100)
pdf_maxwellian(v, vth) = exp(-v^2 / vth^2) / (sqrt(π) * vth)

## Visualization
f = Figure()
ax = Axis(f[1, 1], title = "Maxwellian", xlabel = "vx", ylabel = "PDF")
hist!(ax, vx, normalization = :pdf, bins = 50, label = "Sampled", color = (:blue, 0.5))
lines!(
   ax, v_grid, pdf_maxwellian.(v_grid, vth), label = "Theory", color = :red, linewidth = 2)
axislegend(ax)

f = DisplayAs.PNG(f) #hide

# ## 2. Bi-Maxwellian Distribution
# ```math
# (v_\parallel) \propto \exp\left(-\frac{v_\parallel^2}{v_{th,\parallel}^2}\right)
# ```
# 2D magnitude:
# ```math
# (v_\perp) \propto \frac{v_\perp}{v_{th,\perp}^2} \exp\left(-\frac{v_\perp^2}{v_{th,\perp}^2}\right)
# ```

## Parameters
B = [0.0, 0.0, 1.0]
ppar = 2.0
pperp = 0.5

## Construct distribution
vdf_bi = BiMaxwellian(B, u0, ppar, pperp, n; m)
println(vdf_bi)

## Sample
vs = [sample(vdf_bi) for _ in 1:N]
vpar = [v[3] for v in vs] # Since B is along Z
vperp = [sqrt(v[1]^2 + v[2]^2) for v in vs]

## Theoretical PDFs
vthpar = vdf_bi.vthpar
vthperp = vdf_bi.vthperp
## Bi-Maxwellian Parallel is identical to the 1D Maxwellian
pdf_bi_par(v, vth) = pdf_maxwellian(v, vth)
## Bi-Maxwellian Perpendicular (2D Speed / Rayleigh)
pdf_bi_perp(v, vth) = (2v / vth^2) * exp(-v^2 / vth^2)

## Visualization
f = Figure(size = (800, 400))
ax1 = Axis(f[1, 1], title = "Bi-Maxwellian Parallel (Z)", xlabel = "v_par", ylabel = "PDF")
hist!(ax1, vpar, normalization = :pdf, bins = 50, label = "Sampled", color = (:blue, 0.5))
lines!(
   ax1, v_grid, pdf_bi_par.(v_grid, vthpar), label = "Theory", color = :red, linewidth = 2)

ax2 = Axis(
   f[1, 2], title = "Bi-Maxwellian Perpendicular", xlabel = "v_perp", ylabel = "PDF")
hist!(ax2, vperp, normalization = :pdf, bins = 50, label = "Sampled", color = (:blue, 0.5))
v_grid_perp = range(0, 4vthperp, length = 100)
lines!(ax2, v_grid_perp, pdf_bi_perp.(v_grid_perp, vthperp),
   label = "Theory", color = :red, linewidth = 2)
axislegend(ax2)

f = DisplayAs.PNG(f) #hide

# ## 3. Kappa Distribution
# 1D projection
# ```math
# v \propto \left(1 + \frac{v^2}{(\kappa - 1.5)v_{th}^2}\right)^{-\kappa-1}$
# ```

## Parameters
kappa = 3.0

## Construct distribution
vdf_kappa = Kappa(u0, p, n, kappa; m)
println(vdf_kappa)

## Sample
vs = [sample(vdf_kappa) for _ in 1:N]
vx = [v[1] for v in vs]

## Theoretical PDF (1D projection of 3D Kappa)
## The 1D projection of a 3D Kappa distribution with parameter kappa is a 1D Kappa with parameter kappa.
## f(v) = A * (1 + v^2 / ((kappa - 1.5) * vth^2))^(-kappa)
## Normalization constant A = Gamma(kappa) / (sqrt(pi) * vth * sqrt(kappa - 1.5) * Gamma(kappa - 0.5))
function pdf_kappa(v, vth, k)
   theta2 = (k - 1.5) * vth^2
   A = gamma(k) / (sqrt(π * theta2) * gamma(k - 0.5))
   return A * (1 + v^2 / theta2)^(-k)
end

## Visualization
f = Figure()
ax = Axis(f[1, 1], title = "Kappa (kappa=)", xlabel = "vx", ylabel = "PDF", yscale = log10)
hist!(ax, vx, normalization = :pdf, bins = 30, label = "Sampled", color = (:blue, 0.5))
lines!(ax, v_grid, pdf_kappa.(v_grid, vdf_kappa.vth, kappa),
   label = "Theory", color = :red, linewidth = 2)
axislegend(ax)

f = DisplayAs.PNG(f) #hide

# ## 4. Bi-Kappa Distribution
#
# Similar to Bi-Maxwellian but with power-law tails.

## Parameters
ppar = 1.0
pperp = 1.0
kappa = 4.0

## Construct distribution
vdf_bikappa = BiKappa(B, u0, ppar, pperp, n, kappa; m)
println(vdf_bikappa)

## Sample
vs = [sample(vdf_bikappa) for _ in 1:N]
vpar = [v[3] for v in vs]

## Theoretical PDF (1D projection)
## Similar form to isotropic Kappa 1D projection
pdf_bikappa_par(v, vth, k) = pdf_kappa(v, vth, k)

## Visualization
f = Figure()
ax = Axis(f[1, 1], title = "Bi-Kappa Parallel (Z)",
   xlabel = "v_par", ylabel = "PDF", yscale = log10)
hist!(ax, vpar, normalization = :pdf, bins = 30, label = "Sampled", color = (:blue, 0.5))
lines!(ax, v_grid, pdf_bikappa_par.(v_grid, vdf_bikappa.vthpar, kappa),
   label = "Theory", color = :red, linewidth = 2)
axislegend(ax)

f = DisplayAs.PNG(f) #hide

# ## 5. Self-Similar Distribution
# ```math
# (v) \propto \exp\left(-\left(\frac{|v|}{\alpha}\right)^s\right)
# ```
# The implementation samples independent components from a Generalized Normal distribution.
# ```math
# (v_x) \propto \exp\left(-\left(\frac{|v_x|}{\alpha}\right)^s\right)
# ```

## Parameters
s = 4.0 # Flattop
p = 1.0

## Construct distribution
vdf_ss = SelfSimilar(u0, p, n, s; m)
println(vdf_ss)

## Sample
vs = [sample(vdf_ss) for _ in 1:N]
vx = [v[1] for v in vs]

## Theoretical PDF
## Generalized Normal Distribution
## f(x) = s / (2 * alpha * Gamma(1/s)) * exp(-(|x|/alpha)^s)
## Scale alpha is calculated in the sampler such that variance is vth^2
alpha = vdf_ss.vth * sqrt(gamma(1 / s) / gamma(3 / s))

pdf_gn(x, alpha, s) = s / (2 * alpha * gamma(1 / s)) * exp(-(abs(x) / alpha)^s)

## Visualization
f = Figure()
ax = Axis(f[1, 1], title = "Self-Similar (s=)", xlabel = "vx", ylabel = "PDF")
hist!(ax, vx, normalization = :pdf, bins = 30, label = "Sampled", color = (:blue, 0.5))
lines!(ax, v_grid, pdf_gn.(v_grid, alpha, s), label = "Theory", color = :red, linewidth = 2)
axislegend(ax)

f = DisplayAs.PNG(f) #hide
