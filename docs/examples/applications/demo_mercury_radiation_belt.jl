# # Mercury's Dynamic Electron Radiation Belt
#
# This example models the trapping and solar-wind-driven loss of energetic electrons in
# Mercury's magnetosphere, inspired by
# [Zhao et al. (2026), *Dynamic electron radiation belt at Mercury controlled by solar wind driving*](https://doi.org/10.1038/s41550-026-02925-3).
#
# Because Mercury orbits the Sun on a highly eccentric ellipse ($r_{\mathrm{hel}} \approx
# 0.307$–$0.467\,\mathrm{au}$), the upstream solar wind dynamic pressure $p_{\mathrm{dyn}}
# \propto r_{\mathrm{hel}}^{-2}$ varies by more than a factor of two over each Hermean year:
#
# - **Near aphelion** ($r_{\mathrm{hel}} \approx 0.467\,\mathrm{au}$, weak $p_{\mathrm{dyn}}$),
#   dayside magnetic field lines have a single equatorial $|B|$ minimum. Energetic electrons
#   execute conventional inter-hemispheric bounce motion and remain adiabatically trapped on
#   closed drift shells.
# - **Near perihelion** ($r_{\mathrm{hel}} \approx 0.307\,\mathrm{au}$, strong $p_{\mathrm{dyn}}$),
#   compression of the dayside magnetopause creates a local equatorial $|B|$ maximum flanked
#   by dual off-equatorial minima. Drifting electrons undergo **drift-orbit bifurcation**
#   (Shabansky orbits), jumping between inter- and intra-hemispheric bounce motion. This
#   breaks the second adiabatic invariant, drives rapid radial diffusion and advection, and
#   empties the belt via **magnetopause shadowing**, **planetary precipitation** (enhanced in
#   the southern hemisphere by Mercury's $0.196\,R_M$ northward dipole offset), and
#   **nightside current-sheet curvature scattering** ($\kappa = \sqrt{R_c / r_g} < 3$).
#
# Here we combine Mercury's internal magnetic field from
# [`PlanetaryMagneticFields.jl`](https://github.com/JuliaSpacePhysics/PlanetaryMagneticFields.jl)
# (`anderson2012`) with the dynamic KT17 external magnetospheric current model
# ([Korth et al., 2015](https://doi.org/10.1002/2015JA021022);
# [Korth et al., 2017](https://doi.org/10.1002/2017GL074699)), and trace $240$–$360\,\mathrm{keV}$
# electrons with `TestParticle.jl`'s `RelativisticGuidingCenter` solver.

import DisplayAs #hide
using TestParticle, OrdinaryDiffEq, StaticArrays
import TestParticle as TP
using TestParticle: c, eV, get_magnetic_properties
using PlanetaryMagneticFields: load_model
using LinearAlgebra: norm
using Random, Statistics, Printf, Markdown
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# ## Mercury's Magnetospheric Field (`PlanetaryMagneticFields.jl` + KT17)
#
# We work in Mercury Solar Magnetospheric (**MSM**) coordinates centered on Mercury's
# magnetic dipole, which is offset northward by $z_0 = 0.196\,R_M$ ($R_M = 2440\,\mathrm{km}$)
# relative to the planetary center in MSO coordinates:
#
# ```math
# \mathbf{r}_{\mathrm{MSO}} = \mathbf{r}_{\mathrm{MSM}} + (0,\, 0,\, z_0).
# ```
#
# In the KT17 model, the subsolar magnetopause standoff distance scales with heliocentric
# distance $r_{\mathrm{hel}}$ and the Anderson et al. (2013) disturbance index
# $\mathrm{act} \in [0, 100]$ (fixed at $\mathrm{act} = 50$ for a moderate magnetosphere):
#
# ```math
# R_{ss} = (2.06873 - 0.00279\,\mathrm{act})\,r_{\mathrm{hel}}^{1/3},
# \qquad
# \kappa_{\mathrm{mp}} = \frac{R_0}{R_{ss}},
# \qquad
# R_0 = 1.42\,R_M.
# ```
#
# Because a dipole scales as $r^{-3}$, the KT17 dipole scaling
# $\kappa_{\mathrm{mp}}^3 \mathbf{B}_{\mathrm{dip}}(\kappa_{\mathrm{mp}}\mathbf{r}_{\mathrm{MSM}})$
# reduces identically to $\mathbf{B}_{\mathrm{dip}}(\mathbf{r}_{\mathrm{MSM}})$. We
# therefore evaluate Mercury's internal field directly with
# `PlanetaryMagneticFields.load_model(:mercury, "anderson2012")` (the $n = 1\text{–}4$
# zonal spherical harmonic representation of the offset dipole, agreeing with KT17's dipole
# to $< 0.4\%$ at $r = 1.2\,R_M$ and natively differentiable with `ForwardDiff`) and add the
# KT17 external shielding, tail-disk, and tail-slab currents.

const R_M = 2440.0e3          # [m] Mercury radius
const Z_OFFSET = 0.196        # [R_M] northward dipole offset
const R0_MP = 1.42            # [R_M] reference subsolar magnetopause standoff
const ALFA_MP = 0.5           # Shue (1997) magnetopause flaring parameter
const D0_CS = 0.09            # [R_M] minimal current-sheet half-thickness
const DELTADX_CS = 1.0        # sunward current-sheet expansion factor
const DELTADY_CS = 0.1        # flankward current-sheet expansion factor
const SCALEX_CS = 1.5         # [R_M] sunward expansion scale length
const SCALEY_CS = 9.0         # [R_M] flankward expansion scale length
const ZSHIFT_CS = 3.5         # [R_M] image-sheet offset

## KT17 model coefficients (Korth et al., 2015, 2017)
const F_TAILDISK = SA[
    59048.35734, -135664.4246, -913.4507339, 209989.1008, -213142.937,
]
const B_TAILDISK = SA[
    19.69235037, -18.16704312, 12.69175932, -14.13692134, 14.13449724,
]
const C_TAILDISK = SA[
    7.682813704, 9.663177797, 0.6465427021, 1.274059603, 1.280231032,
]

const C_DIPSHLD = SA[
    7.792407683, 74.37654983, 4.119647072, -131.33086,
    546.6006311, -1077.694401, 52.46268495, 1057.273707,
    -74.91550119, -141.8047123, 3.876004886, 156.2250932,
    -506.6470185, 1439.804381, -64.55225925, -1443.754088,
]
const P_DIPSHLD = SA[0.1412297078, 0.7439847555, 1.042798338, 0.7057116022]

const C_DISKSHLD = SA[
    -398.4670279, -1143.001682, -1836.300383, -73.92180417, -326.3986853, -29.96868107,
    -1157.035602, -604.1846034, -52.04876183, -2030.691236, -1529.120337, -6.382209946,
    2587.666032, 213.8979183, -28.30225993, 630.1309859, 2968.552238, 888.6328623,
    497.3863092, 2304.254471, 858.4176875, 1226.958595, 850.1684953, -20.9011094,
    -203.9184239, -792.6099018, 1115.95569, 527.3226825, 22.47634041, -0.0704405637,
    -1405.093137, -97.20408343, 5.656730182, -138.7129102, -1979.755673, 5.407603749,
]
const P_DISKSHLD = SA[
    1.091088905, 0.6733299808, 0.3266747827, 0.9533161464, 1.362763038, 0.0014515208,
]

const C_SLABSHLD = SA[
    -91.67686636, -87.31240824, 251.8848107, 95.65629983, -80.968107, 198.1447476,
    -283.1968987, -269.1514899, 504.632231, 166.027215, -214.9025413, 623.7920115,
    -35.99544615, -322.864469, 345.710579, 928.8553184, 810.129509, 19.62627762,
    -12.70326428, 490.4662048, -814.0985363, -1781.184984, -1371.261326, 60.3136479,
    116.630551, -178.3347065, 604.0308838, 1155.151174, 770.3896601, -202.8545948,
    298.6337705, 304.7964641, 33.70850254, 393.6080147, 308.1194271, -660.1691658,
]
const P_SLABSHLD = SA[
    1.677629714, 1.292226584, 0.3116253398, -0.4392669057, 0.7578074817, 1.497779521,
]

@inline function kt17_shield(
        coeffs::SVector{M, Float64}, p::SVector{N, Float64}, x, y, z
    ) where {M, N}
    T = promote_type(typeof(x), typeof(y), typeof(z))
    bx = zero(T)
    by = zero(T)
    bz = zero(T)
    sz = ntuple(k -> sin(z * @inbounds(p[k])), Val(N))
    cz = ntuple(k -> cos(z * @inbounds(p[k])), Val(N))
    @inbounds for j in 1:N
        pj = p[j]
        sypj, cypj = sincos(y * pj)
        for k in 1:N
            pk = p[k]
            szpk = sz[k]
            czpk = cz[k]
            sqpp = sqrt(pj * pj + pk * pk)
            epp = exp(x * sqpp)
            cjk = coeffs[(j - 1) * N + k]
            bx -= sqpp * epp * cypj * szpk * cjk
            by += epp * sypj * szpk * pj * cjk
            bz -= epp * cypj * czpk * pk * cjk
        end
    end
    return SA[bx, by, bz]
end

@inline function kt17_taildisk(xmsm, ymsm, zmsm)
    sc = 7.0
    x = (xmsm - 0.3) * sc
    y = ymsm * sc
    z = zmsm * sc
    d0_sc = D0_CS * sc
    dx_sc = DELTADX_CS * sc
    dy_sc = DELTADY_CS * sc

    rho = hypot(x, y)
    if rho > 1.0e-12
        inv_rho = inv(rho)
        drhodx = x * inv_rho
        drhody = y * inv_rho
    else
        drhodx = zero(x)
        drhody = zero(y)
    end

    dex = exp(x / 7.0)
    d = d0_sc + dy_sc * (y / 20.0)^2 + dx_sc * dex
    dddy = dy_sc * y * 0.005
    dddx = (dx_sc / 7.0) * dex

    dzeta = hypot(z, d)
    inv_dzeta = inv(dzeta)
    ddzetadx = d * dddx * inv_dzeta
    ddzetady = d * dddy * inv_dzeta
    ddzetadz = z * inv_dzeta

    T = promote_type(typeof(x), typeof(y), typeof(z))
    bx = zero(T)
    by = zero(T)
    bz = zero(T)

    @inbounds for i in 1:5
        fi = F_TAILDISK[i]
        bi = B_TAILDISK[i]
        ci = C_TAILDISK[i]

        s1 = hypot(rho + bi, dzeta + ci)
        s2 = hypot(rho - bi, dzeta + ci)
        inv_s1 = inv(s1)
        inv_s2 = inv(s2)

        ds1drho = (rho + bi) * inv_s1
        ds2drho = (rho - bi) * inv_s2
        ds1ddz = (dzeta + ci) * inv_s1
        ds2ddz = (dzeta + ci) * inv_s2

        ds1dx = ds1drho * drhodx + ds1ddz * ddzetadx
        ds1dy = ds1drho * drhody + ds1ddz * ddzetady
        ds1dz = ds1ddz * ddzetadz

        ds2dx = ds2drho * drhodx + ds2ddz * ddzetadx
        ds2dy = ds2drho * drhody + ds2ddz * ddzetady
        ds2dz = ds2ddz * ddzetadz

        s1ts2 = s1 * s2
        s1ps2 = s1 + s2
        s1ps2sq = s1ps2^2

        fac1 = sqrt(max(s1ps2sq - (2.0 * bi)^2, 1.0e-30))
        as = fac1 / (s1ts2 * s1ps2sq)
        dasds1 = (inv(fac1 * s2) - as / s1ps2 * (s2^2 + s1 * (3.0 * s1 + 4.0 * s2))) /
            (s1 * s1ps2)
        dasds2 = (inv(fac1 * s1) - as / s1ps2 * (s1^2 + s2 * (3.0 * s2 + 4.0 * s1))) /
            (s2 * s1ps2)

        dasdx = dasds1 * ds1dx + dasds2 * ds2dx
        dasdy = dasds1 * ds1dy + dasds2 * ds2dy
        dasdz = dasds1 * ds1dz + dasds2 * ds2dz

        bx -= fi * x * dasdz
        by -= fi * y * dasdz
        bz += fi * (2.0 * as + x * dasdx + y * dasdy)
    end
    return SA[bx, by, bz]
end

@inline function kt17_tailslab(xmsm, ymsm, zmsm)
    ex = exp(xmsm / SCALEX_CS)
    d = D0_CS + DELTADX_CS * ex + DELTADY_CS * (ymsm / SCALEY_CS)^2
    dddx = (DELTADX_CS / SCALEX_CS) * ex
    inv_d = inv(d)
    zp = zmsm + ZSHIFT_CS
    zm = zmsm - ZSHIFT_CS
    tz = tanh(zmsm * inv_d)
    tzp = tanh(zp * inv_d)
    tzm = tanh(zm * inv_d)
    bx = (tz - 0.5 * (tzm + tzp)) * inv_d
    bz = (zmsm * tz - 0.5 * (zm * tzm + zp * tzp)) * dddx * inv_d^2
    return SA[bx, zero(bx), bz]
end

"Load Mercury's `anderson2012` internal field model from PlanetaryMagneticFields.jl."
const MERCURY_INTERNAL = load_model(:mercury, "anderson2012")

struct MercuryKT17Field{M}
    internal::M
    r_hel::Float64
    rss::Float64
    κ_mp::Float64
    κ_mp3::Float64
    tamp1::Float64
    tamp2::Float64
end

function MercuryKT17Field(r_hel::Real; act::Real = 50.0, internal = MERCURY_INTERNAL)
    f = 2.06873 - 0.00279 * Float64(act)
    rss = f * cbrt(Float64(r_hel))
    κ_mp = R0_MP / rss
    tamp1 = 6.495 + 0.0229 * Float64(act)
    tamp2 = 1.6245 + 0.0088 * Float64(act)
    return MercuryKT17Field(internal, Float64(r_hel), rss, κ_mp, κ_mp^3, tamp1, tamp2)
end

@inline function (m::MercuryKT17Field)(r)
    inv_RM = inv(R_M)
    xm = r[1] * inv_RM
    ym = r[2] * inv_RM
    zm = r[3] * inv_RM

    ## Internal field from PlanetaryMagneticFields.jl (evaluated in planet-centered MSO)
    B_int = SVector{3}(
        m.internal(xm, ym, zm + Z_OFFSET; from = :cartesian, to = :cartesian)
    )

    ## KT17 external shielding and tail current fields (evaluated in scaled MSM)
    xs = m.κ_mp * xm
    ys = m.κ_mp * ym
    zs = m.κ_mp * zm

    B_dsk = kt17_taildisk(xs, ys, zs) + kt17_shield(C_DISKSHLD, P_DISKSHLD, xs, ys, zs)
    B_slb = kt17_tailslab(xs, ys, zs) + kt17_shield(C_SLABSHLD, P_SLABSHLD, xs, ys, zs)
    B_ext = m.κ_mp3 * kt17_shield(C_DIPSHLD, P_DIPSHLD, xs, ys, zs) +
        m.tamp1 * B_dsk + m.tamp2 * B_slb

    return (B_int + B_ext) * 1.0e-9 # [nT] -> [T]
end

@inline function is_inside_planet(r)
    inv_RM = inv(R_M)
    xm = r[1] * inv_RM
    ym = r[2] * inv_RM
    zm = r[3] * inv_RM + Z_OFFSET
    return xm^2 + ym^2 + zm^2 < 1.0
end

@inline function is_outside_mp(r, rss::Float64)
    inv_RM = inv(R_M)
    xm = r[1] * inv_RM
    ym = r[2] * inv_RM
    zm = r[3] * inv_RM
    rmag = sqrt(xm^2 + ym^2 + zm^2)
    rmag < 1.0e-12 && return false
    rmag > 4.0 && return true
    ct = xm / rmag
    1.0 + ct <= 1.0e-12 && return false
    r_mp = rss / sqrt(ALFA_MP * (1.0 + ct))
    return rmag > r_mp
end

# ### Dayside Magnetic Field Compression: Aphelion vs. Perihelion
#
# Let us visualize the equatorial magnetic field magnitude $|B|$ and trace a dayside field
# line anchored at $(X, Y, Z) = (1.2\,R_M, 0, 0)$ at aphelion ($r_{\mathrm{hel}} = 0.467\,\mathrm{au}$)
# and perihelion ($r_{\mathrm{hel}} = 0.307\,\mathrm{au}$). At perihelion, the enhanced
# magnetopause current produces a **local equatorial $|B|$ maximum** on the dayside, which
# splits the single equatorial magnetic mirror into two off-equatorial pockets.

B_aphelion = MercuryKT17Field(0.467)
B_perihelion = MercuryKT17Field(0.307)

function trace_dayside_fieldline(Bmodel::MercuryKT17Field, x0_RM::Float64)
    function rhs(u, p, s)
        B = Bmodel(u)
        return B / norm(B)
    end
    cb = DiscreteCallback(
        (u, s, integrator) -> is_inside_planet(u) || is_outside_mp(u, Bmodel.rss),
        terminate!;
        save_positions = (false, false),
    )
    u0 = SA[x0_RM * R_M, 0.0, 0.0]
    prob_fwd = ODEProblem(rhs, u0, (0.0, 4.0 * R_M))
    prob_bwd = ODEProblem((u, p, s) -> -rhs(u, p, s), u0, (0.0, 4.0 * R_M))
    sol_fwd = solve(
        prob_fwd, Vern7();
        callback = cb, dtmax = 0.05 * R_M, reltol = 1.0e-7, abstol = 1.0e-7,
    )
    sol_bwd = solve(
        prob_bwd, Vern7();
        callback = cb, dtmax = 0.05 * R_M, reltol = 1.0e-7, abstol = 1.0e-7,
    )

    s_vals = vcat(-reverse(sol_bwd.t), sol_fwd.t[2:end]) ./ R_M
    pts = vcat(reverse(sol_bwd.u), sol_fwd.u[2:end])
    B_nT = [norm(Bmodel(pt)) * 1.0e9 for pt in pts]
    return s_vals, B_nT, pts
end

f1 = Figure(size = (1260, 430), fontsize = 14)

xs_grid = range(-2.0, 2.0; length = 140)
ys_grid = range(-2.0, 2.0; length = 140)
θ_mp = range(-0.88π, 0.88π; length = 200)
θ_planet = range(0, 2π; length = 120)
r_eq_planet = sqrt(1.0 - Z_OFFSET^2) # planetary cross-section at Z_MSM = 0

rss_aph = round(B_aphelion.rss; digits = 2)
rss_per = round(B_perihelion.rss; digits = 2)
for (col, (Bmodel, title_str)) in enumerate(
        [
            (B_aphelion, "Aphelion (r_hel = 0.467 au, R_ss = $(rss_aph) R_M)"),
            (B_perihelion, "Perihelion (r_hel = 0.307 au, R_ss = $(rss_per) R_M)"),
        ]
    )
    ax = Axis(
        f1[1, col];
        title = title_str,
        xlabel = L"X\ [R_M]",
        ylabel = col == 1 ? L"Y\ [R_M]" : "",
        aspect = DataAspect(),
        xreversed = true,
        limits = ((-2.0, 2.0), (-2.0, 2.0)),
    )
    Bmag_eq = [
        begin
            r = SA[x * R_M, y * R_M, 0.0]
            if is_inside_planet(r) || is_outside_mp(r, Bmodel.rss)
                NaN
            else
                norm(Bmodel(r)) * 1.0e9
            end
        end for x in xs_grid, y in ys_grid
    ]
    hm = heatmap!(
        ax, xs_grid, ys_grid, Bmag_eq;
        colormap = :turbo, colorscale = log10,
        colorrange = (15.0, 300.0), nan_color = :white,
    )
    r_mp_curve = [Bmodel.rss / sqrt(ALFA_MP * (1.0 + cos(θ))) for θ in θ_mp]
    lines!(
        ax, r_mp_curve .* cos.(θ_mp), r_mp_curve .* sin.(θ_mp);
        color = :black, linewidth = 2, label = "Magnetopause",
    )
    poly!(
        ax, Point2f.(r_eq_planet .* cos.(θ_planet), r_eq_planet .* sin.(θ_planet));
        color = :gray35, strokecolor = :black, strokewidth = 1.5,
    )
    scatter!(
        ax, [1.2], [0.0];
        color = :white, strokecolor = :black, strokewidth = 1.5, markersize = 9,
    )
    if col == 2
        Colorbar(f1[1, 3], hm; label = L"|B|\ [\mathrm{nT}]")
    end
end

ax_prof = Axis(
    f1[1, 4];
    title = L"\mathrm{Dayside\ Field\ Line}\ (X_0 = 1.2\,R_M)",
    xlabel = L"\mathrm{Distance\ along\ field\ line}\ s\ [R_M]",
    ylabel = L"|B(s)|\ [\mathrm{nT}]",
    limits = ((-1.2, 1.2), (100.0, 320.0)),
)
s_aph, B_aph, _ = trace_dayside_fieldline(B_aphelion, 1.2)
s_per, B_per, _ = trace_dayside_fieldline(B_perihelion, 1.2)
lines!(
    ax_prof, s_aph, B_aph;
    color = :royalblue, linewidth = 2.5, label = "Aphelion (0.467 au)",
)
lines!(
    ax_prof, s_per, B_per;
    color = :crimson, linewidth = 2.5, label = "Perihelion (0.307 au)",
)
vlines!(ax_prof, [0.0]; color = :gray50, linestyle = :dash, linewidth = 1)
axislegend(ax_prof; position = :ct, framevisible = false)

colsize!(f1.layout, 1, Relative(0.28))
colsize!(f1.layout, 2, Relative(0.28))
colsize!(f1.layout, 4, Relative(0.36))

f1 = DisplayAs.PNG(f1) #hide

# ## Relativistic Guiding-Center Tracing: Conventional vs. Bifurcated Orbits
#
# Following Zhao et al. (2026, Eqs. 3–6), we trace energetic electrons using
# `TestParticle.jl`'s `RelativisticGuidingCenter` formulation (`TraceGCProblem(...,
# relativistic = true)`), whose state is $\mathbf{y} = (X_x, X_y, X_z, u_\parallel)$ with
# parallel proper velocity $u_\parallel = \gamma v_\parallel$ and conserved relativistic
# magnetic moment $\mu = p_\perp^2 / (2 m_e B)$.
#
# Each particle is initialized on the midnight magnetic equator ($X = -r_0, Y = Z = 0$) with
# kinetic energy $E_k$ and equatorial pitch angle $\alpha_0$. We monitor three physical loss
# channels via a termination callback:
#
# 1. **Planetary precipitation** onto Mercury's surface ($X^2 + Y^2 + (Z + 0.196\,R_M)^2 < R_M^2$),
#    which naturally captures the southern-hemisphere bounce-unstable zone.
# 2. **Magnetopause shadowing** outside the Shue magnetopause.
# 3. **Curvature-driven scattering** ([Büchner & Zelenyi, 1989](https://doi.org/10.1029/JA094iA09p11821))
#    when $\kappa_{\mathrm{BZ}} = \sqrt{R_c / r_g} < 3$, where $R_c = 1 / \|\boldsymbol{\kappa}\|$
#    and $r_g = \sqrt{2 \mu m_e B} / (|q_e| B)$.

@enum LossReason::UInt8 begin
    SURVIVED = 0
    PRECIPITATION = 1
    MAGNETOPAUSE = 2
    CURVATURE_SCATTERING = 3
end

@inline function check_loss(u, p, rss::Float64)
    r = SA[u[1], u[2], u[3]]
    is_inside_planet(r) && return PRECIPITATION
    is_outside_mp(r, rss) && return MAGNETOPAUSE
    ## Curvature scattering (κ < 3) occurs in the weak-field nightside current sheet;
    ## skip the expensive ForwardDiff Jacobian on the strong-field dayside/flanks.
    r[1] > 0.0 && return SURVIVED
    B = p.B(r, 0.0)
    any(isnan, B) && return MAGNETOPAUSE
    norm(B) * 1.0e9 > 90.0 && return SURVIVED
    mp = get_magnetic_properties(r, 0.0, p.B)
    Bmag = mp.Bmag
    (isnan(Bmag) || Bmag <= 0.0) && return CURVATURE_SCATTERING
    κ_mag = norm(mp.κ)
    κ_mag <= 0.0 && return SURVIVED
    rg = sqrt(2.0 * u[5] * Electron.m * Bmag) / (abs(p.q) * Bmag)
    ## κ_BZ = sqrt(R_c / r_g) < 3  <=>  r_g * |κ| > 1/9
    rg * κ_mag > 1.0 / 9.0 && return CURVATURE_SCATTERING
    return SURVIVED
end

"""
Precompute `Bmodel` on a 3D Cartesian grid and wrap it in a C² cubic B-spline
interpolator (`order = 3`) for fast, `ForwardDiff`-compatible ensemble tracing.
Points deep inside the planet (`r_MSO < 0.7 R_M`) are clamped to `r_MSO = 0.7 R_M`
so the central dipole singularity does not affect the spline outside Mercury's surface.
"""
function build_mercury_interpolator(Bmodel::MercuryKT17Field)
    xrange = range(-2.5 * R_M, 1.8 * R_M; length = 55)
    yrange = range(-2.2 * R_M, 2.2 * R_M; length = 55)
    zrange = range(-1.6 * R_M, 1.6 * R_M; length = 41)
    Bgrid = Array{Float64, 4}(undef, 3, length(xrange), length(yrange), length(zrange))
    r_core = 0.7 * R_M
    z_off_m = Z_OFFSET * R_M
    Threads.@threads for k in eachindex(zrange)
        z = zrange[k]
        z_mso = z + z_off_m
        @inbounds for j in eachindex(yrange), i in eachindex(xrange)
            x = xrange[i]
            y = yrange[j]
            r_mso = sqrt(x * x + y * y + z_mso * z_mso)
            if r_mso < r_core
                if r_mso > 0.0
                    s = r_core / r_mso
                    B = Bmodel(SA[x * s, y * s, z_mso * s - z_off_m])
                else
                    B = SA[0.0, 0.0, 0.0]
                end
            else
                B = Bmodel(SA[x, y, z])
            end
            Bgrid[1, i, j, k] = B[1]
            Bgrid[2, i, j, k] = B[2]
            Bgrid[3, i, j, k] = B[3]
        end
    end
    return build_interpolator(
        Bgrid, (xrange, yrange, zrange);
        order = 3, bc = TP.FillExtrap(NaN),
    )
end

function init_equatorial_electron(
        Bfield::F, r0_RM::Real, E_keV::Real, α_deg::Real;
        ϕ_rad::Real = Float64(π), sim_id::Int = 1,
    ) where {F}
    m = Electron.m
    q = Electron.q
    x0 = SA[r0_RM * cos(ϕ_rad) * R_M, r0_RM * sin(ϕ_rad) * R_M, 0.0]
    B_wrap = Bfield isa MercuryKT17Field ? TP.Field(Bfield) : Bfield
    B0 = norm(B_wrap(x0, 0.0))
    γ0 = 1.0 + (E_keV * 1.0e3 * eV) / (m * c^2)
    u_tot = c * sqrt(γ0^2 - 1.0)
    sα, cα = sincos(deg2rad(α_deg))
    u_par = u_tot * cα
    u_perp = u_tot * sα
    μ = 0.5 * m * u_perp^2 / B0
    u0_gc = SA[x0[1], x0[2], x0[3], u_par, μ]
    p_gc = (; q, q2m = q / m, c, E = TP.ZeroField(), B = B_wrap, sim_id)
    return u0_gc, p_gc
end

const GC_ABSTOL = SA[1.0e2, 1.0e2, 1.0e2, 1.0e3, Inf] # μ is constant

function trace_single_orbit(
        Bmodel::MercuryKT17Field;
        r0_RM::Float64 = 1.2, E_keV::Float64 = 300.0, α_deg::Float64 = 43.0,
        tspan::Tuple{Float64, Float64} = (0.0, 22.0),
    )
    u0_gc, p_gc = init_equatorial_electron(Bmodel, r0_RM, E_keV, α_deg)
    reason = Ref(SURVIVED)
    cb = DiscreteCallback(
        (u, t, integrator) -> begin
            r = check_loss(u, integrator.p, Bmodel.rss)
            if r != SURVIVED
                reason[] = r
                return true
            end
            return false
        end,
        terminate!;
        save_positions = (false, false),
    )
    prob = TraceGCProblem(u0_gc, tspan, p_gc; relativistic = true)
    sol = solve(
        prob, Vern7();
        callback = cb, reltol = 1.0e-6, abstol = GC_ABSTOL, maxiters = 200_000,
    )
    return sol, reason[]
end

## Compare a conventional orbit at aphelion vs. a bifurcated (Shabansky) orbit at perihelion
sol_conv, loss_conv = trace_single_orbit(
    B_aphelion; E_keV = 300.0, α_deg = 43.0, tspan = (0.0, 18.0)
)
sol_bifu, loss_bifu = trace_single_orbit(
    B_perihelion; E_keV = 300.0, α_deg = 65.0, tspan = (0.0, 18.0)
)

f2 = Figure(size = (1120, 520), fontsize = 14)

## Sphere mesh for Mercury (centered at (0, 0, -Z_OFFSET) in MSM)
u_sph = range(0, 2π; length = 48)
v_sph = range(0, π; length = 24)
xs_sph = [cos(u) * sin(v) for u in u_sph, v in v_sph]
ys_sph = [sin(u) * sin(v) for u in u_sph, v in v_sph]
zs_sph = [cos(v) - Z_OFFSET for u in u_sph, v in v_sph]

for (col, (sol, Bmodel, title_str, col_line)) in enumerate(
        [
            (
                sol_conv, B_aphelion,
                "Conventional Drift Orbit (Aphelion, α = 43°)", :royalblue,
            ),
            (
                sol_bifu, B_perihelion,
                "Bifurcated Shabansky Orbit (Perihelion, α = 65°)", :crimson,
            ),
        ]
    )
    ax = Axis3(
        f2[1, col];
        title = title_str,
        xlabel = L"X\ [R_M]",
        ylabel = L"Y\ [R_M]",
        zlabel = L"Z\ [R_M]",
        aspect = :data,
        limits = ((-1.8, 1.8), (-1.8, 1.8), (-1.4, 1.4)),
        elevation = 0.28,
        azimuth = 0.65,
    )
    surface!(
        ax, xs_sph, ys_sph, zs_sph;
        color = fill(0.45, size(xs_sph)), colormap = :grays, colorrange = (0, 1),
        shading = true,
    )

    ## Reference magnetic field lines on noon and midnight meridians
    for x_seed in (-1.2, 1.2)
        _, _, pts = trace_dayside_fieldline(Bmodel, x_seed)
        lines!(
            ax,
            [p[1] / R_M for p in pts],
            [p[2] / R_M for p in pts],
            [p[3] / R_M for p in pts];
            color = (:gray40, 0.55), linewidth = 1.5, linestyle = :dash,
        )
    end

    lines!(
        ax,
        [u[1] / R_M for u in sol.u],
        [u[2] / R_M for u in sol.u],
        [u[3] / R_M for u in sol.u];
        color = col_line, linewidth = 1.6,
    )
end

f2 = DisplayAs.PNG(f2) #hide

# At aphelion (left), the $300\,\mathrm{keV}$ electron bounces symmetrically across the
# magnetic equator between northern and southern mirror points and drifts around Mercury on
# a closed, adiabatic shell. At perihelion (right), once the electron drifts onto the
# compressed dayside, the equatorial $|B|$ maximum traps it in a single hemisphere
# (intra-hemispheric Shabansky bounce), breaking the second adiabatic invariant and driving
# it radially outward toward loss.

# ## Dynamical Regimes: Adiabatic, Diffusion, and Advection
#
# Following Fig. 4a–c of Zhao et al. (2026), we launch small ensembles of $300 \pm 60\,\mathrm{keV}$
# electrons from midnight ($r_0 = 1.2\,R_M$) and record their equatorial ($Z = 0$) crossing
# footprints over one drift cycle ($t \le 30\,\mathrm{s}$) across three representative regimes:
#
# 1. **Adiabatic** ($r_{\mathrm{hel}} = 0.467\,\mathrm{au}$, $\alpha_0 = 43^\circ$): narrow
#    closed ring with negligible radial spread.
# 2. **Diffusion-dominated** ($r_{\mathrm{hel}} = 0.387\,\mathrm{au}$, $\alpha_0 = 60^\circ$):
#    bifurcation crossings broaden the radial shell symmetrically.
# 3. **Advection-dominated** ($r_{\mathrm{hel}} = 0.387\,\mathrm{au}$, $\alpha_0 = 85^\circ$):
#    large-pitch-angle electrons experience strong outward radial advection into the nightside
#    tail and magnetopause.

function sample_equatorial_footprints(
        Bfield::F, rss::Float64, α_deg::Float64;
        n_particles::Int = 36, r0_RM::Float64 = 1.2, t_max::Float64 = 30.0, seed::Int = 42,
    ) where {F}
    rng = Xoshiro(seed)
    xs_all = Float64[]
    ys_all = Float64[]
    for _ in 1:n_particles
        E_keV = 240.0 + 120.0 * rand(rng)
        ϕ_rad = π + deg2rad(-15.0 + 30.0 * rand(rng))
        u0_gc, p_gc = init_equatorial_electron(Bfield, r0_RM, E_keV, α_deg; ϕ_rad)
        cb_loss = DiscreteCallback(
            (u, t, integrator) -> check_loss(u, integrator.p, rss) != SURVIVED,
            terminate!;
            save_positions = (false, false),
        )
        prob = TraceGCProblem(u0_gc, (0.0, t_max), p_gc; relativistic = true)
        sol = solve(
            prob, Vern6();
            callback = cb_loss, reltol = 1.0e-5, abstol = GC_ABSTOL, maxiters = 100_000,
        )
        ## Linearly interpolate Z = 0 crossings along the trajectory
        for k in 1:(length(sol.u) - 1)
            z1, z2 = sol.u[k][3], sol.u[k + 1][3]
            if (z1 <= 0.0 && z2 > 0.0) || (z1 >= 0.0 && z2 < 0.0)
                w = abs(z1) / (abs(z1) + abs(z2))
                push!(xs_all, ((1.0 - w) * sol.u[k][1] + w * sol.u[k + 1][1]) / R_M)
                push!(ys_all, ((1.0 - w) * sol.u[k][2] + w * sol.u[k + 1][2]) / R_M)
            end
        end
    end
    return xs_all, ys_all
end

## Precompute cubic B-spline interpolators for the ensemble runs (reused below)
B_mid = MercuryKT17Field(0.387)
B_itp_aph = build_mercury_interpolator(B_aphelion)
B_itp_mid = build_mercury_interpolator(B_mid)

regimes = [
    ("Adiabatic\n(0.467 au, α = 43°)", B_itp_aph, B_aphelion.rss, 43.0, :royalblue),
    ("Diffusion\n(0.387 au, α = 60°)", B_itp_mid, B_mid.rss, 60.0, :darkorange),
    ("Advection\n(0.387 au, α = 85°)", B_itp_mid, B_mid.rss, 85.0, :crimson),
]

f3 = Figure(size = (1180, 400), fontsize = 14)

for (col, (title_str, Bfield, rss, α_deg, col_pt)) in enumerate(regimes)
    xs_fp, ys_fp = sample_equatorial_footprints(Bfield, rss, α_deg; seed = 100 + col)
    ax = Axis(
        f3[1, col];
        title = title_str,
        xlabel = L"X\ [R_M]",
        ylabel = col == 1 ? L"Y\ [R_M]" : "",
        aspect = DataAspect(),
        xreversed = true,
        limits = ((-2.0, 2.0), (-2.0, 2.0)),
    )
    r_mp_curve = [rss / sqrt(ALFA_MP * (1.0 + cos(θ))) for θ in θ_mp]
    lines!(
        ax, r_mp_curve .* cos.(θ_mp), r_mp_curve .* sin.(θ_mp);
        color = :black, linewidth = 1.8,
    )
    poly!(
        ax, Point2f.(r_eq_planet .* cos.(θ_planet), r_eq_planet .* sin.(θ_planet));
        color = :gray35, strokecolor = :black, strokewidth = 1.5,
    )
    scatter!(ax, xs_fp, ys_fp; color = (col_pt, 0.35), markersize = 4)
end

f3 = DisplayAs.PNG(f3) #hide

# ## Radiation-Belt Survival vs. Solar Wind Driving and Pitch Angle
#
# Finally, we run multithreaded `EnsembleProblem` parameter sweeps (cf. Fig. 4e–f of Zhao et
# al., 2026) over $t \in [0, 240]\,\mathrm{s}$ using the precomputed cubic B-spline
# interpolators (`build_interpolator(..., order = 3)`) with `Vern6()` to quantify how the
# remnant (surviving) electron fraction depends on:
#
# - **Heliocentric distance** $r_{\mathrm{hel}} \in \{0.467, 0.427, 0.387, 0.347, 0.307\}\,\mathrm{au}$
#   ($E \in [240, 360]\,\mathrm{keV}$, $\alpha \in [65^\circ, 70^\circ]$, $r_0 = 1.2\,R_M$).
# - **Equatorial pitch angle** $\alpha$ in $8^\circ$ bins from $48^\circ$ to $88^\circ$ at
#   $r_{\mathrm{hel}} = 0.387\,\mathrm{au}$.

const T_END_ENS = 240.0
const N_ENS = 48

function run_survival_ensemble(
        Bfield::F, rss::Float64;
        α_range::Tuple{Float64, Float64} = (65.0, 70.0),
        E_range::Tuple{Float64, Float64} = (240.0, 360.0),
        r0_RM::Float64 = 1.2,
        n_particles::Int = N_ENS,
        seed::Int = 2026,
    ) where {F}
    u0_dummy, p_dummy = init_equatorial_electron(Bfield, r0_RM, 300.0, 67.5)
    base_prob = TraceGCProblem(u0_dummy, (0.0, T_END_ENS), p_dummy; relativistic = true)

    loss_times = fill(Inf, n_particles)
    loss_reasons = fill(SURVIVED, n_particles)

    function prob_func(prob, ctx)
        i = ctx.sim_id
        rng = Xoshiro(seed + 7919 * i)
        E_keV = E_range[1] + (E_range[2] - E_range[1]) * rand(rng)
        α_deg = α_range[1] + (α_range[2] - α_range[1]) * rand(rng)
        ϕ_rad = π + deg2rad(-20.0 + 40.0 * rand(rng))
        u0_i, p_i = init_equatorial_electron(Bfield, r0_RM, E_keV, α_deg; ϕ_rad, sim_id = i)
        return remake(prob; u0 = u0_i, p = p_i)
    end

    cb_loss = DiscreteCallback(
        (u, t, integrator) -> begin
            reason = check_loss(u, integrator.p, rss)
            if reason != SURVIVED
                idx = integrator.p.sim_id
                loss_times[idx] = t
                loss_reasons[idx] = reason
                return true
            end
            return false
        end,
        terminate!;
        save_positions = (false, false),
    )

    ens_prob = EnsembleProblem(base_prob; prob_func, safetycopy = false)
    solve(
        ens_prob, Vern6(), EnsembleThreads();
        trajectories = n_particles,
        callback = cb_loss,
        save_everystep = false, save_start = false,
        reltol = 1.0e-5, abstol = GC_ABSTOL, maxiters = 250_000,
    )
    return loss_times, loss_reasons
end

t_eval = range(0.0, T_END_ENS; length = 121)
function remnant_curve(loss_times)
    n = length(loss_times)
    return [100.0 * count(>(t), loss_times) / n for t in t_eval]
end

B_models_rhel = Dict(
    0.467 => (B_itp_aph, B_aphelion.rss),
    0.427 => let m = MercuryKT17Field(0.427)
        (build_mercury_interpolator(m), m.rss)
    end,
    0.387 => (B_itp_mid, B_mid.rss),
    0.347 => let m = MercuryKT17Field(0.347)
        (build_mercury_interpolator(m), m.rss)
    end,
    0.307 => (build_mercury_interpolator(B_perihelion), B_perihelion.rss),
)

r_hel_list = [0.467, 0.427, 0.387, 0.347, 0.307]
rhel_results = [
    let (Bitp, rss) = B_models_rhel[r_hel]
        (r_hel, run_survival_ensemble(Bitp, rss; seed = 1000 + k)...)
    end
        for (k, r_hel) in enumerate(r_hel_list)
]

α_bins = [(48.0, 56.0), (56.0, 64.0), (64.0, 72.0), (72.0, 80.0), (80.0, 88.0)]
α_results = [
    (α_bin, run_survival_ensemble(B_itp_mid, B_mid.rss; α_range = α_bin, seed = 2000 + k)...)
        for (k, α_bin) in enumerate(α_bins)
]

f4 = Figure(size = (1120, 500), fontsize = 14)
colors_5 = [:navy, :dodgerblue3, :mediumseagreen, :darkorange, :crimson]

ax_rhel = Axis(
    f4[1, 1];
    title = "Heliocentric Distance Dependence (α = 65°–70°)",
    xlabel = "Time [s]",
    ylabel = "Remnant fraction [%]",
    limits = ((0, T_END_ENS), (0, 104)),
)
for ((r_hel, lt, _), col) in zip(rhel_results, colors_5)
    lines!(
        ax_rhel, t_eval, remnant_curve(lt);
        color = col, linewidth = 2.4, label = "$(r_hel) au",
    )
end
Legend(f4[2, 1], ax_rhel; nbanks = 3, orientation = :horizontal, framevisible = false)

ax_alpha = Axis(
    f4[1, 2];
    title = "Equatorial Pitch Angle Dependence (r_hel = 0.387 au)",
    xlabel = "Time [s]",
    ylabel = "Remnant fraction [%]",
    limits = ((0, T_END_ENS), (0, 104)),
)
for (((a1, a2), lt, _), col) in zip(α_results, colors_5)
    lines!(
        ax_alpha, t_eval, remnant_curve(lt);
        color = col, linewidth = 2.4, label = "$(Int(a1))°–$(Int(a2))°",
    )
end
Legend(f4[2, 2], ax_alpha; nbanks = 3, orientation = :horizontal, framevisible = false)

f4 = DisplayAs.PNG(f4) #hide

# ### Loss-Channel Breakdown vs. Heliocentric Distance
#
# As Mercury moves from aphelion ($0.467\,\mathrm{au}$) to perihelion ($0.307\,\mathrm{au}$),
# solar-wind compression drives rapid drift-orbit bifurcation, shifting particles from
# long-term adiabatic trapping into magnetopause shadowing and nightside current-sheet
# curvature scattering ($\kappa < 3$):

io = IOBuffer() #hide
println( #hide
    io, #hide
    "| `r_hel` [au] | `R_ss` [`R_M`] | Survived [%] | Precipitation [%] | Magnetopause [%] | Curvature (κ < 3) [%] |", #hide
) #hide
println(io, "| :---: | :---: | :---: | :---: | :---: | :---: |") #hide
for (r_hel, _, reasons) in rhel_results #hide
    n = length(reasons) #hide
    rss = MercuryKT17Field(r_hel).rss #hide
    surv = 100.0 * count(==(SURVIVED), reasons) / n #hide
    prec = 100.0 * count(==(PRECIPITATION), reasons) / n #hide
    mp = 100.0 * count(==(MAGNETOPAUSE), reasons) / n #hide
    curv = 100.0 * count(==(CURVATURE_SCATTERING), reasons) / n #hide
    @printf( #hide
        io, "| %.3f | %.2f | %.1f | %.1f | %.1f | %.1f |\n", #hide
        r_hel, rss, surv, prec, mp, curv, #hide
    ) #hide
end #hide
Markdown.parse(String(take!(io))) #hide
