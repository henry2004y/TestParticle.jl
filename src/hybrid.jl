# Adaptive hybrid GC-FO solver

@noinline function _verner_missing()
    throw(
        ArgumentError(
            "The default guiding-center solver of AdaptiveHybrid requires " *
            "OrdinaryDiffEqVerner. Run `using OrdinaryDiffEqVerner` first, or pass " *
            "`gc_alg` explicitly to AdaptiveHybrid."
        )
    )
end

default_gc_alg(::Any) = _verner_missing()
default_gc_alg() = default_gc_alg(Val(:default))

"""
    AdaptiveHybrid{T, A}

Hybrid algorithm that switches between adaptive Boris (FO) and a guiding-center ODE
solver `gc_alg` (defaulting to `Vern7()` when `OrdinaryDiffEqVerner` is loaded), used as
`solve(TraceProblem(u0, tspan, p), alg)` with the full-orbit state `u0`.

Mode switching uses hysteresis with two adiabaticity thresholds to avoid
rapid GC ↔ FO chattering near a single bound:
- `threshold_gc_to_fo` (α): switch GC → FO when `ε ≥ α`.
- `threshold_fo_to_gc` (β): switch FO → GC when `ε < β`.
With `α ≥ β` the band `β ≤ ε < α` forms a buffer where the current mode is
retained, so a particle hovering near the boundary does not oscillate.

The adiabaticity parameter `ε` compared against the thresholds is selected by
the `adiabaticity` field, which may be:
- `:curvature` (default) → `ε_curv = ρ_L / R_c` (curvature drift; reproduces
  the legacy `get_adiabaticity` behavior exactly),
- `:gradB` → `ε_gradB = ρ_L / L_B`, with `L_B = |B| / |∇B|`,
- `:both` → OR of the two criteria: the solver switches to the full orbit
  whenever *either* `ε_curv ≥ α` *or* `ε_gradB ≥ α` (equivalently
  `max(ε_curv, ε_gradB) ≥ α`),
- `:jacobian` → `ε_jac = ρ_L · ‖JB‖_F / |B|`, CHIMP's original all-in-one
  criterion (Bates & Fuentes, 2020), where `JB` is the Jacobian of `B` and
  `‖JB‖_F` its Frobenius norm. It captures curvature, grad-B, torsion and
  shear at once, so it switches to full orbit in non-adiabatic regions that
  `:curvature`/`:gradB` miss (e.g. a field whose direction rotates in space
  while keeping constant magnitude).

Both modes adapt their step size to the local field. The full-orbit mode uses a
Boris step `dt = 2π * safety_fo / (|q/m| * |B|)` clipped to `[dtmin, dtmax]` and
recomputed after every accepted step; when the step changes, the staggered
velocity is re-centered on the new step so that the leapfrog scheme stays
time-reversible. The guiding-center mode steps with `gc_alg`, so its steps may
span many gyroperiods.

Because a guiding-center step can be long compared with the scale on which `ε`
varies, the adiabaticity check is triggered by whichever comes first:
`check_interval` accepted steps, or `check_dt` in simulation time since the
previous check. A step is also clipped so that it never crosses the next check
deadline, which makes `check_dt` a hard bound on the simulated time a particle
can travel between examinations rather than a bound that a long step could
overshoot. `check_dt` defaults to `100 * dtmax`: small enough to catch a
transition that a long guiding-center step would otherwise jump over, and large
enough that the clamp does not turn `check_dt` into a step-size cap that costs
more than the monitoring saves. Lower it where `ε` changes on a shorter scale.

While the solver is in the full-orbit mode, `check_margin_max > 1` lengthens
both bounds by up to that factor in proportion to how far `ε` sits above
`threshold_fo_to_gc`: `ε` changes by `O(ε)` over one curvature radius, so an `ε`
that is a factor `k` above the exit threshold needs of order `k` such excursions
before the mode can flip, and checking in between only costs time. The interval
collapses back to its nominal value as `ε` approaches the threshold, and the
guiding-center mode always uses the nominal interval, so a transition *into* a
non-adiabatic region is detected as early as with a uniform interval.

It is off by default (`check_margin_max = 1`). Stretching reorders the mode
history enough to matter: on the magnetic-bottle benchmark `4` makes the run
about 40% faster but the RMS position error 50% larger, because the full-orbit
excursions become fewer and longer and the full-orbit step is the coarser one.
Raise it only where the full-orbit mode is cheap relative to the guiding-center
mode.

On the switch to the guiding center, the state is refined by up to
`switch_correction` fixed-point iterations of `X = x - ρ(X)` instead of the
zeroth-order `X = x - ρ(x)` map, which removes the `O(ε)` jump that the switch
would otherwise introduce. That map is the exact inverse of the reconstruction
back to the full orbit, so the round trip through a switch returns the particle
where it was. Set `switch_correction = 0` for the legacy behavior.

Each check point records only the adiabaticity value used by the selected
`adiabaticity` mode (a single scalar per check) together with the active mode
in `sol.stats.adiabaticity`. This diagnostic is only saved when
`save_adiabaticity = true` (the default); set it to `false` to skip the
buffers and reduce memory/allocation overhead.
"""
struct AdaptiveHybrid{T, A}
    threshold_gc_to_fo::T
    threshold_fo_to_gc::T
    dtmin::T
    dtmax::T
    safety_fo::T
    abstol::T
    reltol::T
    check_dt::T
    check_margin_max::T
    maxiters::Int
    check_interval::Int
    adiabaticity::Symbol
    save_adiabaticity::Bool
    switch_correction::Int
    gc_alg::A
end

function AdaptiveHybrid(;
        threshold = 0.1,
        threshold_gc_to_fo = threshold,
        threshold_fo_to_gc = threshold,
        dtmax, dtmin = 1.0e-2 * dtmax, check_dt = 100 * dtmax,
        check_margin_max = 1.0, safety_fo = 0.1,
        abstol = 1.0e-6, reltol = 1.0e-6, maxiters = 10000, check_interval = 10,
        adiabaticity = :curvature, save_adiabaticity = true, switch_correction = 2,
        gc_alg = nothing
    )
    dtmax > 0 || throw(ArgumentError("dtmax must be positive."))
    dtmin > 0 || throw(ArgumentError("dtmin must be positive."))
    dtmin <= dtmax || throw(ArgumentError("dtmin must be <= dtmax."))
    check_dt > 0 || throw(ArgumentError("check_dt must be positive."))
    check_margin_max >= 1 || throw(ArgumentError("check_margin_max must be >= 1."))
    safety_fo > 0 || throw(ArgumentError("safety_fo must be positive."))
    maxiters > 0 || throw(ArgumentError("maxiters must be positive."))
    threshold_fo_to_gc >= 0 ||
        throw(ArgumentError("threshold_fo_to_gc must be non-negative."))
    abstol > 0 || throw(ArgumentError("abstol must be positive."))
    reltol > 0 || throw(ArgumentError("reltol must be positive."))
    check_interval > 0 || throw(ArgumentError("check_interval must be positive."))
    switch_correction >= 0 ||
        throw(ArgumentError("switch_correction must be non-negative."))
    threshold_gc_to_fo >= threshold_fo_to_gc ||
        throw(
        ArgumentError(
            "threshold_gc_to_fo must be >= threshold_fo_to_gc for hysteresis."
        )
    )
    adiabaticity in (:curvature, :gradB, :both, :jacobian) ||
        throw(
        ArgumentError(
            "adiabaticity must be one of :curvature, :gradB, :both, :jacobian."
        )
    )
    T = promote_type(
        typeof(threshold_gc_to_fo), typeof(threshold_fo_to_gc),
        typeof(dtmin), typeof(dtmax), typeof(check_dt), typeof(check_margin_max),
        typeof(safety_fo), typeof(abstol), typeof(reltol)
    )
    return AdaptiveHybrid{T, typeof(gc_alg)}(
        T(threshold_gc_to_fo), T(threshold_fo_to_gc), T(dtmin), T(dtmax),
        T(safety_fo), T(abstol), T(reltol), T(check_dt), T(check_margin_max),
        maxiters, check_interval, adiabaticity, save_adiabaticity,
        switch_correction, gc_alg
    )
end

function init_gc_integrator(
        xv::SVector{4, T}, tspan, p, dt, alg::AdaptiveHybrid
    ) where {T}
    gc_alg = alg.gc_alg === nothing ? default_gc_alg() : alg.gc_alg
    prob = _gc_ode_problem(xv, (T(tspan[1]), T(tspan[2])), p)
    return SciMLBase.__init(
        prob, gc_alg;
        dt = T(dt), abstol = alg.abstol, reltol = alg.reltol, maxiters = alg.maxiters,
        save_start = false, save_end = false, save_everystep = false, dense = false
    )
end

function reinit_gc_integrator!(
        integrator::SciMLBase.DEIntegrator, xv::SVector{4, T}, tspan, p, dt
    ) where {T}
    t0 = typeof(integrator.t)(tspan[1])
    tf = typeof(integrator.t)(tspan[2])
    dt_t = typeof(integrator.t)(dt)
    integrator.p = p
    integrator.dt = dt_t
    integrator.opts.dtmax = tf - t0
    integrator.force_stepfail = false
    integrator.last_stepfail = false
    integrator.accept_step = false
    integrator.isout = false
    integrator.reeval_fsal = false
    SciMLBase.set_proposed_dt!(integrator, dt_t)
    SciMLBase.reinit!(
        integrator, xv;
        t0, tf,
        erase_sol = false, reset_dt = false,
        reinit_callbacks = false, reinit_retcode = false, reinit_cache = false
    )
    if integrator.fsalfirst !== nothing
        integrator.fsalfirst = integrator.f(xv, p, t0)
    end
    return integrator
end

function gc_step!(
        integrator::SciMLBase.DEIntegrator, t_limit, isoutside::F, p, tspan_end
    ) where {F}
    t_before = integrator.t
    dt_max = t_limit - t_before
    dt_next = integrator.iter == 0 ? integrator.dt : SciMLBase.get_proposed_dt(integrator)
    if zero(dt_max) < dt_max < dt_next && t_limit < tspan_end
        if integrator.iter == 0
            integrator.dt = dt_max
        else
            SciMLBase.set_proposed_dt!(integrator, dt_max)
        end
    end
    SciMLBase.step!(integrator)

    y = SVector{4}(integrator.u)
    isoutside(y, p, integrator.t) && return ReturnCode.Terminated

    code = SciMLBase.check_error(integrator)
    code != ReturnCode.Success && return code
    return integrator.t == t_before ? ReturnCode.DtLessThanMin : ReturnCode.Success
end
