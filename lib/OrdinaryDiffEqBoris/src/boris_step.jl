# The step API.
#
# A Boris step is a map from one staggered state to the next, so it can be
# written without a cache: everything it needs is passed in and everything it
# produces is returned. `perform_step!` and a GPU kernel can then share one
# implementation, which a cache cannot give them, since a kernel cannot capture
# a mutable struct.
#
# The state is the pair `(r, v_half)`: the position at a node and the velocity
# half a step behind it. Keeping the two apart is what makes a step cheap.
# Advancing costs a single field evaluation, and the node velocity, which is
# what one reports, is reconstructed only when something asks for it.

"""
    advance_boris(v_half, r, dt, t, p, alg) -> (r_new, v_half_new)

Advance one step of size `dt` from the node time `t`.

`r` is the position at `t` and `v_half` the velocity half a step behind it, at
`t - dt/2`. The fields are evaluated once, at `t + dt/2`, and the pair returned
is the position at `t + dt` with the velocity half a step ahead of it, which is
the same staggering as the input and can be fed straight into the next step.
"""
@inline @muladd function advance_boris(v_half, r, dt, t, p, alg)
    v_half_new = update_velocity(v_half, r, dt, t + 0.5 * dt, p, alg)
    r_new = r + v_half_new * dt

    return r_new, v_half_new
end

"""
    update_velocity_node(v_half, r, dt, t, p, alg) -> v_node

Reconstruct the velocity at the node `t` from the velocity half a step behind
it, pushing it forward by `dt/2` with the fields evaluated at `(r, t)`.

This is the velocity to report, and the one a saved state holds. It costs a
field evaluation a step would not otherwise make, so a caller that saves rarely
should call it rarely.
"""
@inline function update_velocity_node(v_half, r, dt, t, p, alg)
    return update_velocity(v_half, r, 0.5 * dt, t, p, alg)
end

"""
    update_velocity_half(v, r, dt, t, p, alg) -> v_half

Move the velocity at the node `t` back by `dt/2`, giving the half-step velocity
the integration carries. This is how a trajectory is started, from an initial
condition given at a node.
"""
@inline function update_velocity_half(v, r, dt, t, p, alg)
    return update_velocity(v, r, -0.5 * dt, t, p, alg)
end

"""
    update_velocity_resync(v_half, r, dt_prev, dt, t, p, alg) -> v_half_new

Move a half-step velocity centred on `dt_prev` onto the half step of `dt`, both
at the node `t`, by going through the node velocity.

Changing the step size re-centres the velocity, which is what keeps the scheme
time-reversible under adaptive stepping: without it a change of step would leave
the velocity staggered against the position.
"""
@inline function update_velocity_resync(v_half, r, dt_prev, dt, t, p, alg)
    v_node = update_velocity_node(v_half, r, dt_prev, t, p, alg)

    return update_velocity_half(v_node, r, dt, t, p, alg)
end
