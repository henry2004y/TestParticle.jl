#ifndef TESTPARTICLE_H
#define TESTPARTICLE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================= */
/* Single-Particle Pushers (Node-Synchronized & Staggered)                    */
/* ========================================================================= */

/**
 * Single-particle Boris push for node-synchronized velocities (Float64).
 * Advances position and velocity from node t to node t + dt in-place.
 */
void testparticle_boris_push_f64(
    double* x,
    double* v,
    double q2m,
    double dt,
    const double* E,
    const double* B
);

/**
 * Single-particle Boris push for node-synchronized velocities (Float32).
 */
void testparticle_boris_push_f32(
    float* x,
    float* v,
    float q2m,
    float dt,
    const float* E,
    const float* B
);

/**
 * Single-particle Boris push for staggered velocities (Float64).
 * Advances r from t to t + dt, and v_half from t - dt/2 to t + dt/2.
 */
void testparticle_boris_step_staggered_f64(
    double* x,
    double* v_half,
    double q2m,
    double dt,
    const double* E,
    const double* B
);

/**
 * Single-particle Boris push for staggered velocities (Float32).
 */
void testparticle_boris_step_staggered_f32(
    float* x,
    float* v_half,
    float q2m,
    float dt,
    const float* E,
    const float* B
);

/**
 * Pure Boris velocity rotation kernel (Float64).
 * Updates velocity vector v in-place using electric and magnetic fields.
 */
void testparticle_update_velocity_boris_f64(
    double* v,
    const double* E,
    const double* B,
    double qdt2m
);

/**
 * Pure Boris velocity rotation kernel (Float32).
 */
void testparticle_update_velocity_boris_f32(
    float* v,
    const float* E,
    const float* B,
    float qdt2m
);

/* ========================================================================= */
/* High-Order & Multistep Hyper-Boris Pushers (Zenitani & Kato 2025)         */
/* ========================================================================= */

/**
 * Multistep Boris push with n subcycles and order correction N in {2, 4, 6}.
 * Synchronized (Float64).
 */
void testparticle_multistep_boris_push_f64(
    double* x,
    double* v,
    double q2m,
    double dt,
    const double* E,
    const double* B,
    int32_t n,
    int32_t order
);

/**
 * Multistep Boris push with n subcycles and order correction N in {2, 4, 6}.
 * Synchronized (Float32).
 */
void testparticle_multistep_boris_push_f32(
    float* x,
    float* v,
    float q2m,
    float dt,
    const float* E,
    const float* B,
    int32_t n,
    int32_t order
);

/**
 * Multistep Boris push with n subcycles and order correction N in {2, 4, 6}.
 * Staggered (Float64).
 */
void testparticle_multistep_boris_step_staggered_f64(
    double* x,
    double* v_half,
    double q2m,
    double dt,
    const double* E,
    const double* B,
    int32_t n,
    int32_t order
);

/**
 * Multistep Boris push with n subcycles and order correction N in {2, 4, 6}.
 * Staggered (Float32).
 */
void testparticle_multistep_boris_step_staggered_f32(
    float* x,
    float* v_half,
    float q2m,
    float dt,
    const float* E,
    const float* B,
    int32_t n,
    int32_t order
);

/**
 * Hyper-Boris 4th-order gyrophase-corrected pusher (Float64).
 */
void testparticle_hyper_boris4_push_f64(
    double* x,
    double* v,
    double q2m,
    double dt,
    const double* E,
    const double* B,
    int32_t n
);

/**
 * Hyper-Boris 4th-order gyrophase-corrected pusher (Float32).
 */
void testparticle_hyper_boris4_push_f32(
    float* x,
    float* v,
    float q2m,
    float dt,
    const float* E,
    const float* B,
    int32_t n
);

/**
 * Hyper-Boris 6th-order gyrophase-corrected pusher (Float64).
 */
void testparticle_hyper_boris6_push_f64(
    double* x,
    double* v,
    double q2m,
    double dt,
    const double* E,
    const double* B,
    int32_t n
);

/**
 * Hyper-Boris 6th-order gyrophase-corrected pusher (Float32).
 */
void testparticle_hyper_boris6_push_f32(
    float* x,
    float* v,
    float q2m,
    float dt,
    const float* E,
    const float* B,
    int32_t n
);

/* ========================================================================= */
/* Adaptive Boris & Velocity Resynchronization                               */
/* ========================================================================= */

/**
 * Propose gyroperiod-tracking time step size dt = safety * 2π / (|q2m| * |B|).
 */
double testparticle_boris_propose_dt_f64(
    double q2m,
    const double* B,
    double safety
);

float testparticle_boris_propose_dt_f32(
    float q2m,
    const float* B,
    float safety
);

/**
 * Resynchronize staggered velocity when time step changes from dt_prev to dt_next.
 * Preserves time-reversibility and 2nd-order accuracy under variable dt.
 */
void testparticle_boris_resync_velocity_half_f64(
    double* v_half,
    double q2m,
    double dt_prev,
    double dt_next,
    const double* E,
    const double* B
);

void testparticle_boris_resync_velocity_half_f32(
    float* v_half,
    float q2m,
    float dt_prev,
    float dt_next,
    const float* E,
    const float* B
);

/**
 * Advance an adaptive Boris step and propose / resynchronize for the next step.
 */
void testparticle_adaptive_boris_step_f64(
    double* x,
    double* v_half,
    double dt_current,
    double* dt_next,
    double q2m,
    const double* E,
    const double* B,
    double safety
);

void testparticle_adaptive_boris_step_f32(
    float* x,
    float* v_half,
    float dt_current,
    float* dt_next,
    float q2m,
    const float* E,
    const float* B,
    float safety
);

/* ========================================================================= */
/* Batch & Ensemble Pushers (AoS, SoA, and Multi-Species)                    */
/* ========================================================================= */

/**
 * Batch Boris push for AoS layout (Float64).
 * Memory layout: x[3*n], v[3*n], E[3*n], B[3*n].
 */
void testparticle_boris_push_batch_f64(
    size_t n,
    double* x,
    double* v,
    double q2m,
    double dt,
    const double* E,
    const double* B
);

void testparticle_boris_push_batch_f32(
    size_t n,
    float* x,
    float* v,
    float q2m,
    float dt,
    const float* E,
    const float* B
);

/**
 * Batch Boris push with n_sub subcycles and order correction (Float64).
 */
void testparticle_multistep_boris_push_batch_f64(
    size_t n_particles,
    double* x,
    double* v,
    double q2m,
    double dt,
    const double* E,
    const double* B,
    int32_t n_sub,
    int32_t order
);

void testparticle_multistep_boris_push_batch_f32(
    size_t n_particles,
    float* x,
    float* v,
    float q2m,
    float dt,
    const float* E,
    const float* B,
    int32_t n_sub,
    int32_t order
);

/**
 * Batch Boris push for Structure-of-Arrays (SoA) layout (Float64).
 * Continuous component arrays: x[n], y[n], z[n], vx[n], vy[n], vz[n].
 * If q2m_arr != NULL, uses per-particle q2m; otherwise uses q2m_scalar.
 */
void testparticle_boris_push_batch_soa_f64(
    size_t n,
    double* x, double* y, double* z,
    double* vx, double* vy, double* vz,
    const double* q2m_arr,
    double q2m_scalar,
    double dt,
    const double* Ex, const double* Ey, const double* Ez,
    const double* Bx, const double* By, const double* Bz
);

void testparticle_boris_push_batch_soa_f32(
    size_t n,
    float* x, float* y, float* z,
    float* vx, float* vy, float* vz,
    const float* q2m_arr,
    float q2m_scalar,
    float dt,
    const float* Ex, const float* Ey, const float* Ez,
    const float* Bx, const float* By, const float* Bz
);

/**
 * Staggered leapfrog batch push for Structure-of-Arrays (SoA) layout (Float64).
 */
void testparticle_boris_step_staggered_batch_soa_f64(
    size_t n,
    double* x, double* y, double* z,
    double* vhx, double* vhy, double* vhz,
    const double* q2m_arr,
    double q2m_scalar,
    double dt,
    const double* Ex, const double* Ey, const double* Ez,
    const double* Bx, const double* By, const double* Bz
);

/**
 * Multi-species AoS batch push with individual q2m array (Float64).
 */
void testparticle_boris_push_batch_multispecies_f64(
    size_t n,
    double* x,
    double* v,
    const double* q2m_arr,
    double dt,
    const double* E,
    const double* B
);

/* ========================================================================= */
/* Trajectory Tracing & Dynamic Field Callbacks                              */
/* ========================================================================= */

/**
 * Callback function to evaluate electric and magnetic fields at position x and time t.
 */
typedef void (*testparticle_field_callback_f64)(
    double t,
    const double* x,
    double* E,
    double* B,
    void* user_data
);

typedef void (*testparticle_field_callback_f32)(
    float t,
    const float* x,
    float* E,
    float* B,
    void* user_data
);

/**
 * Integrate a particle orbit over [t_start, t_end] with dynamic field evaluations.
 * If trajectory_out != NULL, saves interleaved (t, x, y, z, vx, vy, vz) per step.
 * Returns 0 on success, negative value on error.
 */
int32_t testparticle_trace_boris_f64(
    double* x,
    double* v,
    double t_start,
    double t_end,
    double dt,
    double q2m,
    testparticle_field_callback_f64 field_fn,
    void* user_data,
    size_t max_steps,
    size_t* n_steps_taken,
    double* trajectory_out
);

int32_t testparticle_trace_boris_f32(
    float* x,
    float* v,
    float t_start,
    float t_end,
    float dt,
    float q2m,
    testparticle_field_callback_f32 field_fn,
    void* user_data,
    size_t max_steps,
    size_t* n_steps_taken,
    float* trajectory_out
);

#ifdef __cplusplus
}
#endif

#endif /* TESTPARTICLE_H */
