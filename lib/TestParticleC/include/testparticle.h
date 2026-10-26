#ifndef TESTPARTICLE_H
#define TESTPARTICLE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
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

/*
 * Single-particle Boris push for node-synchronized velocities (Float32).
 * Advances position and velocity from node t to node t + dt in-place.
 */
void testparticle_boris_push_f32(
    float* x,
    float* v,
    float q2m,
    float dt,
    const float* E,
    const float* B
);

/*
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

/*
 * Single-particle Boris push for staggered velocities (Float32).
 * Advances r from t to t + dt, and v_half from t - dt/2 to t + dt/2.
 */
void testparticle_boris_step_staggered_f32(
    float* x,
    float* v_half,
    float q2m,
    float dt,
    const float* E,
    const float* B
);

/*
 * Vectorized Boris push for an ensemble/batch of n particles (Float64).
 * Array layouts: x[3*n], v[3*n], E[3*n], B[3*n] where component k of particle i
 * is at index 3 * i + k.
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

/*
 * Vectorized Boris push for an ensemble/batch of n particles (Float32).
 */
void testparticle_boris_push_batch_f32(
    size_t n,
    float* x,
    float* v,
    float q2m,
    float dt,
    const float* E,
    const float* B
);

/*
 * Pure Boris velocity rotation kernel (Float64).
 * Updates velocity vector v in-place using electric and magnetic fields.
 */
void testparticle_update_velocity_boris_f64(
    double* v,
    const double* E,
    const double* B,
    double qdt2m
);

/*
 * Pure Boris velocity rotation kernel (Float32).
 */
void testparticle_update_velocity_boris_f32(
    float* v,
    const float* E,
    const float* B,
    float qdt2m
);

#ifdef __cplusplus
}
#endif

#endif /* TESTPARTICLE_H */
