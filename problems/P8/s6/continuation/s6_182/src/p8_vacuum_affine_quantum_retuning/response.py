"""Actual new stationary response, not an adaptive cancellation of its loop."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_response import quantum as previous
from p8_vacuum_affine_kernel_norm import estimate as kernel_norm

from . import bounds, profile, ward

u = profile.u
rho, P, H = profile.rho, profile.P, ward.H
S = rho + P
eta, w = (s.Function(name)(u) for name in ("prepared_time_shift", "invariant_scale"))
R, Q = (
    s.Function(name)(u)
    for name in ("actual_pure_scale_rho_response", "actual_pure_scale_P_response")
)


def simplify_response(value):
    value = s.expand(value.doit()).subs(
        s.diff(R, u), -3 * H * (R + Q) - 3 * S * s.diff(w, u)
    )
    return ward.reduce_ward(value)


@cache
def data():
    n = s.diff(eta, u)
    v = w + H * eta
    drho = R + s.diff(rho, u) * eta
    dP = Q + s.diff(P, u) * eta
    phys = s.Matrix([-drho - S * n, 3 * dP + 3 * S * n])
    adapted = s.Matrix([-ward.D(phys[0]) + H * phys[1], phys[1]])
    eta_local = (
        S * s.diff(eta, u, 2)
        + s.diff(P, u) * n
        - 3 * s.diff(H, u) * S * eta
        - 3 * S * s.diff(w, u)
    )
    scale_local = 3 * s.diff(P, u) * eta + 3 * S * n
    full = s.Matrix([eta_local, 3 * Q + scale_local])
    density = (
        -S * n * n / 2
        - s.Rational(3, 2) * s.diff(H, u) * S * eta * eta
        - 3 * S * eta * s.diff(w, u)
    )
    euler = lambda f: s.diff(density, f) - ward.D(s.diff(density, s.diff(f, u)))
    old_map = {
        previous.H: H,
        previous.rho: rho,
        previous.P: P,
        previous.eta: eta,
        previous.w: w,
        previous.R: R,
        previous.Q: Q,
    }
    # Rename the derivative variable FIRST. Simultaneously inserting
    # functions of u into derivatives with respect to the old time
    # would spuriously differentiate those functions as constants.
    old_map = {key.subs(previous.u, u): value for key, value in old_map.items()}
    old = (
        previous.data()["complete_adapted_vector_response"]
        .subs(previous.u, u)
        .subs(old_map, simultaneous=True)
        .doit()
    )
    profile_phys = ward.data()["profile_current_Hessian"] * s.Matrix([n, v])
    profile_adapted = s.Matrix(
        [-ward.D(profile_phys[0]) + H * profile_phys[1], profile_phys[1]]
    )
    # Literal fixed-profile L2 and its Euler currents give a second route.
    added_density = (
        -(rho + P) * n * n / 2 + 3 * rho * n * v - s.Rational(9, 2) * P * v * v
    )
    direct_added = s.Matrix(
        [
            s.diff(added_density, f) - ward.D(s.diff(added_density, s.diff(f, u)))
            for f in (eta, w)
        ]
    )
    checks = {
        "new_time_channel_from_actual_stationary_current": simplify_response(
            adapted[0] - full[0]
        ),
        "new_scale_channel_from_actual_stationary_current": s.expand(
            adapted[1] - full[1]
        ),
        "new_local_time_density": ward.reduce_ward(euler(eta) - eta_local),
        "new_local_scale_density": ward.reduce_ward(euler(w) - scale_local),
        "old_uncancelled_current_plus_literal_new_profile": s.ImmutableMatrix(
            (old + profile_adapted - full).applyfunc(simplify_response)
        ),
        "profile_from_literal_quadratic_action": s.ImmutableMatrix(
            (direct_added - profile_adapted).applyfunc(ward.reduce_ward)
        ),
        "new_second_time_coefficient_retained": s.diff(eta_local, s.diff(eta, u, 2))
        - S,
        "new_first_time_coefficient_uses_Ward": ward.reduce_ward(
            s.diff(P, u) - (s.diff(S, u) + 3 * H * S)
        ),
        "same_nonlocal_scale_only": s.diff(full[0], R) + s.diff(full[0], Q),
        "same_full_nonlocal_scale_normalization": s.diff(full[1], Q) - 3,
        "no_new_fourth_time_coefficient": s.diff(full[0], s.diff(eta, u, 4)),
        "no_new_third_scale_coefficient": s.diff(full[0], s.diff(w, u, 3)),
        "pure_scale_local_nine_P_contact_now_cancelled": s.expand(
            full[1].subs({eta: 0, s.diff(eta, u): 0}) - 3 * Q
        ),
    }
    # Quantify the first-output derivative of I4 of the ACTUAL local row.
    eps = profile.EPS
    envelopes = {
        (0, 2): {0: 2 * eps, 1: 2 * eps, 2: 2 * eps},
        (0, 1): {0: eps, 1: eps},
        (0, 0): {0: 24 * eps},
        (1, 1): {0: 6 * eps, 1: 6 * eps},
    }
    quantum = 0
    for (_, order), row in envelopes.items():
        for k, value in row.items():
            degree = 3 - order + k
            if degree >= 1:
                quantum += s.binomial(order, k) * value / s.factorial(degree - 1)
    checks["actual_new_first_row_quantum_envelope"] = quantum - s.Rational(59, 2) * eps
    return {
        "normalization": "rho,P and R,Q are the actual physical conditional vector stresses and pure-scale variations divided by kappa. This simply absorbs the old gamma times64pi^2 normalization.",
        "new_complete_stationary_adapted_response": s.ImmutableMatrix(full),
        "new_local_quantum_plus_fixed_profile_density": density,
        "new_local_time_channel": eta_local,
        "new_local_scale_addition": scale_local,
        "literal_fixed_profile_quadratic_density": added_density,
        "actual_new_first_row_coefficient_envelopes": envelopes,
        "actual_new_first_row_quantum_upper": quantum,
        "retained_inverse": "The S6.180 actual classical tree is unchanged as the reference part; the explicit new scalar action is included in the local response above. The fourth-order A=-6delta^2 and the matched scalar F_m are unchanged. Only lower local kernels in the Volterra remainder change, so the same proof gives a unique smooth prepared retained causal response on the slab.",
        "same_massive_kernel_L1_upper": kernel_norm.K_UPPER,
        "new_local_lapse_recovery_C1_upper": bounds.data()[
            "new_retained_first_row_C1_upper"
        ],
        "uncomputed": "The new complete curved weak-log majorant and a useful/small full inverse norm remain unevaluated. Mean stationarity does not imply stability or nonlinear response control.",
        "checks": checks,
        "gates": {
            "same_scalar_kernel_norm_bound": kernel_norm.K_UPPER == 300000000000,
            "actual_new_C1_below_15000": bounds.data()[
                "new_retained_first_row_C1_upper"
            ]
            < 15000,
        },
    }
