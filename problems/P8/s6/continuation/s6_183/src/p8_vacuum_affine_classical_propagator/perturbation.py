"""First formal loop coefficient with the actual nonadaptive stationary response."""

from functools import cache

import sympy as s
from p8_vacuum_affine_metric_response import bounds as physical
from p8_vacuum_affine_quantum_retuning import profile, response, ward

from . import jet_bounds as jets
from . import real_bounds as real

Q_UPPER = 7 * s.Rational(1, 10**770)
BQ_UPPER = 4 * s.Rational(1, 10**765)
FIRST_UPPER = 4 * s.Rational(1, 10**669)


@cache
def data():
    n, v, dr, dp = s.symbols("physical_n physical_v delta_r delta_p", real=True)
    r, p = profile.rho, profile.P
    S = r + p
    old = s.Matrix([-dr - 3 * r * v, 3 * dp + 3 * p * (n + 3 * v)])
    fixed = ward.data()["profile_current_Hessian"] * s.Matrix([n, v])
    actual = s.Matrix([-dr - S * n, 3 * dp + 3 * S * n])
    # This physical calculation is independent of the adapted Ward route.
    eta, w = response.eta, response.w
    mapped = actual.subs(
        {
            n: s.diff(eta, profile.u),
            v: w + ward.H * eta,
            dr: response.R + s.diff(r, profile.u) * eta,
            dp: response.Q + s.diff(p, profile.u) * eta,
        }
    )
    adapted = s.Matrix([-ward.D(mapped[0]) + ward.H * mapped[1], mapped[1]])
    prior = response.data()["new_complete_stationary_adapted_response"]
    R = max(
        physical.data()[
            "actual_conditional_physical_stress_C10_to_C0_over_kappa"
        ].values()
    )
    eps = profile.EPS
    rows = (R + 2 * eps, 3 * R + 6 * eps)
    beta = real.C0 * Q_UPPER
    first = beta * jets.C10
    checks = {
        "literal_fixed_profile_plus_unretuned_vector": s.ImmutableMatrix(
            (old + fixed - actual).applyfunc(s.expand)
        ),
        "actual_stationary_adapted_parent_bridge": s.ImmutableMatrix(
            (adapted - prior).applyfunc(response.simplify_response)
        ),
        "retained_nonlocal_rho_current": s.diff(actual[0], dr) + 1,
        "retained_nonlocal_pressure_current": s.diff(actual[1], dp) - 3,
        "fixed_profile_not_live_response": s.diff(fixed[0], dr) + s.diff(fixed[1], dp),
        "phase_followed_by_quantum_envelope": beta - real.C0 * Q_UPPER,
        "first_formal_coefficient_envelope": first - real.C0 * Q_UPPER * jets.C10,
        "same_anchor_epsilon": eps - s.Rational(1, 10**770),
    }
    return {
        "new_named_parent": "CD-REG-AFFINE-ISO-QG1",
        "actual_complete_stationary_physical_response": s.ImmutableMatrix(actual),
        "actual_source_pinned_physical_stress_response_upper": R,
        "physical_stationary_current_row_upper": rows,
        "stationary_quantum_C10_to_C0_upper": Q_UPPER,
        "classical_after_quantum_C10_to_C0_upper": BQ_UPPER,
        "first_formal_response_C10_force_to_C0_metric_upper": FIRST_UPPER,
        "formal_equations": "(T0+alpha Qstar)(h0+alpha h1+...)=g; T0 h0=g, T0 h1=-Qstar h0, prepared data at every coefficient",
        "marker_scope": "alpha multiplies BOTH the fixed new scalar coefficient and the specified conditional Gaussian vector functional; reference mean is stationary at every formal alpha",
        "future_extension": "A prepared response need not vanish at the right observation endpoint. A smooth future extension and later cutoff cannot alter earlier mode ODEs, local jets or retarded stress; all bounds use only the observation/past interval.",
        "not_an_error_bound": "h1 is a uniquely defined formal coefficient. No differentiable exact solution branch in alpha, finite-alpha remainder, C0 contraction, pole deletion or stability theorem is asserted.",
        "checks": checks,
        "gates": {
            "actual_stress_below_five_e_minus_795": 0 < R < 5 * s.Rational(1, 10**795),
            "complete_stationary_current_below_Q": max(rows) < Q_UPPER,
            "classical_quantum_composition_below_display": beta < BQ_UPPER,
            "first_formal_coefficient_below_display": first < FIRST_UPPER,
        },
    }
