"""Actual full-R chart jets and nonstationary mean-current pullback."""

from functools import cache

import sympy as s
from p8_affine_vacuum_domain import family
from p8_vector_state import wkb

u = wkb.u


@cache
def data():
    h = (1 + u * u) ** 3
    N = s.symbols("N", positive=True)
    # g_physical = R^(-1/2) g_hat + (1-R^(-1/2)) du du/X.
    # N is unchanged in clock gauge, v = vhat - log R(u,N^-2)/4.
    R = family.data()["R"]
    full_jet = [
        s.simplify(s.diff(R, family.X, j).subs(family.X, 1)).subs(family.u, u)
        for j in range(3)
    ]
    Rn = -2 * full_jet[1]
    Rnn = 6 * full_jet[1] + 4 * full_jet[2]
    w1 = s.factor(-Rn / 4)
    w2 = s.factor(-(Rnn - Rn * Rn) / 4)
    rho, P, drho, dP, n, vhat = s.symbols(
        "rho pressure delta_rho delta_pressure n vhat", real=True
    )
    J = s.Matrix([[1, 0], [w1, 1]])
    physical_contact = s.Matrix([[0, -3 * rho], [3 * P, 9 * P]])
    second = s.diag(3 * P * w2, 0)
    contact = s.simplify(J.T * physical_contact * J + second)
    out = s.simplify(J.T * s.diag(-1, 3))
    q = s.Matrix([n, vhat])
    delta_stress = s.Matrix([drho, dP])
    response = out * delta_stress + contact * q
    v = vhat + w1 * n
    direct = s.Matrix(
        [
            -drho - 3 * rho * v + w1 * (3 * dP + 3 * P * (n + 3 * v)) + 3 * P * w2 * n,
            3 * dP + 3 * P * (n + 3 * v),
        ]
    )
    expected = s.Matrix(
        [
            [3 * P * (w1 + 3 * w1 * w1 + w2) - 3 * rho * w1, -3 * rho + 9 * w1 * P],
            [3 * P + 9 * w1 * P, 9 * P],
        ]
    )
    r = s.symbols("reciprocal_h", real=True)
    omega_clock = -s.log(1 + (N**-2 - 1) / h) / 4
    checks = {
        "full_R_clock_value": full_jet[0] - 1,
        "full_R_clock_first_X": s.factor(full_jet[1] - 1 / h),
        "full_R_clock_second_X": full_jet[2],
        "full_chart_first_N": s.factor(w1 - 1 / (2 * h)),
        "full_chart_second_N": s.factor(w2 + 3 / (2 * h) - 1 / h**2),
        "clock_jet_not_global_chart_replacement_first": s.factor(
            s.diff(omega_clock, N).subs(N, 1) - w1
        ),
        "clock_jet_not_global_chart_replacement_second": s.factor(
            s.diff(omega_clock, N, 2).subs(N, 1) - w2
        ),
        "full_nonstationary_current_chain_rule": s.ImmutableMatrix(
            (response - direct).applyfunc(s.factor)
        ),
        "full_contact_matrix": s.ImmutableMatrix(
            (contact - expected).applyfunc(s.factor)
        ),
        "omega_second_sharp_lower_square": s.expand(
            r * r - 3 * r / 2 + s.Rational(9, 16) - (r - s.Rational(3, 4)) ** 2
        ),
        "omega_first_plus_second_combination": s.factor(w1 + w2 - (1 / h**2 - 1 / h)),
    }
    return {
        "physical_metric_map": "g_physical=R^(-1/2)g_hat+(1-R^(-1/2))du du/X, not R-1/2. The complete analytic R is unchanged; only its exact reference-clock jets enter this linear response.",
        "omega_N": w1,
        "omega_NN": w2,
        "linear_input_map": s.ImmutableMatrix(J),
        "normalized_physical_stress_to_clock_current_map": s.ImmutableMatrix(out),
        "physical_volume_readout_contacts": s.ImmutableMatrix(physical_contact),
        "nonzero_second_metric_map_contact": s.ImmutableMatrix(second),
        "complete_nonstationary_clock_contact": s.ImmutableMatrix(contact),
        "clock_current_background": s.ImmutableMatrix(J.T * s.Matrix([-rho, 3 * P])),
        "clock_retarded_response": s.ImmutableMatrix(response),
        "normalization": "Currents and their variations are divided by the FIXED reference a0(t)^3, then by kappa for the bounds. This is not a background-cancelled current, moving-volume physical stress, or an ordinary symmetric retarded-action Hessian.",
        "checks": checks,
    }


@cache
def envelopes():
    w1 = data()["omega_N"]
    coefficients = {j: wkb.box_bound(s.factor(s.diff(w1, u, j))) for j in range(11)}
    rows = {
        j: 1
        + sum(
            s.binomial(j, k) * coefficients[j - k]["absolute_upper"]
            for k in range(j + 1)
        )
        for j in range(11)
    }
    M = max(s.S.One, *rows.values())
    checks = {
        "coefficient_reconstruction_" + str(j): v["reconstruction"]
        for j, v in coefficients.items()
    }
    checks["continuous_C10_input_constant"] = M - s.Rational(46090764897, 8)
    return {
        "coefficient_derivative_envelopes": coefficients,
        "Leibniz_C10_rows": rows,
        "C10_input_upper": M,
        "C0_stress_output_upper": s.Integer(3),
        "C0_physical_source_upper": s.Rational(3, 2),
        "omega_N_upper": s.Rational(1, 2),
        "omega_NN_absolute_upper": s.Rational(9, 16),
        "checks": checks,
        "gates": {
            "full_slab_chart_coefficients_finite": all(
                v["absolute_upper"] > 0 for v in coefficients.values()
            ),
            "C10_input_bound_positive": M > 0,
        },
    }
