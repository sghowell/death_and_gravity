"""Continuous rational coefficient and prepared real-time Gronwall bounds."""

from functools import cache

import sympy as s

from . import hamiltonian as ham

C0 = s.Integer(53000)
Y0 = s.Integer(15000)


@cache
def data():
    inv = 20
    H = th = s.Rational(8, 5)
    d = s.Rational(1, 2)
    e = s.Rational(1, 10)
    w = s.Rational(1, 20)
    c = s.Rational(3, 200)
    A11 = th * c * inv
    A12 = th**2 * inv
    A21 = c**2 * inv + 9 * e**2
    A22 = 3 * H + A11
    rows = (A11 + A12 / 10, 10 * A21 + A22)
    force = (th * inv * (1 + d), 10 * (c * inv * (1 + d) + 1))
    M = s.Integer(7)
    B = s.Integer(48)
    y = B * (3**7 - 1) / M
    n = inv * ((th / 10 + c) * Y0 + (1 + d))
    v = Y0 + d * n
    sigma = w * n + 3 * e * Y0
    Jlower = s.Rational(1215, 800) / s.Rational(5, 4) ** 18
    checks = {
        "actual_real_J_denominator": ham.J * 800 * (1 + ham.u**2) ** 18 - ham.P,
        "L_absolute_bound": c - 3 * e * w,
        "row_one_bound": rows[0] - s.Rational(28, 5),
        "row_two_bound": rows[1] - s.Rational(249, 40),
        "force_row_one": force[0] - 48,
        "force_row_two": force[1] - s.Rational(29, 2),
        "prepared_state_exponential_majorant": y - s.Rational(104928, 7),
        "lapse_reconstruction_envelope": n - 52530,
        "physical_scale_reconstruction_envelope": v - 41265,
        "original_matter_rate_envelope": sigma - s.Rational(14253, 2),
    }
    return {
        "interval": (s.Rational(-1, 2), s.Rational(1, 2)),
        "norm": "Joint maximum of both physical force components; phase max(abs(vhat),10abs(p))",
        "continuous_basic_envelopes": {
            "H": H,
            "theta": th,
            "delta": d,
            "ell": e,
            "w0": w,
            "L": c,
            "one_over_two_J": inv,
        },
        "J_lower": Jlower,
        "weighted_phase_row_upper": rows,
        "weighted_force_row_upper": force,
        "phase_coefficient_upper": M,
        "phase_force_upper": B,
        "state_C0_upper": Y0,
        "physical_lapse_C0_upper": n,
        "physical_log_scale_C0_upper": v,
        "matter_rate_and_prepared_field_C0_upper": sigma,
        "classical_physical_C0_to_C0_upper": C0,
        "checks": {key: s.factor(value) for key, value in checks.items()},
        "gates": {
            "actual_J_above_one_fortieth": Jlower > s.Rational(1, 40),
            "inverse_two_J_below_twenty": 1 / (2 * Jlower) < 20,
            "subtracted_A12_absolute_enclosed": s.Rational(1, 6) < A12,
            "weighted_phase_norm_below_seven": max(rows) < M,
            "weighted_force_norm_at_most_forty_eight": max(force) <= B,
            "strict_state_bound": y < Y0,
            "strict_physical_C0_bound": max(n, v) < C0,
            "matter_rate_below_seventy_two_hundred": sigma < 7200,
        },
    }
