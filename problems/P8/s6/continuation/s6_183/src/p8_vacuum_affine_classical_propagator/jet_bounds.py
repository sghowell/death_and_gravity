"""Uniform complex discs and an explicit no-loss classical C10 bound."""

from functools import cache

import sympy as s

from . import hamiltonian as ham
from . import real_bounds as real

RADIUS = s.Rational(1, 10**8)
C10 = s.Integer(10**96)


@cache
def data():
    u = ham.u
    coeffsum = sum(power[0] * abs(value) for power, value in s.Poly(ham.P, u).terms())
    variation = RADIUS * coeffsum
    inv = s.Rational(2, 5) * s.Rational(13, 10) ** 18
    th, H, d, e, c = (
        s.Integer(6),
        s.Integer(3),
        s.Rational(3, 2),
        s.Integer(1),
        s.Integer(18),
    )
    z = s.Rational(51, 100)
    denom = s.Rational(7, 10)
    A11 = th * c * 50
    A12 = th**2 * 50 + s.Rational(1, 6)
    A21 = c**2 * 50 + 9 * e**2
    A22 = 3 * H + A11
    rows = (A11 + A12 / 10, 10 * A21 + A22)
    frows = (th * 50 * (1 + d), 10 * (c * 50 * (1 + d) + 1))
    M, B = s.Integer(200000), s.Integer(30000)
    cr, dr = s.Integer(1000), s.Integer(125)
    if not all(value > 0 for value in s.Poly(ham.P, u).all_coeffs() if value != 0):
        raise ValueError("The real positive polynomial proof no longer applies")
    state = {j: real.Y0 * s.factorial(j) / RADIUS**j for j in range(11)}
    lapse = {
        j: (j + 1) * s.factorial(j) / RADIUS**j * (cr * real.Y0 + dr) for j in range(11)
    }
    physical = {
        j: s.factorial(j)
        / RADIUS**j
        * (real.Y0 + s.Rational(3, 4) * (cr * real.Y0 + dr) * (j + 1) * (j + 2))
        for j in range(11)
    }
    checks = {
        "positive_polynomial_derivative_sum": coeffsum - 1363114784,
        "disc_polynomial_variation": variation - s.Rational(42597337, 3125000),
        "complex_inverse_two_J": inv - 400 * s.Rational(13, 10) ** 18 / 1000,
        "weighted_reconstruction_crude_bound": 50 * (c + th / 10) - 930,
        "direct_force_reconstruction_bound": 50 * (1 + d) - 125,
        "Leibniz_triangle_sum": sum(j + 1 for j in range(11)) - 66,
    }
    # This recurrence bounds the scaled derivatives without differentiating
    # expanded high-degree rational functions ten times.
    scaled = [real.Y0]
    for j in range(10):
        upper = (
            RADIUS
            / (j + 1)
            * (
                M * sum(scaled)
                + B * sum(RADIUS**i / s.factorial(i) for i in range(j + 1))
            )
        )
        if not upper < real.Y0:
            raise ValueError(
                "Scaled derivative induction failed at order " + str(j + 1)
            )
        scaled.append(real.Y0)
    return {
        "complex_disc_radius": RADIUS,
        "positive_polynomial_degree_coefficient_sum": coeffsum,
        "polynomial_disc_variation_upper": variation,
        "polynomial_disc_modulus_lower": s.Integer(1000),
        "complex_two_J_inverse_upper": s.Integer(50),
        "complex_weighted_phase_rows": rows,
        "complex_weighted_force_rows": frows,
        "complex_phase_matrix_upper": M,
        "complex_force_matrix_upper": B,
        "complex_lapse_state_row_upper": cr,
        "complex_lapse_force_row_upper": dr,
        "phase_derivative_upper": state,
        "physical_lapse_derivative_upper": lapse,
        "physical_log_scale_derivative_upper": physical,
        "classical_physical_C10_to_C10_upper": C10,
        "scaled_derivative_recurrence": "x[j+1] <= radius/(j+1)*(M*sum(x[0:j+1])+B*G*sum(radius^i/i!,i=0..j)); x0<=15000G",
        "checks": checks,
        "gates": {
            "disc_within_fifty_one_hundredths": s.Rational(1, 2) + RADIUS < z,
            "complex_denominator_lower": 1 - z * z > denom,
            "complex_denominator_upper": 1 + z * z < s.Rational(13, 10),
            "disc_polynomial_modulus_above_one_thousand": 1215 - variation > 1000,
            "complex_inverse_two_J_below_fifty": inv < 50,
            "complex_H_below_three": 4 * z / denom < 3,
            "complex_delta_below_three_halves": 1 / (2 * denom**3) < d,
            "complex_ell_below_one": 1 / (10 * denom**6) < e,
            "complex_w_below_six": e * (3 * d + 1) < 6,
            "complex_theta_below_six": 3 + z / denom**4 < th,
            "complex_weighted_phase_below_majorant": max(rows) < M,
            "complex_weighted_force_below_majorant": max(frows) < B,
            "complex_lapse_state_below_majorant": 50 * (c + th / 10) < cr,
            "scaled_derivative_induction_closes": RADIUS * (M * real.Y0 + B) < real.Y0,
            "all_lapse_jets_below_C10_bound": max(lapse.values()) < C10,
            "all_physical_scale_jets_below_C10_bound": max(physical.values()) < C10,
        },
    }
