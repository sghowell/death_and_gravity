"""Corner-subtracted sunset and a complex-regulator finite-part enclosure."""

from functools import cache

import sympy as s


@cache
def data():
    r, t, e = s.symbols("r t epsilon", positive=True)
    a, b, c = s.symbols("mass_a_squared mass_b_squared mass_c_squared", positive=True)
    A = 1 + r + r * t
    V = 1 + t + r * t
    B = a + r * b + r * t * c
    x, y, z = 1 / A, r / A, r * t / A
    U = s.factor(x * y + x * z + y * z)
    jac = s.factor(
        s.det(s.Matrix([[s.diff(y, r), s.diff(y, t)], [s.diff(z, r), s.diff(z, t)]]))
    )
    rho = s.Rational(1, 16)
    J = (2 ** (e - 1) - 1) / (e - 1)
    return {
        "sector_coordinates": {
            "x": x,
            "y": y,
            "z": z,
            "jacobian": jac,
            "U": U,
            "mass_form": B / A,
        },
        "one_ordered_sector_integrand": r ** (e - 1) * V ** (e - 2) * B ** (1 - 2 * e),
        "corner_reference_integral": a ** (1 - 2 * e) * J / e,
        "continued_sector": "a^(1-2e) J(e)/e + integral_0^1 dr dt r^(e-1)[V^(e-2)B^(1-2e)-(1+t)^(e-2)a^(1-2e)]. Sum all six labeled mass permutations.",
        "sunset_normalization": "S_D=Q^(-2) exp(2 gamma_E epsilon) mF^(4 epsilon) Gamma(2 epsilon-1) A(epsilon). Q=16pi^2. The subtracted compact representation continues to Re(epsilon)>-1.",
        "regulator_circle_radius": rho,
        "circle_bounds": {
            "sector_sum_multiplier": s.Integer(120),
            "gamma_and_exponential_multiplier": s.Integer(40),
            "sunset_multiplier": s.Integer(4800),
            "tadpole_multiplier": s.Integer(70),
            "mixed_bubble_multiplier": s.Integer(64),
        },
        "checks": {
            "simplex_coordinates_sum_to_one": s.factor(x + y + z - 1),
            "sector_simplex_jacobian": s.factor(jac - r / A**3),
            "sector_U_formula": s.factor(U - r * V / A**2),
            "mass_form_transform": s.factor(a * x + b * y + c * z - B / A),
            "all_A_powers_cancel": (4 - 2 * e) + (-1 + 2 * e) - 3,
            "all_r_powers_combine": (-2 + e) + 1 - (e - 1),
            "corner_t_integral": s.diff((1 + t) ** (e - 1) / (e - 1), t)
            - (1 + t) ** (e - 2),
            "corner_at_zero_epsilon": s.limit(J, e, 0) - s.Rational(1, 2),
            "derivative_majorant_coefficient": (2 + rho)
            + s.Rational(2, 3) * (1 + 2 * rho)
            - s.Rational(45, 16),
            "strict_derivative_majorant_gap": 3
            - s.Rational(45, 16)
            - s.Rational(3, 16),
            "six_sector_remainder_multiplier": 6 * 3 / (1 - rho) - s.Rational(96, 5),
            "pole_and_remainder_below_120": 120
            - (6 / rho + s.Rational(96, 5))
            - s.Rational(24, 5),
            "gamma_shift_denominator_lower": 2 * rho * (1 - 2 * rho)
            - s.Rational(7, 64),
            "sunset_multiplier_product": 40 * 120 - 4800,
            "tadpole_multiplier_gap": 70 - 4 / (rho * (1 - rho)) - s.Rational(26, 15),
        },
        "scope": "A Cauchy bound on the finite Laurent coefficient of the specified scalar forest. No finite part is set to zero; pole times positive-epsilon products are bounded together. No physical dimensional or momentum cutoff is introduced.",
    }
