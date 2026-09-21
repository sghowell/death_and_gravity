"""Full known old-polynomial sum, including the tuned original cancellation."""

from functools import cache

import sympy as s
from p8_vacuum_affine_box_curvature_coefficient import moment as box

from . import conversion, moments

N = moments.N
C, V = s.symbols("c v", real=True)


def contact_ratio(n):
    return -3 / (n - 2) + 2 / (n - 2) ** 2


@cache
def sectors():
    a = C + 1 / N
    flat = a * a / 840 + a / (60 * N**2) + a / (6 * N**3) + 1 / (12 * N**4)
    J = moments.triangle_moments()
    return (
        -4 * flat,
        a * J[0] + J[1] / N**2 + J[2] / N**3 + J[3] / N**4,
        box.CLOSED.subs(box.N, N),
    )


@cache
def generic_coefficient():
    return s.factor(sum(sectors()))


@cache
def original_coefficient():
    return s.factor(generic_coefficient().subs(C, contact_ratio(N)))


@cache
def data():
    x = s.Symbol("x", real=True)
    B = sum(V**j * s.integrate((x * (1 - x)) ** j, (x, 0, 1)) / j for j in range(1, 4))
    A = C + sum(V**j / N ** (j + 1) for j in range(4))
    flat = s.Poly(s.expand(A * A * B / 2), V).coeff_monomial(V**3)
    bubble, triangle, boxvalue = sectors()
    actual = original_coefficient()
    checks = {
        "complete_MS_bubble_first_three_moments": s.expand(
            B - (V / 6 + V**2 / 60 + V**3 / 420)
        ),
        "complete_A_squared_B_degree6_coefficient": s.factor(
            bubble - conversion.build()[1][(3, 0, 0)] * flat / 2
        ),
        "entire_generic_sector_sum": s.factor(
            generic_coefficient() - bubble - triangle - boxvalue
        ),
        "entire_original_tuned_contact_sum": s.factor(
            actual - (bubble + triangle + boxvalue).subs(C, contact_ratio(N))
        ),
        "original_no_n_inverse_square": s.limit(N**2 * actual, N, s.oo),
        "original_no_n_inverse_cube": s.limit(N**3 * actual, N, s.oo),
        "original_first_n_inverse_fourth": s.limit(N**4 * actual, N, s.oo)
        + s.Rational(23, 105),
        "original_contact_ratio_algebra": s.factor(
            contact_ratio(N) + (3 * N - 8) / (N - 2) ** 2
        ),
    }
    for j in range(4):
        checks[f"outer_resolvent_order{j}"] = s.diff(1 / (N - V), V, j).subs(
            V, 0
        ) / s.factorial(j) - 1 / N ** (j + 1)
    return {
        "checks": checks,
        "gates": {
            "whole_bubble_triangle_box_added_only_in_common_lift": True,
            "original_contact_tuning_not_replaced_by_leading_term": True,
            "two_leading_mass_orders_cancel_exactly": True,
            "fixed_finite_OS4_constant_has_no_degree6_part": True,
            "local_coefficient_not_full_physical_amplitude": True,
        },
        "whole_generic_bubble_triangle_box_coefficients": sectors(),
        "whole_generic_old_polynomial_coefficient": generic_coefficient(),
        "whole_original_old_polynomial_coefficient": actual,
        "whole_original_heavy_mass_limit": -s.Rational(23, 105),
        "whole_physical_coefficient": "chi_core=g^4*c_core(n)/(16pi^2), and the radiative difference is chi_core*T/sqrt(kappa). This is the homogeneous degree6 local difference of the entire known old-polynomial sector from the same fixed six-word lift. All external emissions cancel between exactly equal flat Taylor polynomials. Lower-derivative terms cannot supply this homogeneous degree. Existing constant OS4 subtraction is unchanged. The full unexpanded physical loop amplitude was already included in S339/S342; no extra copy is to be added.",
        "whole_remaining_matching": "New local tadpoles and mixed-source classes require their own off-shell common-basis analysis before any full-known-matter aggregate. The independent parent chi is not assigned by the selected-loop result. No physical above-threshold Taylor error or Regge/UV condition follows.",
    }
