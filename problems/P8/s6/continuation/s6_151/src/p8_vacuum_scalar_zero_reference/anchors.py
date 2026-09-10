"""Full regulated recursive and proper-MS local anchors, including pole products."""

from functools import cache

import sympy as s


@cache
def data():
    e, ell, j, Q, mass = s.symbols("epsilon ell j Q scalar_mass_squared")
    E = 1 + 2 * ell * e + 2 * ell**2 * e**2
    gamma2 = 1 / (2 * e) + s.zeta(2) * e
    A = 3 / e + 3 + 6 * j * e
    bare = s.series(E * gamma2 * A / 3, e, 0, 1).removeO().expand()
    I0 = 1 / e + ell + (ell**2 / 2 + s.zeta(2) / 2) * e
    I0sq = s.series(I0**2, e, 0, 1).removeO().expand()
    recursive = s.expand(bare - I0sq)
    properMS = s.expand(bare - I0 / e)
    recursive = s.series(recursive, e, 0, 1).removeO().expand()
    properMS = s.series(properMS, e, 0, 1).removeO().expand()
    fin_recursive = j + ell - ell**2
    fin_MS = j + ell + ell**2 / 2 + s.pi**2 / 12
    sunset = mass ** (1 - 2 * e) * s.gamma(2 * e) / (2 * e - 1)
    derivative = -s.diff(sunset, mass) / 3
    return {
        "symbols": {"epsilon": e, "ell": ell, "j": j, "Q": Q},
        "Q_squared_bare_zero_wineglass": bare,
        "Q_squared_old_inner_subtracted_J0": recursive,
        "Q_squared_proper_MS_zero_reference": properMS,
        "finite_recursive_J0": fin_recursive / Q**2,
        "finite_proper_MS_J0": fin_MS / Q**2,
        "required_I0_evanescent_coefficient": ell**2 / 2 + s.pi**2 / 12,
        "checks": {
            "equal_mass_derivative_and_Gamma_recurrence": s.factor(
                derivative - mass ** (-2 * e) * s.gamma(2 * e) / 3
            ),
            "bare_double_simple_and_finite_terms": s.expand(
                bare
                - 1 / (2 * e**2)
                - (ell + s.Rational(1, 2)) / e
                - ell**2
                - ell
                - j
                - s.pi**2 / 6
            ),
            "recursive_double_simple_and_finite_terms": s.expand(
                recursive
                + 1 / (2 * e**2)
                - (s.Rational(1, 2) - ell) / e
                - fin_recursive
            ),
            "proper_MS_double_simple_and_finite_terms": s.expand(
                properMS + 1 / (2 * e**2) - 1 / (2 * e) - fin_MS
            ),
            "proper_MS_poles_independent_of_scale": s.diff(
                -1 / (2 * e**2) + 1 / (2 * e), ell
            ),
            "finite_inner_scheme_difference": s.expand(
                fin_MS - fin_recursive - 3 * ell**2 / 2 - s.pi**2 / 12
            ),
            "I0_square_finite_evanescent_product_retained": s.expand(
                I0sq - 1 / e**2 - 2 * ell / e - 2 * ell**2 - s.pi**2 / 6
            ),
            "early_I0_truncation_defect": s.expand(
                fin_MS
                - (bare - (1 / e + ell) / e).coeff(e, 0)
                + ell**2 / 2
                + s.pi**2 / 12
            ),
            "exact_common_scale_factor": s.diff(s.exp(2 * ell * e), e, 2).subs(e, 0) / 2
            - 2 * ell**2,
        },
        "scope": "The proper-MS anchor subtracts I_MS I0, not I0 squared. Its finite pi-squared/12 and half-log-squared term retain the epsilon coefficient of I0 in a pole product. The old recursive reference and the actual MS reference are distinguished explicitly.",
    }
