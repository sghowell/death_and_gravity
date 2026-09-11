"""Leading-source below-pole estimate and the generic real-time norm limitation."""

from functools import cache

import sympy as s

from .vacuum import KAPPA, MAX_ZETA, N


def require_finite_band(radius, zeta):
    for v in (radius, zeta):
        if isinstance(v, bool) or not isinstance(v, (int, s.Integer, s.Rational)):
            raise TypeError("Finite exact real radius and zeta required")
    if radius < 0 or zeta <= 0 or zeta * radius**2 >= 1:
        raise ValueError("A strict below-pole Euclidean momentum ball is required")
    return s.Rational(radius), s.Rational(zeta)


def inverse_upper(radius, zeta):
    b, z = require_finite_band(radius, zeta)
    return (1 + z * b * b) / (1 - z * b * b)


def remainder_upper(radius, zeta):
    b, z = require_finite_band(radius, zeta)
    return 2 * z * b * b / (1 - z * b * b)


@cache
def data():
    inv = inverse_upper(4, MAX_ZETA)
    remainder = remainder_upper(4, MAX_ZETA)
    action = KAPPA * remainder * (64 * N / KAPPA**2) ** 2 / 2
    z, epsilon = s.symbols("zeta epsilon", positive=True)
    pole_eigenvalue = 1 - z * (1 / z + epsilon)
    return {
        "leading_source_support_radius": 4,
        "leading_inverse_Euclidean_component_norm_upper": inv,
        "leading_I_minus_inverse_norm_upper": remainder,
        "leading_degree_eight_Lorentz_source_action_absolute_upper_divided_by_U_squared": action,
        "leading_decimal_upper": "1e-2392",
        "estimate": "The complete flat inverse is used on S4; |p_M^2|<=|p_E|^2<=16 and ||p_sharp tensor p_flat||_2<=16. No derivative expansion is used for this finite-germ estimate.",
        "generic_real_time_limit": "Transverse source packets with p_M^2-mass^2 in [epsilon,2epsilon] have inverse norm at least 1/(2 zeta epsilon). Such smooth compact Fourier packets can be reflected to give real Schwartz sources. Thus there is no global unweighted L2 inverse bound over arbitrary sources.",
        "not_an_image_no_go": "The generic packet argument does not prove the actual nonlinear image S(Phi) contains those packets, does not exclude this parent, and does not replace a weighted/state/contour bound for the actual full source.",
        "checks": {
            "leading_inverse_constant": inv - s.Rational(63, 62),
            "leading_remainder_constant": remainder - s.Rational(1, 62),
            "leading_source_action_constant": action
            - s.Rational(1024, 31) * N**2 / KAPPA**3,
            "transverse_pole_eigenvalue": pole_eigenvalue + z * epsilon,
            "strict_leading_pole_margin": 1 - 16 * MAX_ZETA - s.Rational(124, 125),
        },
        "gates": {
            "leading_band_strictly_below_minimum_vector_mass": bool(16 * MAX_ZETA < 1),
            "leading_Lorentz_action_numerical_bound": bool(
                action < s.Rational(1, 10**2392)
            ),
        },
    }
