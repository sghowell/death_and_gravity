"""Full actual analytic source bounds on the specified small vacuum field class."""

from functools import cache

import sympy as s

N = 1024
KAPPA = s.Integer(10) ** 800
MAX_ZETA = s.Rational(1, 2000)


def require_small_domain(n, kappa, zeta):
    for value in (n, kappa, zeta):
        if isinstance(value, bool) or not isinstance(
            value, (int, s.Integer, s.Rational)
        ):
            raise TypeError("Finite exact integer/rational parameters required")
    if not isinstance(n, (int, s.Integer)) or n < 4 or n % 2:
        raise ValueError("The switch exponent must be even and at least four")
    if kappa < 64 * n or not 0 < zeta <= MAX_ZETA:
        raise ValueError(
            "The certified small-field and positive vector range is required"
        )
    return s.Integer(n), s.Rational(kappa), s.Rational(zeta)


def component_bound(n, kappa):
    return (
        64 * n / kappa**2
        + 384 * n / kappa**3
        + 6144 * n**2 / kappa**4
        + 9216 * n**2 / kappa**5
    )


@cache
def data():
    n, k, x, u = s.symbols("n kappa X u", positive=True)
    rho = 1 / (16 * n)
    # Sharp enough intermediate B bound avoids dropping the factor |1-X|.
    delta_factor = (n + 1) * (1 + rho)
    rx_factor = (n + 1) * rho + 3 * n * (1 + rho)
    normalized = s.cancel(component_bound(n, k) / (n / k**2))
    upper_at_minimum = normalized.subs({n: 4, k: 256}, simultaneous=True)
    euclidean = s.Integer(32768) * N**2 / KAPPA**3
    exact = component_bound(s.Integer(N), KAPPA)
    gradient, box, z = s.symbols("dPhi BoxPhi ZPhi", real=True)
    eps = s.symbols("epsilon")
    Xphi = s.symbols("X_Phi", real=True)
    # r=(R-1)/X^2=-n at the vacuum; terms containing H or Ru start later.
    leading = s.expand(
        eps * gradient * (eps**2 * Xphi * (-n) * eps * box - (-n) * eps**3 * z)
    )
    return {
        "class": "Real Schwartz Phi on R^4, Fourier support in the Euclidean momentum unit ball and Fourier L1 at most one, with U=||Phi||_2; u=Phi/sqrt(kappa). This is a Lorentzian real-spacetime class, not a Wick rotation.",
        "parameters": {"n": N, "kappa": KAPPA, "zeta_upper": MAX_ZETA},
        "small_domain": {"kappa_lower": 64 * n, "X_absolute_upper": 4 / k, "rho": rho},
        "full_function_bounds": {
            "B_upper": (n + 1) * x**2,
            "B_X_absolute_upper": 3 * n * s.Abs(x),
            "R_minus_one_absolute_upper": 2 * n * x**2 / (1 + u * u) ** 3,
            "regular_r_absolute_upper": 2 * n / (1 + u * u) ** 3,
            "R_X_absolute_upper": 4 * n * s.Abs(x) / (1 + u * u) ** 3,
            "R_lower": 1 - 2 * n * rho**2,
            "R_u_identity": "-6u(R-1)/(1+u^2)",
        },
        "per_component_source_L2_upper_divided_by_U": exact,
        "simple_full_source_vector_L2_upper_divided_by_U": 256 * N / KAPPA**2,
        "Euclidean_positive_source_action_upper_divided_by_U_squared": euclidean,
        "Euclidean_decimal_upper": "1e-2389",
        "Euclidean_observable_change": "The independent real Euclidean functional uses X_L=-X_E, Box_L=-Box_E, Z_L=Z_E and its Euclidean one-form source. The same absolute pointwise estimates apply. This positive action bound is NOT a Lorentzian real-time matching bound.",
        "leading_vacuum_source": "S4_mu=(n/kappa^2) partial_mu Phi [Z_Phi-X_Phi Box Phi], degree four with five derivatives",
        "leading_vector_L2_upper_divided_by_U": 64 * N / KAPPA**2,
        "full_source_frequency_support": "The complete analytic source has infinitely many even field degrees; no finite Fourier-support radius is inferred for it. Only S4 has radius at most four on the stated class.",
        "checks": {
            "delta_bound_factor_margin": s.factor(
                2 * n - delta_factor - (16 * n * n - 17 * n - 1) / (16 * n)
            ),
            "RX_bound_factor_margin": s.factor(
                4 * n - rx_factor - (16 * n * n - 4 * n - 1) / (16 * n)
            ),
            "R_lower_margin_above_half": s.factor(
                1 - 2 * n * rho**2 - s.Rational(1, 2) - (64 * n - 1) / (128 * n)
            ),
            "component_bound_expansion": normalized
            - (64 + 384 / k + 6144 * n / k**2 + 9216 * n / k**3),
            "component_bound_at_minimum": upper_at_minimum - s.Rational(269833, 4096),
            "Euclidean_action_coefficient": euclidean
            - KAPPA * (256 * N / KAPPA**2) ** 2 / 2,
            "leading_vacuum_source_from_regular_r": leading.coeff(eps, 4)
            - n * gradient * (z - Xphi * box),
            "leading_source_lower_degrees_zero": sum(
                leading.coeff(eps, i) for i in range(4)
            ),
            "full_source_L2_four_components": 2 * (128 * N / KAPPA**2)
            - 256 * N / KAPPA**2,
        },
        "gates": {
            "n_even_at_least_four": N % 2 == 0 and N >= 4,
            "actual_kappa_in_small_domain": bool(KAPPA >= 64 * N),
            "component_normalized_uniform_bound": bool(upper_at_minimum < 128),
            "Euclidean_full_source_numerical_bound": bool(
                euclidean < s.Rational(1, 10**2389)
            ),
            "actual_parent_X_domain_included": bool(
                s.Rational(1, 16 * N) < s.Rational(1, 4096)
            ),
        },
    }
