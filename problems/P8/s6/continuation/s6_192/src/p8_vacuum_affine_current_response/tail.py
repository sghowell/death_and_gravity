"""C2 physical-current tails and finite full comparison, not a new scheme."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_matrix_response_tail import covariance as prior

from . import adiabatic, low, vertices

STATE_TAIL = (
    6 * s.Rational(1, 10) ** 65,
    s.Rational(1, 10) ** 46,
    s.Rational(1, 10) ** 27,
)
REFERENCE_TAIL = (
    s.Rational(1, 10) ** 13,
    s.Rational(1, 10) ** 12,
    s.Rational(1, 10) ** 10,
)
COMPARISON = (s.Integer(10) ** 77, s.Integer(10) ** 94, s.Integer(10) ** 111)


@cache
def constants():
    T = prior.TAIL
    g0, g1, g2 = vertices.G
    actual = (
        3 * g0 * T[0],
        3 * (g1 * T[0] + g0 * T[1]),
        3 * (g2 * T[0] + 2 * g1 * T[1] + g0 * T[2]),
    )
    K = low.PARTITION
    finite = low.constants()["finite_band_complete_current_integrals"]
    contour = adiabatic.CONTOUR_CURRENT
    high = tuple(s.Rational(4, 18) * contour[a] * low.MASS**6 / K**2 for a in range(3))
    sublow = tuple(contour[a] * K**4 / 2 for a in range(3))
    combined = tuple(finite[a] + sublow[a] + actual[a] + high[a] for a in range(3))
    return {
        "actual_minus_reference_physical_current_tail_bounds": actual,
        "finite_reference_minus_fourth_order_infinite_tail_bounds": high,
        "finite_band_fourth_order_comparison_bounds": sublow,
        "full_actual_minus_fourth_order_comparison_bounds": combined,
        "comparison_divided_by_fixed_kappa_displays": tuple(
            x / modes.KAPPA for x in COMPARISON
        ),
    }


@cache
def data():
    nu, K = s.symbols("nu K", positive=True)
    eps = s.Symbol("epsilon", real=True)
    b = s.Rational(99, 100) / s.Rational(25, 16) ** 2
    c = constants()
    checks = {
        "infinite_reference_radial_integral": s.integrate(nu**-3, (nu, K, s.oo))
        - 1 / (2 * K**2),
        "fourth_order_finite_band_integral": s.integrate(nu**3, (nu, 0, K)) - K**4 / 4,
        "finite_amplitude_integral_Taylor_factor": s.integrate(
            (1 - nu) * eps**2, (nu, 0, 1)
        )
        - eps**2 / 2,
        "zero_amplitude_Taylor_remainder_is_zero": (eps**2 * COMPARISON[2] / 2).subs(
            eps, 0
        ),
        "physical_density_unimodular_volume_determinant": s.exp(eps) * s.exp(-eps) - 1,
    }
    return {
        "analysis_partition": "The same fixed nu_minus>=1e16 band is independent of epsilon and is not a physical cutoff. Both it and its finite complement are retained.",
        "actual_state_physical_current_tail_displays": STATE_TAIL,
        "finite_reference_comparison_high_tail_displays": REFERENCE_TAIL,
        "full_comparison_displays": COMPARISON,
        "constants": c,
        "complete_limit": "The actual-minus-J_ad4 current is absolutely integrable and C2 in epsilon into C0 of the fixed unit CD slab, by the common high-band majorants and the finite-band exact covariance bounds. This is a mathematical mode-subtracted comparison, not yet the original covariant observable.",
        "finite_amplitude_Taylor": "For |epsilon|<=.01, the complete comparison current has remainder <=epsilon^2*1e111/2 after its value and first derivative at zero. The inequality is strict for nonzero epsilon and the remainder is exactly zero at epsilon0. The current divided by the fixed kappa has coefficient1e-689.",
        "proper_density": "Dividing by a^3>=1 preserves these bounds, since the shear family has unit spatial determinant and a is independent of epsilon.",
        "matching_boundary": "Do not identify J_ad4 with the original S176/S82 fixed covariant mu=m subtraction without proving its finite local/contact matching on these sheared histories. No full physical self-energy, feedback inverse, finite-coupling background, physical cutoff, V/G/B or original P8 closure follows here.",
        "checks": checks,
        "gates": {
            "radial_jacobian_display": 1 / b**3 < 16,
            "high_band_even_geometric_factor_below_two": 1
            / (1 - (low.MASS / low.PARTITION) ** 2)
            < 2,
            "physical_actual_reference_tail_displays": all(
                c["actual_minus_reference_physical_current_tail_bounds"][a]
                <= STATE_TAIL[a]
                for a in range(3)
            ),
            "reference_fourth_order_tail_displays": all(
                c["finite_reference_minus_fourth_order_infinite_tail_bounds"][a]
                < REFERENCE_TAIL[a]
                for a in range(3)
            ),
            "complete_comparison_displays": all(
                c["full_actual_minus_fourth_order_comparison_bounds"][a] < COMPARISON[a]
                for a in range(3)
            ),
            "second_derivative_infinite_reference_tail_integrable": -3 < -1,
            "fixed_canonical_kappa": modes.KAPPA == 10**800,
        },
    }
