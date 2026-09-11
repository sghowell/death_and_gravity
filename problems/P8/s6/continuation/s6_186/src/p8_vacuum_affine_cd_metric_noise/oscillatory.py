"""Full pair time smearing, exact Leibniz coefficients and radial majorants."""

from functools import cache

import sympy as s
from p8_vector_state.comparison import AMAX

from . import reference

REF_BOUND = s.Integer(10) ** 25


@cache
def data():
    t = s.Symbol("t", real=True)
    g = s.Function("inverse_phase")(t)
    a = s.Function("amplitude")(t)
    h = s.Function("test")(t)
    L = lambda b: s.diff(g * b, t)
    actual = s.expand(L(L(L(a * h))))
    expected = (
        g**3 * s.diff(a * h, t, 3)
        + 6 * g * g * s.diff(g, t) * s.diff(a * h, t, 2)
        + (7 * g * s.diff(g, t) ** 2 + 4 * g * g * s.diff(g, t, 2)) * s.diff(a * h, t)
        + (
            s.diff(g, t) ** 3
            + 4 * g * s.diff(g, t) * s.diff(g, t, 2)
            + g * g * s.diff(g, t, 3)
        )
        * a
        * h
    )
    r, C, G = s.symbols("radius amplitude_bound inverse_phase_bound", positive=True)
    derivatives = {s.diff(a, t, j): s.factorial(j) * C / r**j for j in range(4)}
    derivatives.update({s.diff(g, t, j): s.factorial(j) * G / r**j for j in range(4)})
    hsyms = s.symbols("h0:4")
    expanded = s.expand(
        actual.subs({s.diff(h, t, j): hsyms[j] for j in range(4)}, simultaneous=True)
    )
    constants = {
        j: s.simplify(
            expanded.coeff(hsyms[j]).subs(derivatives, simultaneous=True)
            / (C * G**3 / r ** (3 - j))
        )
        for j in range(4)
    }
    coefficient = 8 * sum(constants.values()) * reference.PAIR_REF / reference.RADIUS**3
    y = s.Symbol("y", positive=True)
    radial2 = s.integrate(y * y / (1 + y * y) ** 2, (y, 0, s.oo))
    radial5 = s.integrate(y * y / (1 + y * y) ** 5, (y, 0, s.oo))
    m = reference.MASS
    J4 = AMAX**3 / (8 * s.pi * m)
    J10 = 5 * AMAX**3 / (512 * s.pi * m**7)
    nu, mu = s.symbols("nu mu", positive=True)
    checks = {
        "complete_third_integration_by_parts_operator": s.expand(actual - expected),
        **{
            f"full_test_Leibniz_majorant_{j}": constants[j] - [48, 33, 9, 1][j]
            for j in range(4)
        },
        "complete_radial_resolvent_fourth": s.simplify(radial2 - s.pi / 4),
        "complete_radial_resolvent_tenth": s.simplify(radial5 - 5 * s.pi / 256),
        "full_angular_fourth_integral": s.simplify(
            AMAX**3 * radial2 / (2 * s.pi**2 * m) - J4
        ),
        "full_angular_tenth_integral": s.simplify(
            AMAX**3 * radial5 / (2 * s.pi**2 * m**7) - J10
        ),
        "reference_AM_GM_gap": s.expand((nu + mu) ** 2 - 4 * nu * mu - (nu - mu) ** 2),
        "squared_remainder_split_gap": s.factor(
            2 * (nu**-12 + mu**-12) - (nu**-6 + mu**-6) ** 2 - (nu**-6 - mu**-6) ** 2
        ),
    }
    return {
        "compact_time_test": "Arbitrary smooth compact symmetric physical-frame f_ab(t,x), supported strictly inside the unit CD slab; no band restriction",
        "inverse_sum_phase_holomorphic_upper": "2/(nu+mu), from Re W_k>nu/2 and Re W_l>mu/2",
        "third_Leibniz_coefficients": constants,
        "continuous_reference_pair_time_bound": coefficient,
        "reference_pair_display": REF_BOUND,
        "reference_pair_estimate": "|I_ref(k,l)|<=1e25 sqrt(nu mu)/(nu+mu)^3 U3(P), P=k+l",
        "error_pair_estimate": "|I_error(k,l)|<=1e15 sqrt(nu mu)(nu^-6+mu^-6) U0(P)",
        "fourth_resolvent_radial_integral": J4,
        "tenth_resolvent_radial_integral": J10,
        "reference_full_internal_integral_upper": J4 / 4,
        "error_full_internal_integral_upper": "4(1+|P|/(Amax m))*J10 <8 m^-7 (1+|P|^2)",
        "Fourier_convention": "Each spatial integral has d^3k/(2pi)^3; P=k+l preserves this inner/outer convention; Uj(P)^2=sum_(r<=j) integral_I ||partial_t^r fhat(t,P)||F^2 dt",
        "checks": checks,
        "gates": {
            "display_exceeds_complete_reference_pair_bound": coefficient < REF_BOUND,
            "reference_internal_integral_below_inverse_mass": J4 / 4 < 1 / m,
            "error_internal_integral_below_inverse_mass_seventh": J10 < 1 / m**7,
            "unit_time_Cauchy_weight": s.Rational(1, 2) - s.Rational(-1, 2) <= 1,
            "same_fixed_positive_mass": m == 1000,
            "Cauchy_radius_below_one": 0 < reference.RADIUS < 1,
        },
    }
