"""Exact positive-series enclosures without large-argument log cancellation."""

from fractions import Fraction

import sympy as sp
from p8_polynomial_vacuum import model

from . import moments


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require a finite exact rational")
    return sp.Rational(value)


def q2_enclosure(argument, terms=2):
    Z, N = rational(argument), rational(terms)
    if Z <= 1:
        raise ValueError("Require real Legendre argument strictly above one")
    if N.q != 1 or not 1 <= N <= 32:
        raise ValueError("Require one through 32 integer series terms")
    count = int(N)
    lower = sum(
        sp.Rational(2 * n, (2 * n + 1) * (2 * n + 3)) / Z ** (2 * n + 1)
        for n in range(1, count + 1)
    )
    tail = sp.Rational(2, 15) / (Z ** (2 * count + 3) * (1 - 1 / Z**2))
    return {
        "argument": Z,
        "retained_terms": count,
        "strict_lower": lower,
        "strict_upper": lower + tail,
        "tail_upper": tail,
    }


def point(transfer=5):
    T = rational(transfer)
    if not 4 <= T <= 6:
        raise ValueError("Require transfer invariant in the proved low window [4,6]")
    p = model.data()["actual_parameters"]
    M = p["heavy_mass_squared"]
    g = p["cubic_coupling_squared"]
    if T == 4:
        return {
            "physical_transfer": T,
            "light_cut_lower": sp.S.Zero,
            "light_cut_upper": sp.S.Zero,
            "heavy_cut_exact": sp.S.Zero,
            "threshold": True,
        }
    q = q2_enclosure(1 + 2 * M / (T - 4))
    pref = g / (8 * sp.sqrt(T * (T - 4)))
    return {
        "physical_transfer": T,
        "light_cut_lower": pref * q["strict_lower"] / 4,
        "light_cut_upper": pref * q["strict_upper"] / 3,
        "heavy_cut_exact": sp.S.Zero,
        "threshold": False,
        "exact_Q2_enclosure": q,
        "scope": "Actual low-window imaginary spin-two vertex; pi is enclosed by 3<pi<4. No amplitude or Regge coefficient is substituted for it.",
    }


def data():
    q = q2_enclosure(2, 1)
    d = point()
    m = moments.data()
    return {
        "actual_selected_transfer": d,
        "checks": {
            "one_term_positive_Q2_lower": q["strict_lower"] - sp.Rational(2, 15) / 2**3,
            "one_term_Q2_geometric_tail": q["tail_upper"]
            - sp.Rational(2, 15) / (2**5 * (1 - sp.Rational(1, 4))),
            "threshold_imaginary_vertex_zero": point(4)["light_cut_upper"],
            "heavy_cut_absent_on_low_window": d["heavy_cut_exact"],
        },
        "bounds": {
            "actual_point_strictly_positive_light_cut_lower": d["light_cut_lower"] > 0,
            "actual_point_lower_below_upper": d["light_cut_lower"]
            < d["light_cut_upper"],
            "actual_point_Q2_bound_within_uniform_density_bound": d["light_cut_upper"]
            < model.data()["actual_parameters"]["cubic_coupling_squared"]
            / (480 * m["actual_heavy_mass_squared"] ** 3),
            "actual_spectral_fraction_bound_nonzero": m[
                "actual_low_window_fraction_upper"
            ]
            > 0,
        },
    }


def bad_cases():
    invalid = (
        True,
        False,
        1.0,
        sp.Float(1),
        "1",
        sp.I,
        sp.oo,
        -sp.oo,
        sp.zoo,
        sp.nan,
        None,
        sp.Symbol("p"),
    )
    rows = []
    for i, v in enumerate(invalid):
        rows += [
            ("inexact_transfer_" + str(i), point, (v,)),
            ("inexact_argument_" + str(i), q2_enclosure, (v,)),
            ("inexact_term_count_" + str(i), q2_enclosure, (2, v)),
        ]
    rows += [
        ("outside_cut_window_" + str(i), point, (v,))
        for i, v in enumerate((0, -1, 3, 7))
    ]
    rows += [
        ("outside_Q2_domain_" + str(i), q2_enclosure, (v,))
        for i, v in enumerate((0, -1, 1))
    ]
    rows += [
        ("unsupported_series_count_" + str(i), q2_enclosure, (2, v))
        for i, v in enumerate((0, -1, 33, Fraction(3, 2)))
    ]
    return rows
