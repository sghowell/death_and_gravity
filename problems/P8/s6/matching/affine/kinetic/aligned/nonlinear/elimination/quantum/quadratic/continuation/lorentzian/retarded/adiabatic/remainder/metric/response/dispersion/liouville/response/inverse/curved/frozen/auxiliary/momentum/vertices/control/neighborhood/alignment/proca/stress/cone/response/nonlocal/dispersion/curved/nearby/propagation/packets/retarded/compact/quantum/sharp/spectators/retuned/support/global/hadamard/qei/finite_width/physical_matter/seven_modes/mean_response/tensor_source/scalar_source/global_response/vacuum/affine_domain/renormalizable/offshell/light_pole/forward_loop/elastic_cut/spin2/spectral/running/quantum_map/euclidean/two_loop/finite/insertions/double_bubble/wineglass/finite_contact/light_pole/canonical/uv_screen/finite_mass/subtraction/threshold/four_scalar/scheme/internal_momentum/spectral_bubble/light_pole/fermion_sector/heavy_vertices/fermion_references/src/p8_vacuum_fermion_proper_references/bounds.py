"""Exact conservative finite-reference bounds without numerical logarithms."""

from functools import cache

import sympy as sp


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, sp.Integer, sp.Rational)):
        raise TypeError("Require an exact finite rational")
    return sp.Rational(value)


def enclosure(m, Yhi, ahi, Phi_abs, Qlo):
    m, Yhi, ahi, Phi_abs, Qlo = map(rational, (m, Yhi, ahi, Phi_abs, Qlo))
    if m < 2 or Yhi < 0 or ahi < 0 or not 0 <= Phi_abs <= sp.Rational(1, 2) or Qlo <= 0:
        raise ValueError(
            "Need m>=2, nonnegative couplings, |w| bound at most one half, Qlo>0"
        )
    Cf = sp.Rational(4, 3)
    zhi = (3 * Yhi / 4 + ahi * Cf / 2) / Qlo
    etaabs = (Yhi + 2 * ahi * Cf) / Qlo
    upsabs = (2 * Yhi + 6 * ahi * Cf) / Qlo
    return {
        "scalar_boson_to_fermion_mass_squared_ratio": 1 / (m * m),
        "fermion_kinetic_interval": [1, 1 + zhi],
        "finite_fermion_kinetic_increment_upper": zhi,
        "finite_mass_relative_anchor_absolute_upper": etaabs,
        "finite_Yukawa_relative_anchor_absolute_upper": upsabs,
        "complete_Phi_kinetic_increment_absolute_upper": Phi_abs,
        "formal_mass_conversion_absolute_upper": etaabs + zhi,
        "formal_Yukawa_conversion_absolute_upper": upsabs + zhi + Phi_abs / 2,
        "selected_mass_ratio_minus_one_absolute_upper": etaabs + zhi,
        "selected_Yukawa_ratio_minus_one_absolute_upper": 2 * (upsabs + zhi + Phi_abs),
        "positive_selected_mass_reference": bool(etaabs < 1),
        "positive_selected_Yukawa_reference": bool(upsabs < 1),
        "scope": "Selected one-loop local-reference dictionary, not exact charged-particle poles or completed two-loop matching.",
    }


@cache
def data():
    j0, R = sp.symbols("J0 R", real=True)
    z, eta, ups, w, t = sp.symbols("z eta upsilon w sqrt_one_plus_w", real=True)
    mass = (1 + eta) / (1 + z)
    ratio = (1 + ups) / ((1 + z) * t)
    m = sp.Symbol("m", positive=True)
    return {
        "anchor_box": "J0 in [-1,0], J1 in [-3/4,0], R in [1/2,1], hence J0+2R in [0,2]",
        "rational_scalar_field_inverse_sqrt_bound": "For |w|<=1/2, t=sqrt(1+w)>1/2, 1/t<2 and |1/t-1|<=2|w|.",
        "checks": {
            "mass_ratio_difference": sp.factor(mass - 1 - (eta - z) / (1 + z)),
            "Yukawa_ratio_difference_decomposition": sp.factor(
                ratio - 1 - ((ups - z) / ((1 + z) * t) + (1 / t - 1))
            ),
            "inverse_square_root_difference_identity": sp.factor(
                (1 / t - 1) + (t * t - 1) / (t * (t + 1))
            ),
            "positive_square_root_quarter_gap": (1 + w)
            - sp.Rational(1, 4)
            - (sp.Rational(3, 4) + w),
            "sqrt_denominator_above_three_quarters": sp.expand(
                t * (t + 1)
                - sp.Rational(3, 4)
                - (t - sp.Rational(1, 2)) * (t + sp.Rational(3, 2))
            ),
            "positive_fermion_inverse_upper": sp.factor(1 - 1 / (1 + z) - z / (1 + z)),
            "scalar_mass_ratio_within_quarter": sp.factor(
                sp.Rational(1, 4) - 1 / m**2 - (m * m - 4) / (4 * m * m)
            ),
            "scalar_Yukawa_anchor_lower_from_box": j0
            + 2 * R
            - (j0 + 1)
            - 2 * (R - sp.Rational(1, 2)),
            "scalar_Yukawa_anchor_upper_from_box": 2 - j0 - 2 * R - (-j0 + 2 * (1 - R)),
            "no_couplings_no_fermion_reference_increment": enclosure(2, 0, 0, 0, 144)[
                "selected_mass_ratio_minus_one_absolute_upper"
            ],
        },
    }
