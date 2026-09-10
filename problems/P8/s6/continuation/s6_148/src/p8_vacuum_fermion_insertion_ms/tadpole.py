"""Regulated spectral tadpole, two exact beta anchors and its finite remainder."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_outer_ms import laurent
from p8_vacuum_fermion_self_energy_chord.tail import rational


def enclosure(m, Y, Q):
    m, Y, Q = map(rational, (m, Y, Q))
    if m < 2 or Y < 0 or not 0 < Q <= 144:
        raise ValueError("Need m>=2,Y>=0 and 0<Q_lower<=144")
    T = 4 * m * m
    pref = 12 * Y / Q**2
    return {
        "threshold": T,
        "C_over_Q_upper": pref,
        "leading_tadpole_absolute_upper": pref * s.Rational(21, 4) * T,
        "subleading_tadpole_absolute_upper": pref * 8,
        "finite_tadpole_remainder_absolute_upper": pref * 24 / T,
        "complete_tadpole_absolute_upper": pref * (s.Rational(21, 4) * T + 8 + 24 / T),
    }


@cache
def data():
    e = s.Symbol("epsilon", real=True)
    k = s.Symbol("k", positive=True)
    previous = laurent.data()
    F = previous["dimensionless_leading_regulated_F"]
    R = (s.Rational(3, 2) + e) / (1 - 2 * e)
    normal = s.simplify(e * e * R * F)
    reduce = lambda expr: s.simplify(
        expr.xreplace(
            {s.polygamma(1, s.Rational(5, 2)): s.pi**2 / 2 - s.Rational(40, 9)}
        )
    )
    jets = tuple(
        reduce(s.diff(normal, e, n).subs(e, 0) / s.factorial(n)) for n in range(3)
    )
    nextjets = tuple(
        reduce(s.diff(-2 * e * e * F, e, n).subs(e, 0) / s.factorial(n))
        for n in range(3)
    )
    pref_normal = (
        -s.exp(2 * s.EulerGamma * e)
        * 4 ** (-e)
        * s.sqrt(s.pi)
        * s.gamma(1 + e)
        / (2 * (1 - e) * s.gamma(s.Rational(3, 2) - e))
    )
    pf = reduce(s.diff(pref_normal, e).subs(e, 0))
    h2 = (1 - k) ** (-2) - 1 - 2 * k
    gamma_direct = (
        s.exp(2 * s.EulerGamma * e)
        * 4 ** (-e)
        * s.sqrt(s.pi)
        / 2
        * (s.Rational(3, 2) - e)
        * s.gamma(e - 1)
        * s.gamma(2 * e - 1)
        / s.gamma(s.Rational(3, 2) + e)
    )
    gamma_reduced = gamma_direct.xreplace(
        {
            s.gamma(e - 1): s.gamma(1 + e) / (e * (e - 1)),
            s.gamma(2 * e - 1): s.gamma(1 + 2 * e) / (2 * e * (2 * e - 1)),
            s.gamma(s.Rational(3, 2) + e): s.gamma(s.Rational(5, 2) + e)
            / (s.Rational(3, 2) + e),
        }
    )
    ffinite = previous["leading_finite_MS_coefficient"]
    return {
        "exact_leading_tadpole_per_T_over_C_over_Q": R * F,
        "exact_subleading_tadpole_over_C_over_Q": -2 * F,
        "leading_double_pole_simple_pole_finite": jets,
        "subleading_double_pole_simple_pole_finite": nextjets,
        "finite_remainder_prefactor": pf,
        "finite_tadpole_remainder": "T integral_0^1 (1-z)^(3/2)/z^2 h2(z/T) [4log2-3-2logz+log(1-z)] dz, h2(k)=(1-k)^(-2)-1-2k.",
        "remainder_pole": "-T integral_0^1 (1-z)^(3/2)/z^2 h2(z/T) dz / epsilon; nonzero, removed by MS before retaining its finite product.",
        "complete_finite_tadpole": "T_insert,MS=(C/Q){T(13/4+pi^2/8)-47/9-pi^2/6+delta_T}, |delta_T|<=24/T, C=2NY/Q,T=4m^2.",
        "checks": {
            "leading_beta_Gamma_normalization": s.factor(gamma_reduced - R * F),
            "leading_double_pole": jets[0] - s.Rational(3, 4),
            "leading_simple_pole": jets[1] - s.Rational(1, 4),
            "leading_finite": jets[2] - s.Rational(13, 4) - s.pi**2 / 8,
            "subleading_double_pole": nextjets[0] + 1,
            "subleading_simple_pole": nextjets[1] - s.Rational(7, 3),
            "subleading_finite": nextjets[2] + s.Rational(47, 9) + s.pi**2 / 6,
            "subleading_is_minus_twice_prior_F_finite": nextjets[2] + 2 * ffinite,
            "remainder_residue": pref_normal.subs(e, 0) + 1,
            "remainder_finite_prefactor": pf - 4 * s.log(2) + 3,
            "exact_remainder_two_power_factor": s.factor(
                h2 - k * k * (3 - 2 * k) / (1 - k) ** 2
            ),
            "remainder_below_four_k_squared": s.factor(
                4 * k * k - h2 - k * k * (1 - 6 * k + 4 * k * k) / (1 - k) ** 2
            ),
            "remainder_positive_gap_on_domain": 1
            - 6 * k
            + 4 * k * k
            - (s.Rational(5, 8) + 6 * (s.Rational(1, 16) - k) + 4 * k * k),
            "finite_remainder_constant": 4 * (3 + 2 + 1) - 24,
            "leading_finite_pi_upper": s.Rational(13, 4)
            + 16 / s.Integer(8)
            - s.Rational(21, 4),
            "subleading_finite_pi_upper_gap": 8
            - (s.Rational(47, 9) + 16 / s.Integer(6))
            - s.Rational(1, 9),
            "zero_Y_tadpole": enclosure(2, 0, 144)["complete_tadpole_absolute_upper"],
        },
    }
