"""Actual mass-one scalar vacuum with the whole inner on-shell forest retained."""

from functools import cache

import sympy as s
from p8_vacuum_fermion_insertion_ms import tadpole
from p8_vacuum_fermion_outer_ms import laurent
from p8_vacuum_fermion_self_energy_chord.tail import rational


def enclosure(m, Y, Q):
    m, Y, Q = map(rational, (m, Y, Q))
    if m < 2 or Y < 0 or not 0 < Q <= 144:
        raise ValueError("Need m>=2,Y>=0 and 0<Q_lower<=144")
    T = 4 * m * m
    pref = 6 * Y / Q**2
    return {
        "threshold": T,
        "massless_scalar_vacuum_absolute_upper": pref * 19 * m**4,
        "actual_scalar_vacuum_mass_correction_absolute_upper": pref
        * (s.Rational(21, 4) * T + 4 + 12 / T),
        "finite_scalar_vacuum_remainder_absolute_upper": pref * 12 / T,
        "complete_scalar_vacuum_absolute_upper": pref
        * (19 * m**4 + s.Rational(21, 4) * T + 4 + 12 / T),
    }


@cache
def data():
    e = s.Symbol("epsilon", real=True)
    k = s.Symbol("k", positive=True)
    F = laurent.data()["dimensionless_leading_regulated_F"]
    A = s.exp(s.EulerGamma * e) * s.gamma(1 + e)
    reduced = F.xreplace(
        {
            s.gamma(1 + 2 * e): 2 ** (2 * e)
            * s.gamma(s.Rational(1, 2) + e)
            * s.gamma(1 + e)
            / s.sqrt(s.pi),
            s.gamma(s.Rational(5, 2) + e): (e + s.Rational(1, 2))
            * (e + s.Rational(3, 2))
            * s.gamma(e + s.Rational(1, 2)),
            4 ** (-e): 2 ** (-2 * e),
        }
    )
    Frat = A * A / e**2 * (3 - 2 * e) / (2 * (1 - e) * (1 + 2 * e) * (3 + 2 * e))
    K0 = (
        -F
        * (e + s.Rational(1, 2))
        * (e + s.Rational(3, 2))
        / ((2 * e - 2) * (2 * e - 1))
    )
    vs = (3 - 2 * e) / ((1 - e) ** 2 * (2 * e - 1))
    h3 = 1 / (1 - k) - 1 - k - k * k
    t = tadpole.data()
    full_first = (s.Rational(3, 2) + e) / (1 - 2 * e) * F
    x, v = s.symbols("q_squared v", positive=True)
    return {
        "complete_scalar_OS_vacuum_spectral_identity": "V_OS,D=1/2 integral rho_D(v) Tad_D(v)/(v-1) dv. The polynomial trace is scaleless; the remaining integral is continued from 1<Re(epsilon)<3/2.",
        "regulated_scalar_OS_vacuum": "V_OS,D/(NY/Q^2)=T^2 Gamma(epsilon-1)H(epsilon) integral_0^1 z^(2epsilon-3)(1-z)^(3/2-epsilon)/(1-z/T) dz, H=e^(2gamma epsilon)4^(-epsilon)sqrt(pi)/(2Gamma(3/2-epsilon)).",
        "leading_beta_reference_per_T_squared": K0,
        "second_beta_reference_per_T": full_first,
        "third_beta_reference": -F,
        "finite_actual_scalar_reference": "V_scalar,paired=(NY/Q^2){-19m^4+T(13/4+pi^2/8)-47/18-pi^2/12+delta_V}. Includes the proper fermion MS counterterm and the full physical scalar mass/residue forest.",
        "finite_remainder": "delta_V=T^2 integral_0^1 (1-z)^(3/2)/z^3 h3(z/T)[4log2-3-2logz+log(1-z)] dz, h3(k)=k^3/(1-k); |delta_V|<=12/T.",
        "remainder_pole": "The compact h3 integral has a nonzero value at epsilon=0. Multiplying by Gamma(epsilon-1)H gives a nonzero simple pole and a finite prefactor product; subtract the pole only.",
        "checks": {
            "explicit_Gamma_duplication_reduction": s.simplify(reduced - Frat),
            "leading_beta_is_full_raw_massless_vacuum": s.factor(
                -16
                * Frat
                * (e + s.Rational(1, 2))
                * (e + s.Rational(3, 2))
                / ((2 * e - 2) * (2 * e - 1))
                - A * A * vs / e**2
            ),
            "second_beta_is_previous_tadpole_leading": s.simplify(
                full_first - t["exact_leading_tadpole_per_T_over_C_over_Q"]
            ),
            "third_beta_is_half_previous_subleading": s.simplify(
                -F - t["exact_subleading_tadpole_over_C_over_Q"] / 2
            ),
            "three_terms_and_remainder_exact": s.factor(h3 - k**3 / (1 - k)),
            "remainder_below_two_k_cubed": s.factor(
                2 * k**3 - h3 - k**3 * (1 - 2 * k) / (1 - k)
            ),
            "positive_gap_for_remainder_majorant": 1
            - 2 * k
            - (s.Rational(7, 8) + 2 * (s.Rational(1, 16) - k)),
            "complete_remainder_log_moment_constant": 2 * (3 + 2 + 1) - 12,
            "spectral_contact_decomposition": s.factor(
                (x + 1) / (x + v) - 1 - (1 - v) / (x + v)
            ),
            "vacuum_half_trace_after_contact_removal": -s.Rational(1, 2)
            * (1 - v)
            / (v - 1) ** 2
            - s.Rational(1, 2) / (v - 1),
            "zero_Y_scalar_vacuum": enclosure(2, 0, 144)[
                "complete_scalar_vacuum_absolute_upper"
            ],
        },
    }
