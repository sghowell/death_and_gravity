"""Complete one-mode local contact occupation correction."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import reference
from p8_vector_state.comparison import AMAX

from . import prefactor


@cache
def constants():
    m = reference.MASS
    coefficient = s.Rational(3, 2) * reference.FIELD_NORM**2 * prefactor.BETA**2
    J11 = 8 * AMAX**3 / (2835 * m**8)
    return {
        "complete_three_mode_contact_coefficient": coefficient,
        "J11_pi_lower_bound": J11,
        "contact_L2_bound": coefficient * J11,
        "full_band_K8_tail_numerator": coefficient * AMAX**11 / 144,
        "full_band_K_tail_numerator": coefficient * AMAX**11 / (144 * m**7),
    }


@cache
def data():
    b = s.Symbol("occupation", nonnegative=True)
    v = s.Matrix(s.symbols("v0:10", complex=True))
    a = s.sqrt(1 + b)
    D = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]])
    G = s.Matrix([[2, -1, 1], [-1, -3, 2], [1, 2, 1]]) / 7
    e, h = s.symbols("epsilon eta", real=True)
    X = e * D + h * G
    H = (D * G + G * D) / 2
    t, K = s.symbols("t K", positive=True)
    c = constants()
    checks = {
        "complete_ten_field_covariance_occupation_correction": a
        * v
        * (a * v).conjugate().T
        - v * v.conjugate().T
        - b * v * v.conjugate().T,
        "positive_full_noncommuting_contact": (s.eye(3) + X + X * X / 2).applyfunc(
            lambda x: s.expand(x).coeff(e, 1).coeff(h, 1)
        )
        - H,
        "full_three_mode_half_Hamiltonian_factor": c[
            "complete_three_mode_contact_coefficient"
        ]
        - s.Rational(3, 2) * 4000**2 * (2 * 10**6) ** 2,
        "complete_eleventh_massive_radial_integral": s.integrate(
            t * t / (1 + t * t) ** s.Rational(11, 2), (t, 0, s.oo)
        )
        - s.Rational(16, 315),
        "full_band_eleventh_radial_tail": s.integrate(t**-9, (t, K, s.oo))
        - 1 / (8 * K**8),
    }
    return {
        "definition": "R_alpha,contact is the full local second-metric contact with alpha_initial W8 minus the same contact with unit W8. Each one-mode covariance differs exactly by b_k v_k v_k^dagger.",
        "full_contact": "The noncommuting product H=(D Gamma+Gamma D)/2 has positive mass sign; in Hamiltonian feature order the full matrix is diag(H,H,H,0). Its norm is<=||D||F||Gamma||F, including every physical mode and the temporal constraint.",
        "pointwise": "The complete half-Hamiltonian sum over all3 modes is bounded by(3/2)4000^2 B^2 nu^-11 ||D||F||Gamma||F. No reference contact or longitudinal contribution is discarded.",
        "integral": "J11=integral d^3k/(2pi)^3 nu^-11=8 Amax^3/(315pi^2 m^8). Thus |R_alpha,contact|<1e-3||D||L2||Gamma||L2 on the fixed slab.",
        "tail": "The original ONE-mode removed band gives J11_tail(K)<=Amax^11/(16pi^2 K^8)<Amax^11/(144K^8). For K>=m its contact tail is below||D||L2||Gamma||L2/K. This is not the memory overlap region.",
        "constants": c,
        "checks": checks,
        "gates": {
            "full_contact_strict_display": c["contact_L2_bound"] < s.Rational(1, 1000),
            "full_contact_K_tail_display": c["full_band_K_tail_numerator"] < 1,
            "J11_below_inverse_mass_eighth": 0
            < c["J11_pi_lower_bound"]
            < 1 / reference.MASS**8,
            "noncommuting_contact_retained": D * G != G * D,
        },
    }
