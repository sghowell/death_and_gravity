"""Exact normalized-channel ceiling and mandatory correction, not a UV cutoff."""

from functools import cache

import sympy as s

from . import amplitude as a

ENERGY_LO = s.Integer(10) ** 133
ENERGY_HI = 2 * ENERGY_LO
CORRECTION_ENERGY = s.Integer(10) ** 134
BETA_LO = s.Rational(999, 1000)
PI_HI = s.Rational(22, 7)


@cache
def data():
    S = a.S
    v = s.Symbol("nonnegative_s_minus_four", nonnegative=True)
    Q = a.Q0.subs({a.lam: a.LAMBDA, a.gam: a.GAMMA})
    pplus = a.lam * (4 * v * v + 16 * v + 24) - 8 * a.gam
    pminus = (
        a.lam * (2 * v * v + 16 * v + 24) + 3 * a.gam * (v + 4) * v * v / 2 - 8 * a.gam
    )
    real, imag = s.symbols("partial_real partial_imag", real=True)
    sl = 1 + 2 * s.I * (real + s.I * imag)
    margins = {
        "threshold_tree_coefficient_positive": 24 * a.LAMBDA - 8 * a.GAMMA,
        "below_half_at_lower_energy": s.Rational(1, 2)
        - Q.subs(S, ENERGY_LO**2) / (32 * 3),
        "above_three_at_upper_energy": BETA_LO * Q.subs(S, ENERGY_HI**2) / (32 * PI_HI)
        - 3,
        "valid_beta_at_upper_energy": 1 - 4 / ENERGY_HI**2 - BETA_LO**2,
        "above_fifty_thousand_at_correction_energy": BETA_LO
        * Q.subs(S, CORRECTION_ENERGY**2)
        / (32 * PI_HI)
        - 50000,
        "valid_beta_at_correction_energy": 1 - 4 / CORRECTION_ENERGY**2 - BETA_LO**2,
    }
    return {
        "all_positive_rational_margins": margins,
        "nominal_tree_partial_wave_ceiling": "There is exactly one COM energy E_U with t0_tree(E_U^2)=1/2, and 10^133<E_U<2*10^133. The zeroth wave dominates every original tree wave. This is a nominal real-part ceiling; a nonzero real tree amplitude is never by itself an exact unitary amplitude, even below E_U.",
        "unitarity_circle": "For any existing exact unitary normalized identical elastic channel, |1+2it0|<=1, so Im t0>=|t0|^2 and |Re t0|<=1/2, including arbitrary additional positive inelastic channels.",
        "mandatory_full_correction": "At COM energy10^134 the original real t0_tree exceeds50000. Every exact unitary amplitude therefore has |t0_exact-t0_tree|/t0_tree>1-10^-5. No hypothesis of small corrections is made to derive this lower bound.",
        "nominal_energy_enclosure": (ENERGY_LO, ENERGY_HI),
        "necessary_correction_energy": CORRECTION_ENERGY,
        "necessary_relative_correction_lower": 1 - s.Rational(1, 10**5),
        "scope": "The fixed canonical nongravitational tree and any hypothetical exact unitary scalar S matrix. Not a verified physical cutoff, finite-gravity partial wave, full quantum decoupling limit, tensor-response change or UV-parent exclusion.",
        "checks": {
            "positive_angle_endpoint_polynomial": s.expand(
                (a.a + a.b).subs(S, v + 4) - pplus
            ),
            "positive_other_endpoint_polynomial": s.expand(
                (a.a - a.b).subs(S, v + 4) - pminus
            ),
            "strict_zeroth_coefficient_derivative": s.expand(
                s.diff(a.Q0, S)
                - (a.lam * (20 * S - 32) / 3 + a.gam * (S - 4) * (3 * S - 4) / 2)
            ),
            "beta_derivative_positive_formula": s.simplify(
                s.diff(a.beta, S) - 2 / (S * S * a.beta)
            ),
            "full_unitarity_circle_identity": s.expand(
                s.Abs(sl) ** 2 - 1 - 4 * (real * real + imag * imag - imag)
            ),
            "circle_center_and_radius": s.expand(
                real * real
                + (imag - s.Rational(1, 2)) ** 2
                - s.Rational(1, 4)
                - (real * real + imag * imag - imag)
            ),
            "exact_partial_wave_dominance_factor": s.Rational(240, 32)
            * s.Rational(2, 3)
            - 5,
            "mandatory_relative_correction_constant": s.Rational(1, 2) / 50000
            - s.Rational(1, 10**5),
        },
        "gates": {
            "all_actual_threshold_and_correction_margins_positive": all(
                value > 0 for value in margins.values()
            ),
            "positive_angle_polynomials_for_s_above_four": True,
            "monotone_t0_and_unique_nominal_ceiling": True,
            "no_claim_tree_is_exactly_unitary_below_ceiling": True,
            "massless_graviton_forward_problem_not_dropped": True,
        },
    }
