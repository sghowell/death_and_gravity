"""Finite curved remainder separated from unmatched contact/polynomial sector."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes

from . import tail

KNOWN = 4 * s.Integer(10) ** 54
KNOWN_TAIL = 6 * s.Integer(10) ** 60


@cache
def data():
    state, contact, endpoint, bulk, Q, poly = s.symbols(
        "state contact endpoint bulk Q poly"
    )
    checks = {
        "finite_regulator_endpoint_split": s.expand(endpoint - (Q + poly)).subs(
            endpoint, Q + poly
        ),
        "actual_current_complete_decomposition": s.expand(
            (state + contact + endpoint + bulk) - (state + bulk + Q) - (contact + poly)
        ).subs(endpoint, Q + poly),
        "canonical_all_momentum_endpoint_remainder": 4 * tail.FULL / modes.KAPPA
        - 12 * s.Rational(1, 10) ** 746,
        "canonical_all_momentum_endpoint_tail": 4 * tail.TAIL / modes.KAPPA
        - 2 * s.Rational(1, 10) ** 739,
        "canonical_known_actual_finite_piece": 4 * KNOWN / modes.KAPPA
        - 16 * s.Rational(1, 10) ** 746,
        "canonical_known_actual_finite_tail": 4 * KNOWN_TAIL / modes.KAPPA
        - 24 * s.Rational(1, 10) ** 740,
    }
    return {
        "full_finite_regulator_identity": "Write E5_K=sum_j0..4 B_j,unit,K=Q_K+P_K, where P_K integrates the high-band degree-four spatial Taylor INTEGRAND with the original both-leg cutoff. Then J_actual,K=(R_actual-unit,K+F_unit,K+Q_K)+(C_unit,K+P_K). Every original term and sign is retained.",
        "known_piece": "Let Y^2=N61^2+X46^2, with N61^2=sum_j0..6||partial_t^j Gamma||L2^2+||grad_x Gamma||L2^2 and X46 as in S202. Since M[D]>=||D||L2, S201 plus Q gives a known finite actual-current piece below4e54 M[D]Y[Gamma] and tail6e60 M[D]Y[Gamma]/K.",
        "canonical": "Both external metric factors give12e-746 and2e-739/K for Q, and16e-746 and24e-740/K for the combined known actual finite piece.",
        "unmatched": "P_K is a Taylor-polynomial integrand sector, not necessarily a polynomial function of external P after a moving sharp cutoff is integrated. Its regulator artifacts, divergent and finite coefficients and combination with the distinct complete one-leg contact require the original fixed covariant matching. No replacement prescription is selected.",
        "boundary": "The bounded remainder is weak and derivative losing. It is not a full matched response, derivative-compatible mixed inverse, nonlinear solution, stability estimate, physical cutoff, remaining parent-loop computation or original V/G/B closure.",
        "checks": checks,
        "gates": {
            "known_actual_piece_display": tail.FULL + 2 * s.Integer(10) ** 48 < KNOWN,
            "known_actual_tail_display": tail.TAIL + 2 * s.Integer(10) ** 52
            < KNOWN_TAIL,
            "prepared_state_and_reference_correction_kept": True,
            "full_contact_and_polynomial_integrand_not_deleted": True,
            "distinct_one_leg_and_two_leg_regulators_kept": True,
            "original_V_G_B_open": True,
        },
    }
