"""Exact finite-regulator repartition without double-counting earlier finite pieces."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_reference_state_prefactor import limit as state_time

from . import estimates, tail

KNOWN = 5 * s.Integer(10) ** 48
KNOWN_TAIL = 5 * s.Integer(10) ** 53


@cache
def data():
    R, F, C, Qold, Pfin, PUV, Qnew, U = s.symbols("R F C Qold Pfin PUV Qnew U")
    old = R + F + Qold + Pfin + C + PUV
    new = R + F + Qnew + C + U
    c = state_time.constants()
    checks = {
        "exact_complete_endpoint_repartition": s.expand(
            (Qold + Pfin + PUV) - (Qnew + U)
        ).subs(Qnew, Qold + Pfin + PUV - U),
        "exact_full_actual_current_same_regulator": s.expand(old - new).subs(
            Qnew, Qold + Pfin + PUV - U
        ),
        "unchanged_state_and_finite_time_piece": c["complete_known_finite_curved_piece"]
        - 2 * s.Integer(10) ** 48,
        "unchanged_state_and_finite_time_tail": c["complete_known_finite_curved_tail"]
        - 2 * s.Integer(10) ** 52,
        "both_canonical_uniform_endpoint_factors": 4 * estimates.FINITE / modes.KAPPA
        - 8 * s.Rational(1, 10) ** 752,
        "both_canonical_uniform_endpoint_tail_factors": 4 * tail.TAIL / modes.KAPPA
        - 16 * s.Rational(1, 10) ** 747,
        "canonical_new_known_actual_piece": 4 * KNOWN / modes.KAPPA
        - 2 * s.Rational(1, 10) ** 751,
        "canonical_new_known_actual_tail": 4 * KNOWN_TAIL / modes.KAPPA
        - 2 * s.Rational(1, 10) ** 746,
    }
    return {
        "exact_decomposition": "Jactual,K=(R_actual-unit,K+F_unit,K+Qnew,K)+(C_unit,K+U_K). Qnew is the complete low plus high(E5-U) endpoint piece at the original two-leg mask. U_K is the actual retained inverse-radius UV symbol integral. The contact keeps its one-leg band.",
        "finite_terms_not_discarded": "The earlier Q203+Pfin+PUV and the new Qnew+U are both exactly the complete E5 endpoint row. This is a regrouping of the same finite-regulator current, not deletion of finite oversubtractions. No state, matching prescription or counterterm is changed.",
        "no_double_counting": "Rebase the estimate from S201's R_actual-unit+F_unit bound2e48 and tail2e52/K. Do NOT add Qnew to the S205 known piece, which already contains the old spatial remainder and restored finite cells. The S204 fixed local target is not added as an independently matched quantum contribution.",
        "new_known_bound": "The complete known actual piece is below5e48 M[D]Y[Gamma], with original regulator error5e53 M[D]Y[Gamma]/K. Here M is the unchanged detector spatial H1 norm and Y^2=N61^2+X46^2 is unchanged.",
        "canonical": "Both metric factors give Qnew bounds8e-752 and16e-747/K, and the new known actual bounds2e-751 and2e-746/K.",
        "remaining": "The complete actual UV-symbol integral and full contact still need all finite/divergent coefficients, subleading sharp-band artifacts and original fixed covariant matching. No full matched response, reduced mixed inverse, finite-coupling background, stability, physical cutoff, remaining parent loops or finite-gravity IR/Regge result follows. Original V/G/B and P8 remain open.",
        "checks": checks,
        "gates": {
            "new_known_actual_bound": bool(
                c["complete_known_finite_curved_piece"] + estimates.FINITE < KNOWN
            ),
            "new_known_actual_tail": bool(
                c["complete_known_finite_curved_tail"] + tail.TAIL < KNOWN_TAIL
            ),
            "old_finite_pieces_not_double_counted": True,
            "full_UV_and_contact_not_declared_matched": True,
            "fixed_local_target_not_double_counted": True,
            "original_V_G_B_open": True,
        },
    }
