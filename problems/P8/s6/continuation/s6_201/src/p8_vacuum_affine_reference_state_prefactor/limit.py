"""Actual reference change and unchanged remaining curved matching sector."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_affine_retarded_reference_boundary import bounds as finite
from p8_vacuum_affine_spatial_current.bounds import require_band

from . import contact, memory

DISPLAY = s.Integer(2) * 10**12
TAIL = s.Integer(10) ** 17


def tail_bound(K):
    return TAIL / require_band(K)


@cache
def constants():
    mc, cc = memory.constants(), contact.constants()
    return {
        "full_prefactor_memory_contact_bound": mc["full_memory_H1_coefficient"]
        + cc["contact_L2_bound"],
        "full_prefactor_regulator_tail_numerator": mc["low_external_K_numerator"]
        + mc["high_external_K_numerator"]
        + cc["full_band_K_tail_numerator"],
        "actual_to_unit_reference_display": s.Integer(3) * 10**23,
        "actual_to_unit_reference_tail": s.Integer(2) * 10**27,
        "complete_known_finite_curved_piece": s.Integer(2) * 10**48,
        "complete_known_finite_curved_tail": s.Integer(2) * 10**52,
    }


@cache
def data():
    actual, weighted, unit = s.symbols("J_actual J_alpha_W8 J_unit_W8")
    c = constants()
    checks = {
        "exact_full_current_telescoping": s.expand(
            (actual - weighted) + (weighted - unit) - (actual - unit)
        ),
        "canonical_prefactor_complete_bound": 4 * DISPLAY / modes.KAPPA
        - 8 * s.Rational(1, 10) ** 788,
        "canonical_prefactor_regulator_tail": 4 * TAIL / modes.KAPPA
        - 4 * s.Rational(1, 10) ** 783,
        "canonical_actual_to_unit_reference": 4
        * c["actual_to_unit_reference_display"]
        / modes.KAPPA
        - 12 * s.Rational(1, 10) ** 777,
        "canonical_actual_to_unit_reference_tail": 4
        * c["actual_to_unit_reference_tail"]
        / modes.KAPPA
        - 8 * s.Rational(1, 10) ** 773,
        "canonical_complete_known_finite_piece": 4
        * c["complete_known_finite_curved_piece"]
        / modes.KAPPA
        - 8 * s.Rational(1, 10) ** 752,
        "canonical_complete_known_finite_tail": 4
        * c["complete_known_finite_curved_tail"]
        / modes.KAPPA
        - 8 * s.Rational(1, 10) ** 748,
        "prefactor_regulator_display_at_proof_partition": tail_bound(10**16) - 10,
    }
    return {
        "prefactor_result": "The complete reference change R_alpha has an absolute all-momentum limit with |R_alpha|<2e12 M[D]M[Gamma], and full same-regulator error<1e17 M[D]M[Gamma]/K for K>=1000. Both canonical displays are8e-788 and4e-783/K.",
        "actual_unit_comparison": "At each finite regulator J_actual-J_unit=(J_actual-J_alpha)+(J_alpha-J_unit). Combining S197 gives full actual-to-unit-W8 remainder<3e23 M[D]M[Gamma], with error<2e27 M[D]M[Gamma]/K. No actual state term is changed or omitted.",
        "unit_reference_extraction": "Unit-W8 pair amplitudes satisfy the SAME upper bounds Cref=1e9 and time-analytic inverse-phase estimates as S198, since no |alpha0|<2 enlargement is needed. Repeating that exact six-step identity gives F_unit=B5_unit+bulk6_unit with its existing1e48 bound and1e52/K error.",
        "exact_curved_decomposition": "J_actual,K=R_actual-unit,K+C_unit,K+sum_(j=0)^4 B_j,unit,K+F_unit,K. The full contact remains one-mode regulated and each boundary/memory term retains both created momenta.",
        "complete_known_piece": "The already controlled piece R_actual-unit+F_unit is below2e48 M[D]N61[Gamma], with regulator error2e52 M[D]N61[Gamma]/K and canonical displays8e-752 and8e-748/K. The unknown contact-plus-five-endpoint sector is not included.",
        "reference_boundary": "The remaining unit-W8 reference has canonical Wronskian but is not an exact bisolution or a new physical state. Initial Borel momentum prefactors are removed only from this comparison sector, without assuming or differentiating their momentum analyticity.",
        "remaining": "The joint spatial analytic bounds, full reference spatial Taylor/subtraction expansion and original fixed covariant contact/endpoint matching still require proof. Full mixed inverse, interacting background, stability, cutoff and original V/G/B remain OPEN.",
        "constants": c,
        "checks": checks,
        "gates": {
            "complete_prefactor_bound": c["full_prefactor_memory_contact_bound"]
            < DISPLAY,
            "complete_prefactor_tail": c["full_prefactor_regulator_tail_numerator"]
            < TAIL,
            "actual_to_unit_reference_bound": 2 * 10**23 + DISPLAY
            < c["actual_to_unit_reference_display"],
            "actual_to_unit_reference_error": 10**27 + TAIL
            < c["actual_to_unit_reference_tail"],
            "unchanged_full_finite_piece_bound": finite.DISPLAY
            + c["actual_to_unit_reference_display"]
            < c["complete_known_finite_curved_piece"],
            "unchanged_full_finite_piece_tail": finite.TAIL
            + c["actual_to_unit_reference_tail"]
            < c["complete_known_finite_curved_tail"],
        },
    }
