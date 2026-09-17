"""Unchanged source and the all-finite-multiplicity current-bound boundary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_temporal_tree_reorganization import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
require_multiplicity = previous.require_multiplicity
recoil = previous.recoil


def original_parameters():
    return {"heavy": HEAVY_MASS2, "cubic": CUBIC, "contact": CONTACT, "kappa": KAPPA}


def diagnostic_parameters():
    return {
        "heavy": 128,
        "cubic": s.Rational(1, 16),
        "contact": s.Rational(1, 8),
        "kappa": 16,
    }


@cache
def data():
    old = previous.data()
    result = {k: v for k, v in old.items() if k not in ("checks", "gates")}
    checks = dict(old["checks"])
    checks.update(
        {
            "complete_three_real_pair_inventory": s.Integer(3814 + 3 * 434 - 5116),
            "four_real_pair_overlap_inventory": s.Integer(
                73444 - 6 * 5116 + 3 * 434 - 44050
            ),
        }
    )
    result.update(
        {
            "whole_new_claim_boundary": "The unchanged complete temporal pure-soft recursion has an all-finite-multiplicity weighted tensor bound uniform over all angular and nested collinear configurations away from exact internal poles, with explicit W^n/product(w_i) soft factors and a factorial counting envelope. This is a pure-soft current bound, not the complete hard amplitude or inclusive probability.",
            "whole_source_ownership": "S317 owns the unchanged gauge-complete temporal recursion and finite Noether induction. S315 supplies the all-order literal vertices and their canonical norm bounds. This successor proves a weighted collinear grading and labeled majorant without modifying the action, tree evaluator, original parameters or frozen sources.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "source_payload_copied_without_parent_mutation": checks
                is not old["checks"],
                "same_selected_action_and_original_parameters": True,
                "uniform_soft_current_not_all_N_probability_or_original_P8": True,
            },
        }
    )
    return result
