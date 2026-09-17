"""Unchanged source and the complete-current factorization boundary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_all_multiplicity_tree_source import source as previous

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
            "whole_new_claim_boundary": "At every finite multiplicity, the full one-off-shell graviton root current built from physical on-shell free leaves is conserved. Each fixed two-real cubic branch factors through this completed current and obeys the conditional S311 angular bound. Matching-forest inclusion-exclusion accounts for overlapping branches. No norm bound on the completed remainder or all-N rate follows.",
            "whole_source_ownership": "S315 owns the unchanged arbitrary finite-multiplicity covariant tree source. S311 owns the exact conserved-pair quotient. This successor adds complete-current factorization, graph-group accounting and explicit counterexamples to isolated-cluster and multi-off-shell Ward shortcuts, without editing either frozen source.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "source_payload_copied_without_parent_mutation": checks
                is not old["checks"],
                "same_selected_action_and_original_parameters": True,
                "factorization_not_all_N_probability_or_original_P8": True,
            },
        }
    )
    return result
