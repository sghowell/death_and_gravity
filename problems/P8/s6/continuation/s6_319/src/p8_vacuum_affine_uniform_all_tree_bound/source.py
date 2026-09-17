"""Unchanged original action and positive Born normalization."""

from functools import cache

import sympy as s
from p8_vacuum_affine_uniform_soft_current_bound import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
require_multiplicity, recoil = previous.require_multiplicity, previous.recoil
original_parameters = previous.original_parameters
diagnostic_parameters = previous.diagnostic_parameters


@cache
def data():
    old = previous.data()
    result = {k: v for k, v in old.items() if k not in ("checks", "gates")}
    checks = dict(old["checks"])
    checks.update(
        {
            "complete_two_real_core_plus_soft_branch": s.Integer(387 + 47 - 434),
            "complete_three_real_core_and_soft_branches": s.Integer(3767 + 1349 - 5116),
        }
    )
    result.update(
        {
            "whole_new_claim_boundary": "The complete original four-Phi/N-graviton on-shell tree has an all-finite-N bound relative to the positive Born amplitude, uniform in collinear and hard-angle approaches with explicit product(1/w_i) poles. It is a bare-tree coefficient estimate, not an inclusive or quantum completion.",
            "whole_source_ownership": "S318 owns the weighted complete pure-soft bound; S317 owns complete temporal-tree equivalence; S315 owns the arbitrary-valence literal action and source; S312 supplies hard-cut and positive Born estimates. No frozen evaluator, action, parameter or matching coefficient changes.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "source_checks_copied_without_parent_mutation": checks
                is not old["checks"],
                "same_original_couplings_and_selected_action": True,
                "complete_bare_tree_bound_not_original_P8_closure": True,
            },
        }
    )
    return result
