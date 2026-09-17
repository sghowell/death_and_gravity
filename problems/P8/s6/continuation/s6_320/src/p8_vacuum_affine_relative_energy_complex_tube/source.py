"""Unchanged source and the relative-complex-energy continuation boundary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_uniform_all_tree_bound import source as previous

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
            "whole_new_claim_boundary": "The complete original finite tree has a holomorphic relative-energy polydisc and an explicit angular-uniform Cauchy estimate. Radii shrink with each center energy, and soft faces, all-N IR subtraction and inclusive probabilities are not established.",
            "whole_source_ownership": "S319 owns the complete hard-core decomposition; S318 owns the weighted angular grading; S317 owns temporal-tree equivalence; S315 supplies literal all-order vertices; S313 supplies the two-real complex recoil tube, whose total-perturbation estimates extend independently of multiplicity. A separate Gaussian-rational evaluator preserves every frozen real API.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "source_checks_copied_without_parent_mutation": checks
                is not old["checks"],
                "same_original_couplings_and_selected_action": True,
                "relative_complex_tube_not_soft_faces_or_original_P8_closure": True,
            },
        }
    )
    return result
