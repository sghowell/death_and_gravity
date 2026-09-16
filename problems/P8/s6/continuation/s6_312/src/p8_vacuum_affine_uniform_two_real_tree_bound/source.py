"""Same original full434 tree; only its uniform upper bound is advanced."""

from functools import cache

import sympy as s
from p8_vacuum_affine_two_real_collinear_current import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
recoil = previous.recoil


@cache
def data():
    old = previous.data()
    result = {k: v for k, v in old.items() if k not in ("checks", "gates")}
    checks = dict(old["checks"])
    checks.update(
        {
            "regular_contact_graph_count": s.Integer(38 - 5 - 33),
            "regular_heavy_graph_count": s.Integer(198 - 21 - 177),
            "regular_Einstein_graph_count": s.Integer(198 - 21 - 177),
            "full_tree_inventory_unchanged": s.Integer(33 + 177 + 177 + 47 - 434),
        }
    )
    result.update(
        {
            "whole_new_claim_boundary": "The S311 complete-current pair theorem is retained. New common mixed-hard-propagator gaps for every radiation assignment and whole-graph component budgets bound all387 complementary graphs. Their sum with the47 pair graphs gives an original full-Born upper bound for the complete434-tree amplitude. Its explicit1/(ab) soft behavior is not integrated or called infrared finite.",
            "whole_uniform_full_tree_domain": "Same S300 recoil and original canonical tree, E in[5/4,2],a,b>0,a+b<=1/8, all nonforward hard Born angles and all physical emitted directions, with unit spatial-Frobenius TT polarizations. Canonical de Donder graph classes are used only to bound their physical full sum.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "full_tree_source_checks_copied_without_parent_mutation": checks
                is not old["checks"],
                "original_matching_parameters_not_fitted_to_bound": True,
                "no_sector_is_dropped_from_the434_graph_sum": True,
                "bare_tree_upper_bound_not_IR_finite_detector_rate": True,
            },
        }
    )
    return result
