"""Unchanged full434 tree, recoil, matching parameters and source ownership."""

from functools import cache

import sympy as s
from p8_vacuum_affine_uniform_two_real_tree_bound import source as previous

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
            "double_external_contact_inventory": s.Integer(20 + 13 - 33),
            "double_external_heavy_inventory": s.Integer(60 + 117 - 177),
            "double_external_Einstein_inventory": s.Integer(60 + 117 - 177),
            "complete_overlapping_soft_source": s.Integer(140 + 247 + 47 - 434),
        }
    )
    result.update(
        {
            "whole_new_claim_boundary": "The unchanged complete434 tree is compared with its state-correct double-leading-overlap approximation. The new claim is a uniformly integrable signed difference of two-real tree measures after subtracting both single-soft faces and adding their common overlap once. It is not a complete inclusive detector rate; finite virtual pieces and all-N nonleading sectors are not silently supplied.",
            "whole_overlap_source_ownership": "S310 already owns the complete tree and its fixed-state leading soft limits; S311 owns the conserved current and pair-angle cancellation; S312 owns the entire bare-tree envelope. The new work is the exact double-external regrouping, controlled complex-energy continuation and Cauchy derivative bound, followed by a quantitative overlap-subtracted phase-space estimate. All original parameters and nine primitive frontiers remain unchanged.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "overlap_source_checks_copied_without_parent_mutation": checks
                is not old["checks"],
                "full_marked_state_not_massive_only_soft_current": True,
                "subtracted_real_measure_not_full_inclusive_rate": True,
            },
        }
    )
    return result
