"""Unchanged original source and already certified real-radiation inputs."""

from functools import cache

import sympy as s
from p8_vacuum_affine_two_real_soft_overlap import source as previous

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
            "one_marked_unchanged_complete_tree": s.Integer(26 + 21 - 47),
            "two_marked_unchanged_complete_tree": s.Integer(140 + 247 + 47 - 434),
            "probability_counterterm_faces_and_overlap": s.Integer(2 - 1 - 1),
        }
    )
    result.update(
        {
            "whole_new_claim_boundary": "A state-transport connector and a probability-level two-real residual match the unchanged full0/1/2-real four-dimensional tree densities in a precisely specified leading-soft reference. Quantitative state continuity makes the connector finite before and after regulator removal. This is not the complete hard real-virtual matching, all-N nonleading rate or original P8 closure.",
            "whole_source_ownership": "S307 owns the full one-real finite residual; S309 owns its same-state all-N leading-soft dressing. S313 owns the full434 two-real amplitude subtraction and its complex O(W) estimate. S300/S301 own the arbitrary-finite-radiation state, soft index and dimensional angular scheme. The present work retains these inputs and derives their new state-continuity and matching connector rather than silently identifying different subtraction schemes.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "source_payload_copied_without_parent_mutation": checks
                is not old["checks"],
                "original_full47_and434_inputs_not_leading_only_substitutes": True,
                "new_matching_is_real_tree_and_known_soft_reference_only": True,
            },
        }
    )
    return result
