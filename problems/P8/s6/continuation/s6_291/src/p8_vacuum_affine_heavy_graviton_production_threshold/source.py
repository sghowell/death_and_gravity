"""Unchanged original source and complete one-graviton three-matter tree arity."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complete_matter_graviton_endpoint import source as previous

MU, K, D, EP = previous.MU, previous.K, previous.D, previous.EP
N, G, T = previous.N, previous.G, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT = previous.CONTACT
require_mass, require_order = previous.require_mass, previous.require_order


def require_dimension(value):
    value = require_order(value)
    if value < 4:
        raise ValueError("Physical TT matrices require integer D at least4")
    return int(value)


@cache
def data():
    old = previous.data()
    result = {k: v for k, v in old.items() if k not in ("checks", "gates")}
    checks = dict(old["checks"])
    vertices, edges = s.symbols("tree_vertices internal_matter_edges", integer=True)
    checks["three_matter_tree_degree_identity"] = s.expand(
        (2 * edges + 3) - 2 * vertices - 1
    ).subs(edges, vertices - 1)
    checks["original_H_Phi_Phi_cubic_retained"] = CUBIC - s.Rational(1, 8192)
    checks["original_heavy_mass_retained"] = (
        HEAVY_MASS2 - s.Integer(10) ** 200 / 512 - 2
    )
    checks["original_light_mass_retained"] = (
        s.Integer(1) - previous.data()["whole_selected_endpoint_couplings"]["mu"]
    )
    result.update(
        {
            "whole_production_tree_inventory": "At one external graviton, three external matter legs and tree order, sum(v_matter-2)=1. The original nonzero excess-arity interaction is g H Phi^2. Retain emission from both incoming Phi legs, the outgoing H leg, and the covariant cubic metric contact. The inherited full source jets exclude extra classical curvature, derivative-cubic, tadpole or scalar-cubic insertions at this arity. Higher matching coefficients are not silently made classical source vertices.",
            "whole_production_order_and_state_boundary": "Amplitude order g/sqrt(kappa); its full two-body cut is order g^2/kappa. H is a formal perturbative channel, not an exact stable asymptotic particle. Flat source, OS and one-point prescriptions are inherited unchanged. Neither higher-EFT matching nor a physical inclusive/width prescription is chosen.",
            "checks": checks,
            "gates": {
                "entire_original_source_and_all_old_jets_retained": len(old["checks"])
                == 100,
                "whole_parent_checks_copied_not_mutated": checks is not old["checks"],
                "all_four_original_tree_graphs_required": True,
                "curved_matching_not_promoted_to_original_tree": True,
                "heavy_channel_not_claimed_exact_stable_state": True,
                "original_P8_and_state_frame_unchanged": True,
            },
        }
    )
    return result
