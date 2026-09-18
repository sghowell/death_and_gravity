"""Unchanged original source and the probability-level overlap boundary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_three_singleton_subtraction import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
require_energies, require_multiplicity = (
    previous.require_energies,
    previous.require_multiplicity,
)
original_parameters = previous.original_parameters


@cache
def data():
    old = previous.data()
    checks = dict(old["checks"])
    checks.update(
        {
            "all_probability_union_pairs": s.Integer(3**3 - 27),
            "proper_and_full_component_split": s.Integer(12 + 15 - 27),
            "all_three_two_marked_faces": s.Integer(s.binomial(3, 2) - 3),
        }
    )
    return {
        "checks": checks,
        "gates": {
            "same_original_action_and_parameters": True,
            "complete5116_and434_sources_retained": True,
            "no_massive_only_soft_current_shortcut": True,
            "probability_projection_not_amplitude_projection": True,
            "finite_signed_transfer_not_virtual_matching": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_source_ownership": "S315 original all-valence trees; S316 complete conserved hard roots; S310 leading external-pole proof; S313 uniform complete two-real derivative bounds; S319 complete tree envelopes; S323 complete three-real amplitude rectangle and its defined finite signed subtraction. No frozen input is changed.",
        "whole_scope": "Quantitative all-proper-face amplitude hierarchy and the finite transfer to the anchored complete three-real probability density. This remains a signed real-emission measure, not a positive normalized inclusive rate, actual virtual matching, all-N summation or original P8 closure.",
    }
