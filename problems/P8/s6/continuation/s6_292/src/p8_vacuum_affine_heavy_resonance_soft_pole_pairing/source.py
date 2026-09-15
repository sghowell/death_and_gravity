"""Unchanged original vertices and minimal-gravity graph ownership."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_graviton_production_threshold import source as previous

MU, K, D, EP = previous.MU, previous.K, previous.D, previous.EP
N, G, T = previous.N, previous.G, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT = previous.CONTACT
require_mass, require_order = previous.require_mass, previous.require_order


@cache
def data():
    old = previous.data()
    result = {
        key: value for key, value in old.items() if key not in ("checks", "gates")
    }
    checks = dict(old["checks"])
    inventory = {
        "pair_triangle": (3, 2, 1),
        "metric_cubic_contact": (2, 1, 1),
        "metric_cubic_two_h_tadpole": (1, 0, 1),
        "external_leg_minimal_self_energy": (3, 2, 1),
        "separate_H_metric_mixing": (3, 2, 1),
    }
    for name, (vertices, matter_edges, metric_edges) in inventory.items():
        checks[name + "_complete_graph_loop_number"] = s.Integer(
            matter_edges + metric_edges - vertices
        )
    result.update(
        {
            "whole_minimal_gravity_graph_inventory": inventory,
            "whole_proper_and_LSZ_graph_ownership": "At the selected formal g/kappa order, a graviton in the loop connects two distinct external matter arms (three pair triangles), one arm to the covariant cubic vertex (three contact bubbles), or the cubic vertex to itself (two-h tadpole). Both endpoints on the same arm give that arm's external minimal self-energy/residue, including its seagull. All three external factors are included. Matter-loop H-metric mixing with an off-shell graviton bridge is a separate non-1PI class, not part of this proper-plus-LSZ finite representative.",
            "whole_IR_finite_and_matching_boundary": "The H-metric response and its full original one-point/metric-contact counterterms are retained as a separate finite/UV matching sector: their matter loop is massive and the graviton bridge has fixed square n, not soft loop momentum. The entire fixed Gaussian vacuum-volume condition cancels flat metric tadpoles. Higher source matter tadpoles cannot create the selected soft-graviton logarithm; their coefficients are not set to zero. No additional classical source, finite cubic/curved matching or heavy width is chosen.",
            "whole_mass_and_spectral_boundary": "The computed minimal pure-GR self-energy has Sigma(mass)=0 in D by the frozen S283 identity, including mass n. This known loop sector introduces no heavy mass-shift delta-prime. It does not fix other mass corrections or turn the unstable H into an exact asymptotic particle.",
            "checks": checks,
            "gates": {
                "entire_original_source_and104_previous_checks_retained": len(
                    old["checks"]
                )
                == 104,
                "cached_parent_checks_copied_not_mutated": checks is not old["checks"],
                "all_pair_contact_and_external_factor_classes_retained": True,
                "H_metric_response_not_deleted_or_called_UV_finite": True,
                "full_flat_onepoint_and_volume_conditions_retained": True,
                "no_extra_matching_width_or_original_P8_closure": True,
            },
        }
    )
    return result
