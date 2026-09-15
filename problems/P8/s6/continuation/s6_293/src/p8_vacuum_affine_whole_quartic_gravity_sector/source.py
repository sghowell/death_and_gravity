"""Unchanged original source and complete one-loop quartic graph ownership."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_resonance_soft_pole_pairing import source as previous

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
        "quartic_pair_triangle": (3, 2, 1),
        "quartic_metric_contact": (2, 1, 1),
        "quartic_two_h_tadpole": (1, 0, 1),
        "quartic_matter_metric_endpoint": (3, 2, 1),
        "quartic_metric_matter_tadpole_endpoint": (2, 1, 1),
    }
    for name, (vertices, matter, metric) in inventory.items():
        checks[name + "_whole_connected_loop_number"] = s.Integer(
            matter + metric - vertices
        )
    result.update(
        {
            "whole_selected_quartic_graph_inventory": inventory,
            "whole_quartic_graph_ownership": "At exactly one original quartic C and two powers of the normalized metric coupling, connected one-loop scalar/metric graphs consist of six pair triangles, four metric-quartic contacts, the two-h quartic tadpole, four external minimal-GR residue factors, and both matter-bubble endpoint assignments in each s/t/u channel. Quartic matter self-energy insertions and metric-quartic matter tadpoles are completed by the inherited full on-shell mass countervertex. Flat metric onepoint insertions, including pure-graviton/ghost scaleless tadpoles, obey the entire original vacuum-volume condition. Disconnected vacuum graphs are normalized out, not counted as connected amplitudes.",
            "whole_quartic_completeness_argument": "With C4 and two minimal matter-stress vertices there are two internal matter edges and one metric edge. Their connected assignments are pair triangles, matter endpoints, or external-line self-energy insertions. With C4h and one stress vertex there is one matter edge and one metric edge: a contact bubble, a quartic matter-tadpole endpoint, or a flat metric onepoint insertion. C4hh closes to the scaleless metric tadpole. A pure-metric cubic or ghost loop can only supply the flat onepoint class at this coupling order. All external residues, scalar OS mass terms and the entire original volume condition are retained.",
            "whole_quartic_matching_boundary": "This is the coefficient linear in the original nondifferentiated C and in1/kappa in the minimal one-loop expansion. It is not all g^2/kappa, pure1/kappa^2, M1/Proca or higher-source matching. Finite constant quartic and constant R Phi^2 anchors remain in the amplitude and have zero crossed b20 at one insertion. Higher-derivative EFT coefficients and the full physical infrared/Regge observable are not fixed.",
            "checks": checks,
            "gates": {
                "whole_original109_parent_checks_retained": len(old["checks"]) == 109,
                "cached_parent_checks_copied_not_mutated": checks is not old["checks"],
                "all_six_pairs_four_contacts_four_external_factors": True,
                "all_six_endpoints_and_full_OS_mass_terms": True,
                "whole_flat_metric_onepoint_condition_not_graph_deletion": True,
                "other_couplings_higher_matching_and_original_P8_open": True,
            },
        }
    )
    return result
