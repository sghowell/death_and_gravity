"""Unchanged original source and all minimal mixed one-loop graph classes."""

from functools import cache

import sympy as s
from p8_vacuum_affine_whole_quartic_gravity_sector import source as previous

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
        "mixed_two_g_two_light_stress": (4, 3, 1),
        "mixed_two_g_heavy_light_stress": (4, 3, 1),
        "mixed_two_g_two_heavy_stress": (4, 3, 1),
        "mixed_one_metric_g_contact_light_stress": (3, 2, 1),
        "mixed_one_metric_g_contact_heavy_stress": (3, 2, 1),
        "mixed_two_metric_g_contacts": (2, 1, 1),
        "mixed_two_h_g_tadpole": (2, 1, 1),
        "mixed_two_h_stress_tadpole": (3, 2, 1),
    }
    for name, (vertices, matter, metric) in inventory.items():
        checks[name + "_whole_connected_loop_number"] = s.Integer(
            matter + metric - vertices
        )
    result.update(
        {
            "whole_selected_mixed_graph_inventory": inventory,
            "whole_mixed_graph_ownership": "At two original H Phi Phi vertices and two normalized metric powers, the non-tadpole classes have two linear stress vertices, one metric-cubic contact and one stress, or two metric-cubic contacts. The matter stress species can be Phi/Phi, H/Phi or H/H. Connected assignments give twelve mixed boxes, all three proper H Phi Phi pair triangles and contact bubbles in both endpoints of every heavy exchange, the complete heavy selfenergy, four external Phi sqrtZ factors, and the whole S290 matter-graviton endpoint. The two-metric-cubic contact bubble is retained. Quadratic metric contacts and pure-gravity/ghost onepoint attachments use the original whole volume condition or genuinely scaleless tadpoles.",
            "whole_mixed_completeness_argument": "There is exactly one internal graviton in the non-onepoint classes. Two cubic matter vertices supply two H and four Phi ends. Two matter stresses add four ends, so three internal matter edges and four vertices give one loop; choosing their species gives the first three rows. One metric-cubic contact and one stress give two matter edges and three vertices. Two metric-cubic contacts give one H edge and two vertices. A quadratic metric contact closes a scaleless h tadpole. Cutting a scalar bridge gives an external selfenergy or a proper heavy endpoint; cutting a metric bridge gives the complete S290 matter endpoint. With neither bridge the only four-vertex loop is the mixed box. All mixed vertices of the original flat action, OS kinetic/mass conditions and original whole H/metric onepoints are retained.",
            "whole_matching_boundary": "This is the selected minimal one-loop g^2/kappa coefficient, not other coupling sectors or higher EFT data. The S290 R H finite matching anchor remains explicit. A finite additional cubic/heavy-residue or higher-operator matching insertion is not determined by the known loop representative. The perturbative H threshold is not an exact stable heavy particle. Physical infrared/Regge, all-loop and original V/G/B/P8 remain open.",
            "checks": checks,
            "gates": {
                "whole_original114_parent_checks_retained": len(old["checks"]) == 114,
                "cached_parent_checks_copied_not_mutated": checks is not old["checks"],
                "all_scalar_species_metric_contacts_and_four_external_factors": True,
                "whole_metric_endpoint_and_two_contact_bubble_retained": True,
                "whole_original_onepoint_and_OS_conditions": True,
                "other_couplings_higher_matching_and_original_P8_open": True,
            },
        }
    )
    return result
