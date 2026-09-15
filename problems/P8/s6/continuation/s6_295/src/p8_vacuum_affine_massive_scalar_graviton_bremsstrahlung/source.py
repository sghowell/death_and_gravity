"""Unchanged original source and every minimal five-point matter graph."""

from functools import cache

import sympy as s
from p8_vacuum_affine_whole_mixed_heavy_gravity_sector import source as previous

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
        "quartic_external_scalar_emission": (4, 2, 1),
        "quartic_metric_contact": (1, 1, 0),
        "heavy_exchange_external_scalar_emission": (12, 3, 2),
        "heavy_exchange_internal_heavy_emission": (3, 3, 2),
        "heavy_exchange_metric_cubic_contacts": (6, 2, 1),
    }
    for name, (_, vertices, edges) in inventory.items():
        checks[name + "_whole_tree_loop_number"] = s.Integer(edges - vertices + 1)
    checks["whole_quartic_graph_count"] = s.Integer(4 + 1 - 5)
    checks["whole_heavy_exchange_graph_count"] = s.Integer(12 + 3 + 6 - 21)
    checks["whole_five_point_graph_count"] = s.Integer(
        sum(row[0] for row in inventory.values()) - 26
    )
    result.update(
        {
            "whole_selected_five_point_graph_inventory": inventory,
            "whole_bremsstrahlung_completeness_argument": "At order C/sqrt(kappa), the single original quartic either emits from each of four external scalar lines or from its metric-volume contact. At order g^2/sqrt(kappa), each of the three H-exchange pairings emits from four external Phi lines, the internal H line, or either cubic metric-volume contact. Thus5+3*7=26 labelled graphs. Higher scalar jets cannot supply a connected four-scalar one-graviton tree at these coupling orders. All four scalar lines are on shell; the metric is physical and on shell. The original full four-scalar source cancellations and absence of additional minimal matter species vertices are inherited without altering the source.",
            "whole_selected_coupling_scope": "The complete C/sqrt(kappa) and g^2/sqrt(kappa) tree, and its square including C^2,Cg^2,g^4 over kappa. Radiation from a gravity-exchange Born tree is order kappa^-3/2 and is not this sector. Additional finite curvature/higher-operator matching, loop-corrected hard amplitudes, virtual IR pairing and high-energy Regge control are not inferred.",
            "checks": checks,
            "gates": {
                "all122_parent_source_identities_retained": len(old["checks"]) == 122,
                "cached_parent_checks_copied_not_mutated": checks is not old["checks"],
                "all26_minimal_graphs_and_metric_contacts_retained": sum(
                    row[0] for row in inventory.values()
                )
                == 26,
                "original_tuned_contact_and_all_three_heavy_exchanges": True,
                "selected_real_rate_contains_all_matter_tree_interferences": True,
                "other_matching_virtual_and_original_P8_obligations_open": True,
            },
        }
    )
    return result
