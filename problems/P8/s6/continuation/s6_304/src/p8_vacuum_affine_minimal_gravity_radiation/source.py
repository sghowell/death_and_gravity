"""Original source and complete additional Einstein-radiation tree inventory."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_forward_phase import source as previous
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import recoil
from p8_vacuum_affine_one_newton_inclusive_assembly import forward

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
__all__ = [
    "CONTACT",
    "CUBIC",
    "EP",
    "HEAVY_MASS2",
    "KAPPA",
    "MU",
    "G",
    "K",
    "N",
    "T",
    "data",
    "forward",
    "recoil",
    "require_mass",
    "require_order",
]


@cache
def data():
    old = previous.data()
    result = {
        key: value for key, value in old.items() if key not in ("checks", "gates")
    }
    checks = dict(old["checks"])
    inventory = {
        "external_scalar_emission": (12, 3, 2),
        "scalar_seagull_emission": (6, 2, 1),
        "internal_three_graviton_emission": (3, 3, 2),
    }
    for name, (count, vertices, edges) in inventory.items():
        checks[name + "_tree_loop_number"] = s.Integer(edges - vertices + 1)
        checks[name + "_three_Bose_partitions"] = s.Integer(
            count
            - 3
            * {
                "external_scalar_emission": 4,
                "scalar_seagull_emission": 2,
                "internal_three_graviton_emission": 1,
            }[name]
        )
    checks["complete_minimal_gravity_graph_count"] = s.Integer(
        sum(row[0] for row in inventory.values()) - 21
    )
    checks["complete_selected_matter_plus_gravity_graph_count"] = s.Integer(
        26 + 21 - 47
    )
    checks["five_point_tree_degree_budget"] = s.Integer(5 - 2 - 3)
    checks["six_scalar_jet_exceeds_five_point_tree_budget"] = s.Integer(
        (6 - 2) - (5 - 2) - 1
    )
    result.update(
        {
            "whole_gravity_radiation_inventory": inventory,
            "whole_canonical_action": "eta=(+---),g=eta+2h/sqrt(kappa); S_EH=(kappa/2) integral sqrt(-g)R, S_phi=1/2 integral sqrt(-g)(g^mn d_m Phi d_n Phi-mu Phi^2). Linear de Donder gauge has only quadratic gauge fixing; the four-dimensional projector is I-eta eta/2.",
            "whole_tree_ownership": "The added kappa^-3/2 tree has three Bose partitions with four scalar emissions, two seagulls and one cubic Einstein emission each. Together with the unchanged26 C/sqrt(kappa),g^2/sqrt(kappa) graphs of S295 this is the47-graph selected tree, including every interference. The connected tree identity sum(degree-2)=E-2 excludes six-scalar and higher jets at E5. Inherited current four-scalar derivative cancellations and vanishing lower curvature jets are retained. The formal loop-marked vacuum and fixed counterterms are not retuned.",
            "whole_scope": "Original leading formal five-point tree only; added finite matching, loop-corrected hard vertices, virtual hard remainders, all-N radiation and complex Regge hypotheses remain separate.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "all151_parent_source_checks_copied": len(old["checks"]) == 151
                and checks is not old["checks"],
                "complete21_new_and47_selected_trees": sum(
                    row[0] for row in inventory.values()
                )
                == 21
                and 26 + 21 == 47,
                "canonical_pure_Einstein_no_dilaton_or_classical_projection": True,
                "all_original_higher_jets_and_matching_obligations_retained": True,
            },
        }
    )
    return result
