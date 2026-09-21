"""Original source retained; unknown curvature matching is not a new chosen source."""

from functools import cache

from p8_vacuum_affine_dimensional_real_remainder import source as previous

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
    return {
        "checks": dict(previous.data()["checks"]),
        "gates": {
            "original_action_couplings_and_matching_unchanged": True,
            "S239_flat_not_complete_curved_counterfunctional_credited": True,
            "S297_finite_kappa_flat_graph_equality_credited": True,
            "S334_S335_complete_tree_not_one_loop_hard_matching": True,
            "comparison_coordinate_neither_chosen_nor_given_a_zero_default": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_source_boundary": "An explicit formal one-loop curvature matching direction compares effective actions without assigning its physical coefficient or changing the original classical source. It is invisible to retained flat amplitudes through one loop but visible to four-scalar one-graviton radiation. The comparison is not a claim of a UV-complete parent for arbitrary coefficient.",
    }
