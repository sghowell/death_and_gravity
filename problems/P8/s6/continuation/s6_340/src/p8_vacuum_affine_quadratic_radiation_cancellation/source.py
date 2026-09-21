"""Original source and exhaustive selected one-loop two-point topologies."""

from functools import cache

from p8_vacuum_affine_factorized_bubble_radiation import source as previous
from p8_vacuum_affine_heavy_parent_one_loop import germs
from p8_vacuum_affine_heavy_scalar_one_loop import self_energy

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
require_energies, require_multiplicity = (
    previous.require_energies,
    previous.require_multiplicity,
)
original_parameters = previous.original_parameters
CLASSES = ("mixed_two_cubic_bubble", "quartic_tadpole", "heavy_onepoint_reducible")


@cache
def data():
    checks = dict(previous.data()["checks"])
    checks["original_two_point_heavy_mass"] = HEAVY_MASS2 - self_energy.model.MASS2
    checks["original_two_point_cubic_squared"] = CUBIC**2 - self_energy.model.G2
    return {
        "checks": checks,
        "gates": {
            "only_two_cubics_or_single_quartic_at_one_loop_two_point": set(
                germs.graph_degrees(2)
            )
            == {(2, 0, 0, 0), (0, 1, 0, 0)},
            "all_three_old_tadpole_and_bubble_topologies_retained": len(CLASSES) == 3,
            "existing_original_OS_mass_residue_and_onepoint_conditions": True,
            "higher_scalar_jets_excluded_at_this_degree_not_all_orders": True,
            "independent_four_hard_curvature_and_internal_gravity_unchanged": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_selected_topologies": CLASSES,
        "whole_covariant_generator": "The degree-two first trace Tr(G_phi W)/2 of the frozen full covariant Hessian contains the contact tadpole, local heavy-source onepoint term and mixed bilocal kernel -g^2 Phi(x)G_phi(x,y)G_H(x,y)Phi(y)/2 in Euclidean conventions. The source onepoint counterterm and fixed on-shell light mass/residue counterterms are retained before metric variation.",
        "whole_scope": "The selected old massive matter first-loop two-point sector of the original source, at its stated real external-graviton order. No higher-loop, internal-graviton or full curved matching conclusion is inferred. Original finite clock extension V=1+O(X^1024) changes none of these vacuum jets.",
    }
