"""Original graph ownership and the unchanged massive source."""

from functools import cache

from p8_vacuum_affine_calorimetric_soft_resummation import gamma, poisson
from p8_vacuum_affine_minimal_gravity_radiation import bounds as radiation_bounds
from p8_vacuum_affine_minimal_gravity_radiation import vertices
from p8_vacuum_affine_mixed_gravity_physical_rate import source as previous
from p8_vacuum_affine_physical_virtual_soft_pairing import continuity

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
born = previous.born
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
    "born",
    "continuity",
    "data",
    "gamma",
    "poisson",
    "radiation_bounds",
    "require_mass",
    "require_order",
    "vertices",
]


@cache
def data():
    old = previous.data()
    result = {
        key: value for key, value in old.items() if key not in ("checks", "gates")
    }
    checks = dict(old["checks"])
    checks.update(
        {
            "same_complete_radiation_Newton_coefficient": KAPPA
            - radiation_bounds.source.KAPPA,
            "same_calorimetric_Newton_coefficient": KAPPA - gamma.source.KAPPA,
            "same_physical_conversion_Newton_coefficient": KAPPA
            - continuity.source.KAPPA,
            "same_full_Born_Newton_coefficient": KAPPA - born.source.KAPPA,
        }
    )
    result.update(
        {
            "whole_selected_graph_ownership": "The same21 minimal Einstein and26 massive-matter single-graviton trees are owned by S304/S295. No new diagram enumeration is claimed. This successor regroups the exact external currents and improves the all-angle detector-threshold dependence of their complete real-minus-leading-soft remainder.",
            "whole_hard_and_soft_ownership": "The positive full Born is S297. The all-N leading-soft total-energy sum is S299; its elastic finite physical conversion is S296. S300/S301 control leading-soft kernels of radiative states. Hard finite sectors remain S297/S302/S303/S305/S306 and are not multiplied into a claim about all-N nonleading amplitudes.",
            "whole_reference_convention": "Original mu=nu=1,5/4<=E<=2,all nonforward hard angles,0<x<=1/8. The analytic soft reference is nu=1; the detector threshold is x. The leading IR virtual-real pairing and regulator removal are performed at fixed positive threshold before any x->0 comparison.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "new_source_checks_are_not_parent_alias": checks is not old["checks"],
                "same47_complete_single_real_trees": True,
                "elastic_leading_soft_reference_not_full_detector_rate": True,
                "no_hard_matching_or_multi_real_nonleading_term_chosen": True,
            },
        }
    )
    return result
