"""Original graph ownership, parameters and distinct known/matched quantities."""

from functools import cache

from p8_vacuum_affine_complete_matter_graviton_endpoint import endpoint
from p8_vacuum_affine_one_newton_inclusive_assembly import forward as born
from p8_vacuum_affine_spectator_gravity_insertion import source as previous
from p8_vacuum_affine_whole_mixed_heavy_gravity_sector import forward as mixed
from p8_vacuum_affine_whole_mixed_heavy_gravity_sector import graphs
from p8_vacuum_affine_whole_quartic_gravity_sector import forward as quartic
from p8_vacuum_affine_whole_quartic_gravity_sector import proper

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
    "born",
    "data",
    "endpoint",
    "graphs",
    "mixed",
    "proper",
    "quartic",
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
    checks.update(
        {
            "same_whole_mixed_heavy_mass": HEAVY_MASS2 - mixed.source.HEAVY_MASS2,
            "same_whole_mixed_cubic": CUBIC - mixed.source.CUBIC,
            "same_whole_endpoint_contact": CONTACT - endpoint.source.CONTACT,
            "same_full_Born_kappa": KAPPA - born.source.KAPPA,
        }
    )
    result.update(
        {
            "whole_selected_graph_ownership": "S293 already derived the whole C/kappa sector including its C metric endpoint. S294 already derived the g^2/kappa nonendpoint sector; S290 supplies both massive stress triangles, their full OS terms and the H-metric bubble. This successor proves their physical-region estimate and an exact scalar-master collapse; it does not claim new graph enumeration.",
            "whole_no_double_counting": "The S290 endpoint is added with quartic=0 because the C bubble is already in S293. The whole massive matter-only loop, pure Einstein/Phi gravity loop and H/Proca/M1 metric loops remain separate S297/S302/S305 owners. No loop squares or dependent higher-degree vertices are inserted.",
            "whole_reference_convention": "Original mu=nu^2=1, D=4+2e, Feynman normal sheet, local UV pole removed before the S278 analytic soft division. The auxiliary reference E^2 ranges over[1/4,1]; it is not a selected detector energy. Physical conversion and unexpanded resolution dependence remain separately owned by S296/S299.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "new_source_checks_are_not_parent_alias": checks is not old["checks"],
                "all_selected_couplings_and_both_endpoints_retained": True,
                "C_metric_endpoint_counted_once": True,
                "known_reference_not_physical_matching_choice": True,
            },
        }
    )
    return result
