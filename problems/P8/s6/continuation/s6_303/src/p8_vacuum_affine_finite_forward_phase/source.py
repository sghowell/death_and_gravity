"""Unchanged original source, full known finite reference and positive Born input."""

from functools import cache

from p8_vacuum_affine_minimal_gravity_finite import assembly, masters
from p8_vacuum_affine_minimal_gravity_finite import bounds as compact
from p8_vacuum_affine_minimal_gravity_finite import source as previous
from p8_vacuum_affine_one_newton_inclusive_assembly import forward

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
dimensional, poles = previous.dimensional, previous.poles
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
    "assembly",
    "compact",
    "data",
    "dimensional",
    "forward",
    "masters",
    "poles",
    "require_mass",
    "require_order",
]


@cache
def data():
    old = previous.data()
    result = {
        key: value for key, value in old.items() if key not in ("checks", "gates")
    }
    result["selected_all_angle_known_hard_interference"] = (
        "The same S302 known minimal-gravity finite reference is retained, with "
        "all three physical finite anchors still unassigned. The leading forward "
        "Coulomb phase is separated only to bound the real one-loop interference "
        "against the full positive S297 Born amplitude. This does not remove "
        "the complex transfer logarithm, supply gravity-Born radiation or close Regge."
    )
    result["checks"] = dict(old["checks"])
    result["gates"] = {
        **old["gates"],
        "all151_source_checks_and_three_unassigned_anchors_retained": len(old["checks"])
        == 151,
    }
    return result
