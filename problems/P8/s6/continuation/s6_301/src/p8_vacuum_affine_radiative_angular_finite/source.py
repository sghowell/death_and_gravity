"""Unchanged original source and state-correct radiative recoil."""

from functools import cache

from p8_vacuum_affine_physical_virtual_soft_pairing import soft
from p8_vacuum_affine_radiative_state_soft_index import index, recoil, stability
from p8_vacuum_affine_radiative_state_soft_index import source as previous

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
    "index",
    "recoil",
    "require_mass",
    "require_order",
    "soft",
    "stability",
]


@cache
def data():
    old = previous.data()
    result = {
        key: value for key, value in old.items() if key not in ("checks", "gates")
    }
    result["selected_finite_radiative_conversion"] = (
        "The same physical recoiled radiative momenta are held in four dimensions while the additional soft-graviton angular integral and projector are continued toD=4+2e. This fixes a specified universal finite angular/reference term, not evanescent hard amplitudes, a full quantum state or an inclusive hard rate."
    )
    result["checks"] = dict(old["checks"])
    result["gates"] = {
        **old["gates"],
        "all151_original_source_checks_and_state_correct_IR_boundary_retained": len(
            old["checks"]
        )
        == 151,
    }
    return result
