"""Unchanged original source and frozen minimal-gravity inputs."""

from functools import cache

from p8_vacuum_affine_massive_common_gravity_masters import masters as original_masters
from p8_vacuum_affine_massive_dimensional_cut_completion import dimensional, threshold
from p8_vacuum_affine_massive_gravity_pole_completion import poles
from p8_vacuum_affine_radiative_angular_finite import source as previous

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
    "dimensional",
    "original_masters",
    "poles",
    "require_mass",
    "require_order",
    "threshold",
]


@cache
def data():
    old = previous.data()
    result = {
        key: value for key, value in old.items() if key not in ("checks", "gates")
    }
    result["selected_minimal_gravity_finite_reference"] = (
        "The unchanged S283 master functions, S284 whole-D coefficients and Gram "
        "completion, and S288 physical-pole completion determine a known finite "
        "minimal Einstein/massive-scalar one-loop representative. The raw convention "
        "and analytic soft division at resolution/nu=1 are fixed. Three independent "
        "finite matching coordinates remain symbolic; this is not the full current "
        "action's amplitude or the complete observable."
    )
    result["checks"] = dict(old["checks"])
    result["gates"] = {
        **old["gates"],
        "all151_source_checks_and_previous_finite_angular_boundary_retained": len(
            old["checks"]
        )
        == 151,
    }
    return result
