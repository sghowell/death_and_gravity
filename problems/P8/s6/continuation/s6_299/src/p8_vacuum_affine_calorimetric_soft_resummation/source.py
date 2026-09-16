"""Original source and the already matched full D-dimensional soft kernel."""

from functools import cache

from p8_vacuum_affine_physical_virtual_soft_pairing import soft
from p8_vacuum_affine_unequal_mass_regge_residue import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT = previous.CONTACT
EP = soft.EP
require_mass, require_order = previous.require_mass, previous.require_order


@cache
def data():
    old = previous.data()
    result = {
        key: value for key, value in old.items() if key not in ("checks", "gates")
    }
    result["selected_leading_soft_observable"] = (
        "Four fixed on-shell massive scalar hard momenta; sum over any number of leading soft gravitons whose TOTAL COM energy is at most E. The unchanged analytic virtual convention uses S278/S288 and the complete evanescent angular conversion uses S296. This is not a bound on all non-leading soft or hard terms."
    )
    result["checks"] = dict(old["checks"])
    result["gates"] = {
        **old["gates"],
        "all151_source_checks_and_both_previous_cut_obligations_preserved": len(
            old["checks"]
        )
        == 151,
    }
    return result
