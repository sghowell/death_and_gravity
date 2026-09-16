"""Original fixed source and inherited same-prescription matter graph inventory."""

from functools import cache

from p8_vacuum_affine_one_newton_inclusive_assembly import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT = previous.CONTACT
require_mass, require_order = previous.require_mass, previous.require_order


@cache
def data():
    old = previous.data()
    result = {
        key: value for key, value in old.items() if key not in ("checks", "gates")
    }
    result["selected_Regge_input"] = (
        "Only the generated g-squared matter stress F1 and its PhiPhi/HH cuts are used for the new spin-two comparison. All original source vertices and the finite-kappa graph inventory remain unchanged; no UV tower is asserted."
    )
    result["checks"] = dict(old["checks"])
    result["gates"] = {
        **old["gates"],
        "all151_parent_source_checks_preserved": len(old["checks"]) == 151,
    }
    return result
