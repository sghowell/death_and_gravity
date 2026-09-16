"""Original parameters and all inherited physical-source obligations."""

from functools import cache

from p8_vacuum_affine_calorimetric_soft_resummation import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order


@cache
def data():
    old = previous.data()
    result = {
        key: value for key, value in old.items() if key not in ("checks", "gates")
    }
    result["selected_radiative_state"] = (
        "The soft logarithmic coefficient is evaluated on the exact recoiled four-massive-scalar plus any finite number of positive-energy outgoing gravitons. A continuous finite angular-energy-measure extension is only a kinematic coefficient construction, not the original quantum state or full physical rate."
    )
    result["checks"] = dict(old["checks"])
    result["gates"] = {
        **old["gates"],
        "all151_original_source_checks_and_prior_soft_boundary_retained": len(
            old["checks"]
        )
        == 151,
    }
    return result
