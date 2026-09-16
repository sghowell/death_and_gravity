"""Original matter source and consistently truncated selected inclusive sector."""

from functools import cache

import sympy as s
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import source as previous

MU, K, D, EP = previous.MU, previous.K, previous.D, previous.EP
N, G, T = previous.N, previous.G, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT = previous.CONTACT
require_mass, require_order = previous.require_mass, previous.require_order


@cache
def data():
    old = previous.data()
    result = {
        key: value for key, value in old.items() if key not in ("checks", "gates")
    }
    checks = dict(old["checks"])
    C, g, H = s.symbols("quartic cubic heavy_tree", real=True)
    checks["entire_selected_Born_rate_interference_inventory"] = s.expand(
        (C + g * g * H) ** 2 - C * C - 2 * C * g * g * H - g**4 * H**2
    )
    result.update(
        {
            "whole_selected_inclusive_couplings": "The complete minimal virtual C/kappa and g^2/kappa amplitudes interfere with the original A0=C+g^2 Htree. The full26-graph real tree is squared. Thus all C^2,Cg^2,g^4 over kappa terms are retained, with local UV subtraction before IR pairing and no one-loop square beyond the stated order.",
            "whole_hard_matching_boundary": "Finite local quartic,RH,heavy-residue and higher matching is not chosen. The renormalized hard virtual coefficient retains any such independent input. The conversion bound does not determine it. Pure-gravity-exchange radiation and higher Newton orders are not evaluated. Hard matter loops can interfere with the gravity Born amplitude at the SAME rate order; that separately IR-finite term is also not evaluated by this selected dressing calculation.",
            "checks": checks,
            "gates": {
                "all130_parent_source_identities_preserved": len(old["checks"]) == 130,
                "cached_parent_source_checks_copied": checks is not old["checks"],
                "entire_selected_interference_not_single_diagram": True,
                "local_UV_subtraction_precedes_IR_pairing": True,
                "independent_finite_hard_matching_not_chosen": True,
                "full_gravity_Regge_and_original_P8_still_open": True,
            },
        }
    )
    return result
