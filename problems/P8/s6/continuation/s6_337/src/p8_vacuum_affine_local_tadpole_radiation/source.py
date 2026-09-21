"""Same original couplings and local six-field sector, not full loop matching."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_parent_one_loop import germs
from p8_vacuum_affine_radiative_curvature_matching import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
require_energies, require_multiplicity = (
    previous.require_energies,
    previous.require_multiplicity,
)
original_parameters = previous.original_parameters
LAMBDA = germs.LAMBDA
NAMES = ("phi6", "phi4Y", "phi2Y2", "Y3", "phi2_L3_minus_L4", "Y_L3_minus_L4")


def local_couplings():
    return {name: value for name, value in germs.couplings().items() if name in NAMES}


@cache
def data():
    checks = dict(previous.data()["checks"])
    for name, value in (
        ("KAPPA", germs.KAPPA),
        ("HEAVY_MASS2", germs.MASS2),
        ("CUBIC", germs.G),
        ("CONTACT", germs.CONTACT),
    ):
        checks["same_full_parent_local_source_" + name] = s.factor(
            globals()[name] - value
        )
    checks["same_original_lambda"] = LAMBDA - s.Rational(1, 10**600)
    return {
        "checks": checks,
        "gates": {
            "all_six_original_local_degree_six_vertices_retained": len(
                local_couplings()
            )
            == 6,
            "S297_finite_kappa_graph_transfer_credited": True,
            "mixed_H_Phi2Y_and_other_loop_sectors_not_absorbed": True,
            "scalar_times_R_source_vertex_zero_only_for_linear_physical_TT": True,
            "independent_S336_curvature_matching_not_selected": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_six_local_couplings": local_couplings(),
        "whole_sector_boundary": "Compute only the known complete six local degree-six tadpole contribution to four-scalar one-graviton radiation in the retained S239 prescription. S297 transfers these vertices to the actual finite-kappa source. The scalar-times-R degree-six source has zero linear on-shell TT response, not zero curved response. Mixed H/source, polynomial/heavy, internal-graviton loops, independent physical curvature matching and original P8 remain open.",
    }
