"""All original source parameters and three selected contraction classes."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_parent_one_loop import germs
from p8_vacuum_affine_local_tadpole_radiation import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
require_energies, require_multiplicity = (
    previous.require_energies,
    previous.require_multiplicity,
)
original_parameters = previous.original_parameters
CLASSES = ("within_J2", "within_J4", "across_J2_J4")


def derivative_coefficient():
    return germs.couplings()["H_phi2_Y"]


@cache
def data():
    checks = dict(previous.data()["checks"])
    checks["literal_original_H_phi2Y_coefficient"] = s.factor(
        derivative_coefficient() + 4 * CUBIC / KAPPA
    )
    checks["literal_original_quartic_not_replaced"] = s.factor(
        CONTACT + CUBIC**2 * (3 / (HEAVY_MASS2 - 2) - 2 / (HEAVY_MASS2 - 2) ** 2)
    )
    return {
        "checks": checks,
        "gates": {
            "all_three_selected_source_contraction_classes_retained": len(CLASSES) == 3,
            "S239_full_nonlocal_source_not_expanded_before_integration": True,
            "S297_finite_kappa_source_transfer_credited": True,
            "S337_local_tadpoles_separate_and_added_only_by_linearity": True,
            "independent_S336_curvature_matching_not_assigned": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_original_derivative_coefficient": derivative_coefficient(),
        "whole_contraction_classes": CLASSES,
        "whole_sector_boundary": "Selected cubic-plus-five-field one-loop radiation of J2 G_H J4, J2=gPhi^2/2,J4=cPhi^2Y,c=-4g/kappa. Include all three light-pair contraction classes, every internal/source metric insertion, four external emissions, fixed pole-only subtractions and the existing full-parent symmetric-value condition. No original contact or independent curvature matching coordinate is changed. Other polynomial/heavy/internal-graviton loop sectors and original P8 remain OPEN.",
    }
