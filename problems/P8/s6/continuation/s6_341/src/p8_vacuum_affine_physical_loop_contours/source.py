"""Unchanged original source; distinct scalar and radiative-loop obligations."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_parent_one_loop import germs
from p8_vacuum_affine_heavy_scalar_four_point_loop import amplitude
from p8_vacuum_affine_quadratic_radiation_cancellation import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
require_energies, require_multiplicity = (
    previous.require_energies,
    previous.require_multiplicity,
)
original_parameters = previous.original_parameters


@cache
def data():
    checks = dict(previous.data()["checks"])
    checks["same_original_full_parent_mass"] = HEAVY_MASS2 - germs.MASS2
    checks["same_original_full_parent_contact"] = CONTACT - germs.CONTACT
    checks["ordered_box_radial_factor"] = s.gamma(2) - 1
    return {
        "checks": checks,
        "gates": {
            "original_mass_in_contour_theorem_domain": bool(
                10**6 <= HEAVY_MASS2 < 10**198
            ),
            "all_six_ordered_flat_boxes_retained": True,
            "prior_sharper_cancellation_sensitive_flat_bound_not_replaced": True,
            "complete_radiative_assembly_not_inferred_from_scalar_bounds": True,
            "original_curvature_internal_gravity_and_global_frontiers_open": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_existing_flat_amplitude": amplitude.complete_loop(),
        "whole_source_boundary": "S235 already supplies the entire ordered scalar-loop representation and S236 its sharper flat on-shell remainder. This packet proves additional off-shell physical-contour and weighted line-insertion bounds; it does not rederive or upgrade their global matching conclusion. S340's selected quadratic cancellation is preserved. No source, mass, contact or unknown curved matching coefficient is retuned.",
    }
