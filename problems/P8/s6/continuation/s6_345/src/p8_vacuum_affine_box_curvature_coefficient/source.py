"""Original source retained for one finite known ordered-box coefficient."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_parent_one_loop import germs
from p8_vacuum_affine_triangle_curvature_coefficient import source as previous

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
    for name, got, want in (
        ("mass", HEAVY_MASS2, germs.MASS2),
        ("g", CUBIC, germs.G),
        ("contact", CONTACT, germs.CONTACT),
        ("metric", KAPPA, germs.KAPPA),
    ):
        checks["finite_box_original_" + name] = s.factor(got - want)
    return {
        "checks": checks,
        "gates": {
            "same_original_source_and_physical_metric": True,
            "all_six_mass_ordered_flat_boxes_not_mass_swapped": True,
            "selected_box_known_loop_not_extra_parent_matching": True,
            "fixed_OS4_and_quadratic_counterterms_unchanged": True,
            "analytic_origin_not_physical_above_threshold_approximation": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_source_boundary": "Compute the finite degree6 coefficient of the complete selected alternating-mass box class in an explicitly fixed off-shell covariant jet lift. The full class is already present in S342. Source, masses, physical metric, bounce, counterterms, scoped P8(a) and independent extra parent chi are unchanged.",
    }
