"""Original source retained while one local matching direction is classified."""

from functools import cache

import sympy as s
from p8_vacuum_affine_radiative_curvature_matching import source as original
from p8_vacuum_affine_triangle_box_radiation import source as previous

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
    for name in ("KAPPA", "HEAVY_MASS2", "CUBIC", "CONTACT"):
        checks["curvature_basis_same_original_" + name] = s.factor(
            globals()[name] - getattr(original, name)
        )
    return {
        "checks": checks,
        "gates": {
            "same_massive_source_and_physical_metric": True,
            "same_named_flat_OS4_not_a_curved_matching_condition": True,
            "formal_curvature_coordinate_not_a_new_selected_source": True,
            "selected_matter_poles_not_full_internal_gravity": True,
            "all_original_quantum_and_background_obligations_open": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_source_boundary": "Classify the lowest parity-even explicit-curvature local contact for the unchanged original massive scalar and physical metric. Recompute known selected-matter pole identities; do not change finite counterterms, chi, the source, the bounce, or scoped P8(a).",
    }
