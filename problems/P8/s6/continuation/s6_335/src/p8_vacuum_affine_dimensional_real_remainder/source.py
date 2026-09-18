"""Unchanged source and precise extension of the earlier matter-only limit."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_gravity_radiation import source as previous
from p8_vacuum_affine_one_newton_inclusive_assembly import source as born_source
from p8_vacuum_affine_physical_virtual_soft_pairing import soft as phase

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
        checks["same_positive_Born_source_" + name] = s.simplify(
            globals()[name] - getattr(born_source, name)
        )
    checks["same_original_phase_at_zero"] = phase.phase_normalization(0) - 1
    return {
        "checks": checks,
        "gates": {
            "same_original_action_and_couplings_not_auxiliary_species": True,
            "S296_matter_only_dimensional_integral_credited": True,
            "S304_full_physical_recoil_gap_and_D4_tree_credited": True,
            "S334_full_D_tree_and_canonical_trace_response_used": True,
            "same_D_complete_Born_times_same_D_soft_current_subtracted": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_source_boundary": "Extend the selected47 original tree real-minus-soft integral, not the separate divergent rates or virtual hard matching. Use the same positive D=4+2epsilon convention, original masses/couplings and matching coordinates, exact recoil and full polarization sew. S296 already supplies the matter-only limit and fixed-domain measure; S304 supplies the physical gap; S334 supplies the missing gravity trace response.",
    }
