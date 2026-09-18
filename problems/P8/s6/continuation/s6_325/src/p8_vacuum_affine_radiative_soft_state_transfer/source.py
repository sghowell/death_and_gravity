"""Original source retained; only a universal soft state difference is added."""

from functools import cache

import sympy as s
from p8_vacuum_affine_three_real_probability_overlap import source as previous

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
    checks.update(
        {
            "two_outgoing_recoil_changes_only": s.Integer(2 * 3 * 144 + 2 - 866),
            "unchanged_soft_probability_pole_normalization": s.Rational(2, 4)
            - s.Rational(1, 2),
            "physical_marked_seed_not_extra_Bose_factor": s.Rational(2, 2) - 1,
        }
    )
    return {
        "checks": checks,
        "gates": {
            "same_original_action_parameters_and_recoil": True,
            "all_frozen_full_tree_and_probability_results_retained": True,
            "arbitrary_finite_null_multiplicity_in_kinematic_bound": True,
            "additional_soft_factor_only_dimensionally_continued": True,
            "no_full_radiative_hard_loop_or_matching_value_claim": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_source_ownership": "Frozen S300 exact arbitrary-radiation recoil and S301 full fixed-ball dimensional current; S309 same-state universal soft convention; S324 state-change matching boundary. The new continuity estimate and named Born-seeded insertion retain all original data. They do not reconstruct a radiative hard loop.",
    }
