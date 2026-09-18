"""Frozen original current, fixed reference and unchanged physical boundary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_borel_soft_conversion import current
from p8_vacuum_affine_borel_soft_conversion import source as previous
from p8_vacuum_affine_radiative_soft_state_transfer import continuity

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
    t = s.Rational(1, 128)
    return {
        "checks": {
            **dict(previous.data()["checks"]),
            "actual_positive_index_increment": current.index_increment_upper(t)
            * KAPPA
            / t
            - 1700,
            "actual_Born_index_modulus": continuity.index_change_upper(t) * KAPPA / t
            - 1400,
            "actual_positive_conversion_increment": s.simplify(
                current.conversion_increment_upper(t) * KAPPA / (t * (1 - s.log(t)))
                - 11000
            ),
            "actual_Born_conversion_modulus": s.simplify(
                continuity.finite_conversion_change_upper(t)
                * KAPPA
                / (t * (1 - s.log(t)))
                - 10000
            ),
        },
        "gates": {
            "same_original_action_parameters_and_fixed_hard_E_u": True,
            "same_Born_D4_reference_not_state_dependent_intensity": True,
            "only_evaluate_each_state_on_its_own_kinematic_cut": True,
            "limiting_additional_soft_mark_not_outer_hard_matching": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_source_boundary": "Use S331 positive Borel increments and S325 Born moduli for the SAME original fixed-ball additional-soft convention, fixed E,u and physical D4 leading reference. Marks are evaluated only when their own total energy is <=x<=1/8. No dimensionally continued outer hard amplitude or radiation dynamics is supplied.",
    }
