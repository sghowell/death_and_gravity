"""Original fixed-ball soft convention and unchanged hard-source frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_borel_soft_conversion import current
from p8_vacuum_affine_physical_virtual_soft_pairing import soft
from p8_vacuum_affine_radiative_angular_finite import angular
from p8_vacuum_affine_radiative_soft_state_transfer import continuity
from p8_vacuum_affine_remaining_log_cutoff import source as previous

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
    t, e = s.Rational(1, 128), s.Symbol("epsilon", positive=True)
    return {
        "checks": {
            **dict(previous.data()["checks"]),
            "actual_positive_regulated_quotient_modulus": s.simplify(
                current.regulated_increment_upper(t) * KAPPA / (t * (1 - s.log(t)))
                - 18000
            ),
            "actual_Born_regulated_quotient_modulus": s.simplify(
                continuity.regulator_quotient_upper(t) * KAPPA / (t * (1 - s.log(t)))
                - 15000
            ),
            "actual_phase_times_radial_beta": s.simplify(
                soft.phase_normalization(e) * angular.beta_multiplier(e)
                - (4 * s.pi) ** (-e) / s.gamma(1 + e)
            ),
        },
        "gates": {
            "same_original_action_and_fixed_hard_E_u": True,
            "same_Born_D4_reference_and_own_cut_for_each_mark": True,
            "full_additional_soft_projector_radial_and_phase_retained": True,
            "joint_soft_mark_bound_not_outer_hard_matching": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_source_boundary": "Use the exact S296/S301 phase and fixed-ball additional-soft convention, S325 Born moduli, S331 positive Borel increments and S332 common-space conditioning split. Both limiting procedures act on this SAME named soft insertion; no outer hard D-dimensional amplitude or state-dependent intensity is inferred.",
    }
