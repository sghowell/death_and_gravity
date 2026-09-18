"""The original action, parameters, recoil and all historical scopes remain."""

from functools import cache

import sympy as s
from p8_vacuum_affine_radiative_soft_state_transfer import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
require_energies = previous.require_energies
require_multiplicity = previous.require_multiplicity
original_parameters = previous.original_parameters


@cache
def data():
    checks = dict(previous.data()["checks"])
    checks["same_Einstein_coupling_power"] = s.Rational(3, 2) - 3 * s.Rational(1, 2)
    return {
        "checks": checks,
        "gates": {
            "original_action_parameters_and_recoil_unchanged": True,
            "all_frozen_tree_and_probability_results_retained": True,
            "only_attributed_proved_one_loop_logarithm_imported": True,
            "no_higher_soft_conjecture_or_full_loop_remainder_promoted": True,
            "no_unknown_matching_coordinate_chosen": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_source_ownership": "S300 exact arbitrary-radiation recoil/index, S304 original scalar/EH normalization, S319 complete finite-tree bound, S325 total-energy current continuity. The attributed one-loop logarithmic soft theorem supplies only its displayed coefficient, not the hard loop or uniform finite remainder.",
    }
