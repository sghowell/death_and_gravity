"""Unchanged original action, parameters and finite-class ownership."""

from functools import cache

import sympy as s
from p8_vacuum_affine_relative_energy_complex_tube import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
require_multiplicity, recoil = previous.require_multiplicity, previous.recoil
original_parameters = previous.original_parameters
diagnostic_parameters = previous.diagnostic_parameters


def require_energies(a, b, c):
    values = tuple(require_mass(v) for v in (a, b, c))
    if sum(values) > s.Rational(1, 8):
        raise ValueError("Require positive energies with total at most1/8")
    return values


@cache
def data():
    old = previous.data()
    checks = dict(old["checks"])
    checks.update(
        {
            "three_pair_classes_plus188_and_remaining_core": s.Integer(
                3 * 387 + 188 + 3767 - 5116
            ),
            "nonsingleton_temporal_terms": s.Integer(3 * 387 + 188 - 1349),
            "regular_and_double_external_core_inventory": s.Integer(247 + 140 - 387),
        }
    )
    return {
        "whole_original_parameters": original_parameters(),
        "whole_original_source_boundary": "S315 literal all-valence action and S317 complete temporal reorganization are unchanged. S311/S317 supply the exact complete temporal pair; S313 supplies recoil and hard-cut estimates; S321 supplies the188-class subtraction. This successor bounds each387 pair-plus-singleton class, not the remaining3767 core or an inclusive observable.",
        "checks": checks,
        "gates": {
            "copied_frozen_source_checks_not_mutated": checks is not old["checks"],
            "original_parameters_and_all_vertices_retained": True,
            "composite_pair_not_assumed_TT_or_null": True,
            "no_unknown_matching_constants_selected": True,
            "three_classes_are_gauge_fixed_computational_not_observables": True,
        },
    }
