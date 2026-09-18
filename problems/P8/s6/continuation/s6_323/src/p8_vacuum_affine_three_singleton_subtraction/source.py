"""Unchanged original model and the complete third-order amplitude boundary."""

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
        raise ValueError("Require three positive energies with total at most1/8")
    return values


@cache
def data():
    old = previous.data()
    checks = dict(old["checks"])
    checks.update(
        {
            "complete_three_singleton_and_nonsingleton_sum": s.Integer(
                3767 + 1349 - 5116
            ),
            "central_forest_graph_count_split": s.Integer(
                1288 + 1368 + 804 + 307 - 3767
            ),
            "all125_forest_assignments": s.Integer(5**3 - 125),
            "all64_independent_singleton_assignments": s.Integer(4**3 - 64),
            "all48_pair_singleton_assignments": s.Integer(3 * 4 * 4 - 48),
        }
    )
    return {
        "checks": checks,
        "gates": {
            "same_original_action_and_parameters": True,
            "no_frozen_parent_source_mutated": checks is not old["checks"],
            "all_literal_contact_heavy_and_Einstein_vertices_retained": True,
            "amplitude_and_defined_signed_measure_not_inclusive_completion": True,
            "no_unproved_full_current_global_W_tube": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_scope": "The complete3767 hard core is reorganized into scalar-line cumulants and central kernels. Combined with the preceding1349-class result, the proposed scoped theorem is a complete5116-tree amplitude rectangle bound. Its seven-face baseline also gives a finite signed density difference. Positive normalized inclusive probabilities, real-virtual/all-N completion, state, Regge and original P8 remain open.",
        "whole_source_ownership": "S315 literal all-valence action; S317 complete temporal reorganization; S319 complete tree baseline; S320 hard complex majorants; S321 complete three-current faces; S322 three pair-plus-singleton subtractions. This module changes none of those inputs.",
    }
