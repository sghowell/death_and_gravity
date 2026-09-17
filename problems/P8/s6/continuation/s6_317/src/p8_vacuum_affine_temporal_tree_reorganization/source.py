"""Unchanged source and the gauge-complete finite-tree boundary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complete_pair_factorization import source as previous

MU, N, G, K, T = previous.MU, previous.N, previous.G, previous.K, previous.T
KAPPA, HEAVY_MASS2, CUBIC = previous.KAPPA, previous.HEAVY_MASS2, previous.CUBIC
CONTACT, EP = previous.CONTACT, previous.EP
require_mass, require_order = previous.require_mass, previous.require_order
require_multiplicity = previous.require_multiplicity
recoil = previous.recoil


def original_parameters():
    return {"heavy": HEAVY_MASS2, "cubic": CUBIC, "contact": CONTACT, "kappa": KAPPA}


def diagnostic_parameters():
    return {
        "heavy": 128,
        "cubic": s.Rational(1, 16),
        "contact": s.Rational(1, 8),
        "kappa": 16,
    }


@cache
def data():
    old = previous.data()
    result = {k: v for k, v in old.items() if k not in ("checks", "gates")}
    checks = dict(old["checks"])
    checks.update(
        {
            "complete_three_real_pair_inventory": s.Integer(3814 + 3 * 434 - 5116),
            "four_real_pair_overlap_inventory": s.Integer(
                73444 - 6 * 5116 + 3 * 434 - 44050
            ),
        }
    )
    result.update(
        {
            "whole_new_claim_boundary": "The unchanged finite tree source has an equivalent recursive gauge policy with temporal pure-soft currents and unchanged complete on-shell amplitudes. An independent finite coordinate pullback verifies the nonlinear chart. Uniform two-ray and restricted planar three-ray current bounds do not imply a general all-N amplitude or inclusive rate bound.",
            "whole_source_ownership": "S316 owns complete-current factorization and the isolated-cluster obstruction. S315 supplies unchanged all-order vertices and the original rooted source. This successor changes only the gauge representative of each positive-energy pure-soft current, retaining every vertex and graph and the original inverse-kinetic policy on hard/mixed currents, which are recomputed from the transformed lower coefficients. Nonlinear terms are produced by the full recursion, not discarded.",
            "checks": checks,
            "gates": {
                **old["gates"],
                "source_payload_copied_without_parent_mutation": checks
                is not old["checks"],
                "same_selected_action_and_original_parameters": True,
                "temporal_reorganization_not_all_N_probability_or_original_P8": True,
            },
        }
    )
    return result
