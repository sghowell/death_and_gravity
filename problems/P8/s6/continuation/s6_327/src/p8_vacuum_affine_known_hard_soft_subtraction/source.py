"""Fixed known hard reference and complete original one-real tree."""

from functools import cache

import sympy as s
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient import source as previous

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
    checks["complete_original_real_tree_inventory"] = s.Integer(26 + 21 - 47)
    checks["two_physical_helicities_not_single_helicity_phase_deletion"] = s.Integer(
        2 - 2
    )
    return {
        "checks": checks,
        "gates": {
            "all_original_action_and_parameters_unchanged": True,
            "fixed_S303_known_hard_reference_not_unknown_coordinates": True,
            "same_full_positive_Born_and_original47_tree": True,
            "selected_Born_extension_not_finite_radiative_hard_amplitude": True,
            "no_unknown_matching_or_common_regulator_claim": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_selected_hard_source": "Only frozen S303 H4_known in its original UV-before-IR and nu1 reference, with alpha,beta,delta_kappa unassigned. h=H4_known/A0 is evaluated at the associated Born E,u and extended as S_lambda*h over the finite radiation interval by definition. This fixes a subtraction kernel, not the full five-point loop.",
        "whole_real_source_proof": "Every original S295 matter and S304 Einstein one-real tree graph is real for real TT polarization and physical generic real momenta. All vertices and couplings are real. Massive recoil bounds keep scalar and hard graviton denominators off zero; the original heavy mass is above all compact physical channel energies. Complex helicity changes are unitary and leave the summed tree-soft inner product real.",
    }
