"""Original positive energy increments and fixed-Born D4 cutoff reference."""

from functools import cache

from p8_vacuum_affine_leading_cloud_logarithmic_coefficient import source as previous

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
    return {
        "checks": dict(previous.data()["checks"]),
        "gates": {
            "original_action_parameters_and_recoil_unchanged": True,
            "positive_added_energy_with_fixed_E_and_u": True,
            "same_fixed_Born_D4_leading_reference": True,
            "no_signed_perturbation_or_interacting_completion": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_source_boundary": "Keep E,u fixed and sigma=tau+rho positive with mass sigma<=1/8. The same original S300 recoil and full S328 logarithmic coefficient are used. Cutoff comparisons use only S329 fixed-Born D4 leading probability, not a new radiative intensity or dimensional matching.",
    }
