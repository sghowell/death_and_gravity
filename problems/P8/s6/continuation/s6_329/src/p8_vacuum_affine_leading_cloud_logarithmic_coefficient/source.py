"""Fixed original Born leading measure; no new radiative emission dynamics."""

from functools import cache

from p8_vacuum_affine_continuous_logarithmic_coefficient import source as previous

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
            "unchanged_original_action_and_parameters": True,
            "frozen_Born_angular_measure_not_state_dependent_branching": True,
            "D4_sharp_cutoff_not_new_dimensional_conversion": True,
            "leading_reference_not_full_interacting_probability": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_defined_reference": "Use the S299 D4 Born helicity-summed intensity Lambda(dn)*dw/w for0<w<=1, with integral Lambda=a=K0/(4pi^2 kappa). The Born hard state and intensity remain fixed. The S328 coefficient is a mark depending on the resulting angular-energy measure and its physical recoil, not a replacement of the emission intensity.",
    }
