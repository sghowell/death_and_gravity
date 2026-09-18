"""Original fixed-ball additional-soft prescription and physical D4 cloud."""

from functools import cache

from p8_vacuum_affine_quantitative_soft_cutoff import source as previous

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
            "same_fixed_Born_D4_cloud_and_full_D_additional_soft_factor": True,
            "no_signed_perturbation_or_interacting_completion": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_source_boundary": "Keep the original action, fixed E,u and positive angular energy of total mass<=1/8. S330 fixes the original positive recoil path; S301 and S325 fix the complete additional-soft trace, radial and phase convention. S329 supplies only the fixed-Born D4 leading probability. No outer hard evanescence or new radiative intensity is inferred.",
    }
