"""Original action and named logarithmic coefficient; no hard-loop completion."""

from functools import cache

from p8_vacuum_affine_known_hard_soft_subtraction import source as previous

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
            "complete_named_coefficient_not_individual_current": True,
            "finite_positive_angular_energy_not_quantum_state": True,
            "frozen_S326_rejecting_API_and_sources_not_rewritten": True,
        },
        "whole_original_parameters": original_parameters(),
        "whole_source_boundary": "The full named S326 logarithmic coefficient is reorganized exactly on conserved on-shell radiative states. Its individual null current can lack a collinear point limit while the complete coefficient has a continuous extension. The fixed original scalar/Einstein action, finite hard matching, all-N hard amplitudes and original P8 frontier are unchanged.",
    }
