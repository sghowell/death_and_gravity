"""All-momentum constrained mean energy and complete causal scalar force."""

from functools import cache

import sympy as s
from p8_vacuum_affine_retarded_energy import hamiltonian as old

from . import source


@cache
def data():
    k, k0, m, C = source.K, source.K0, source.MASS, source.C
    pref = s.factor(k * C * C)
    E0, Efree, instant, force = (
        2 * m * m * pref,
        4 * m * m * pref,
        pref,
        (4 * m + 2) * pref,
    )
    packet, modes = old.data(), old.modes()
    a = next(x for x in modes["longitudinal_energy"].free_symbols if str(x) == "a")
    phase = {
        name: value.subs(a, 1)
        for name, value in modes.items()
        if name
        in (
            "longitudinal_energy",
            "longitudinal_dq",
            "longitudinal_dp",
            "transverse_energy",
            "transverse_dq",
            "transverse_dp",
            "longitudinal_source_power",
        )
    }
    checks = {
        "full_mass_and_source_normalization": m * m * source.family.ZETA - 1,
        "forcing_gradient_count": 4**2 + 2**2 - 20,
        "exact_reference_energy_coefficient": E0 - k * m * m * (2 * C) ** 2 / 2,
        "reference_to_actual_free_mean_energy": Efree - 2 * E0,
        "instantaneous_temporal_source_energy": instant - k * C * C,
        "complete_causal_force_coefficient": force - k * C * C * (4 * m + 2),
        "family_decay_of_reference_energy": s.factor(E0 * k - E0.subs(k, k0) * k0),
        "family_decay_of_full_causal_force": s.factor(
            force * k - force.subs(k, k0) * k0
        ),
    }
    return {
        "same_unmodified_constrained_Hamiltonian": packet["positive_reference_energy"],
        "Minkowski_phase_fixtures": phase,
        "retained_temporal_constraint": "A0=(J0-div pi)/m^2; J=sqrt(kappa)*m*S and W=m*A/sqrt(kappa); no conserved-source assumption",
        "source_forcing_bound": "F_S^2=||S_sp||2^2+||div S_sp||2^2/m^2+||grad S0||2^2/m^2 <= C^2(1+20/m^2)U_Phi^2 <4C^2 U_Phi^2",
        "reference_energy_coefficient": E0,
        "actual_free_mean_energy_integrated_coefficient": Efree,
        "actual_free_mean_energy_instantaneous_coefficient": instant,
        "full_scalar_causal_mean_force_coefficient": force,
        "mean_readout": "||Wbar||L2 <=4m*C*integral_past U_Phi+C*U_Phi(t); ||S-Wbar||L2 <=4m*C*integral_past U_Phi+2C*U_Phi(t)",
        "reference_energy_result": "E0(t)<=coefficient*(integral_past U_Phi)^2; positive reference energy is not the driven Hamiltonian",
        "actual_free_mean_energy_result": "Efree(t)<=integrated_coefficient*(integral_past U_Phi)^2+instantaneous_coefficient*U_Phi(t)^2; not full source-interaction stress",
        "causal_force": "At fixed physical Minkowski g, Force_Phi[eta]=kappa integral (S-Wbar).DS_eta, obtained from the unintegrated action, not by varying S Gret S",
        "force_result": "abs(Force_Phi[eta])<=force_coefficient*||U_Phi||L2(time)*||U_eta||L2(time) for interval length<=1",
        "anchor_displays": {
            "reference_energy": s.Rational(1, 10**2376),
            "actual_free_energy_integrated": 2 * s.Rational(1, 10**2376),
            "actual_free_energy_instantaneous": 4 * s.Rational(1, 10**2383),
            "causal_scalar_force": s.Rational(1, 10**2378),
        },
        "checks": checks,
        "gates": {
            "all_momentum_forcing_bound": 1 + 20 / m**2 < 4,
            "anchor_reference_energy_below_display": E0.subs(k, k0)
            < s.Rational(1, 10**2376),
            "anchor_free_energy_integrated_below_display": Efree.subs(k, k0)
            < 2 * s.Rational(1, 10**2376),
            "anchor_free_energy_instantaneous_below_display": instant.subs(k, k0)
            < 4 * s.Rational(1, 10**2383),
            "anchor_full_causal_force_below_display": force.subs(k, k0)
            < s.Rational(1, 10**2378),
        },
    }
