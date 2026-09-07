"""Exact proper-clock, source and Planck dictionaries without observed inputs."""

import sympy as sp
from p8a_maxwell.domain import nonnegative, rational
from p8a_maxwell_focusing import focusing

from . import calibration, geometry


def proper_scales(tau, *, anchored=False):
    tau = nonnegative(tau, positive=True)
    tube = geometry.neighborhood(anchored=anchored)
    return {"tau": tau, "actual_short_history_duration": geometry.RATIO*tau,
            "reference_Hstar": 2/tau,
            "reference_power_law_age_range": [tau/4, tau/3],
            "actual_observed_H0_range": [value/tau for value in tube["observed_H0_times_tau_range"]],
            "anchored_at_observer": anchored,
            "reference_power_law_age_is_actual_quantum_age": False}


def source_budget(*, kappa, hbar, tau, beta_m, cosmological_constant, other_eed_lower):
    beta = rational(beta_m)
    data = focusing.physical_enclosure(kappa=kappa, hbar=hbar, tau=tau,
                                       cosmological_constant=cosmological_constant,
                                       other_eed_lower=other_eed_lower)
    joint = data["rational_delta_upper"]*(1+abs(beta))
    return {**data, "beta_M": beta,
            "exact_weighted_delta": data["exact_delta"]*(1+abs(beta)),
            "rational_weighted_delta_upper": joint,
            "sufficient_quantum_budget": bool(joint <= calibration.WEIGHTED_DELTA_MAX),
            "sufficient_source_budget": bool(data["sigma"] <= calibration.SIGMA_MAX),
            "actual_field_state_source_or_metric_verified": False}


def identities():
    tau, g, hbar, c, seconds, lp2, tp2 = sp.symbols("tau G hbar c tau_seconds ellP2 tP2", positive=True)
    rho, pressure, hubble = sp.symbols("rho pressure Hstar", real=True)
    lam = sp.Symbol("Lambda", real=True)
    delta = g*hbar/(sp.pi*tau*tau)
    return {"proper_Planck_dictionary_c_one": sp.simplify(delta.subs(g*hbar, lp2)-lp2/(sp.pi*tau**2)),
            "seconds_Planck_dictionary": sp.simplify((lp2/(sp.pi*tau**2)).subs({tau: c*seconds, lp2: c*c*tp2})-tp2/(sp.pi*seconds**2)),
            "perfect_fluid_EED": sp.expand(rho-(rho-3*pressure)/2-(rho+3*pressure)/2),
            "ordinary_radiation_EED": sp.simplify(((rho+3*pressure)/2).subs(pressure, rho/3)-rho),
            "ordinary_dust_EED": sp.simplify(((rho+3*pressure)/2).subs(pressure, 0)-rho/2),
            "reference_Lambda_fraction": sp.simplify((lam/(3*hubble**2)).subs(hubble, 2/tau)-lam*tau*tau/12)}


def calibration_data():
    cases = {}
    for lam, lower in ((1, 2), (1, -2), (-1, 2), (-1, -2)):
        # Only source arithmetic, not units for a physical universe.
        cases[f"Lambda={lam},other_lower={lower}"] = source_budget(
            kappa=2, hbar=0, tau=1, beta_m=-7,
            cosmological_constant=lam, other_eed_lower=lower)
    return {"anchored_scales_tau_one": proper_scales(1, anchored=True),
            "unanchored_scales_tau_one": proper_scales(1),
            "Lambda_positive_fraction_of_reference_3Hstar_squared_upper": sp.Rational(5, 12),
            "Lambda_positive_fraction_of_actual_3H0_squared_upper_unanchored": sp.Rational(50000, 118803),
            "signed_source_controls": cases,
            "ordinary_matter_nonnegative_EED_is_separate_assumption": True,
            "pure_Maxwell_specialization_has_other_EED_zero": True,
            "observed_cosmological_parameters_used": False}
