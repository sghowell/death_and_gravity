"""Conditional Einstein/extra-source dictionary, not a focusing example."""

import sympy as sp

from .bounds import qsei_coefficients
from .domain import nonnegative, rational


def see_ricci(maxwell_eed, other_eed, kappa, cosmological_constant):
    """G_FK + Lambda*g_FK = -kappa*(T_Maxwell+T_other)."""
    maxwell, other, kappa, lam = map(sp.sympify, (maxwell_eed, other_eed, kappa, cosmological_constant))
    return lam-kappa*(maxwell+other)


def geometric_constants(h0, h1, h2, h3, duration, *, beta_m, kappa, hbar,
                        cosmological_constant, other_eed_lower):
    """Conditional uniform constants for every permitted normal-geodesic window.

other_eed_lower is a signed pointwise lower bound, not an assertion of
SEC for arbitrary extra matter. Positive Lambda and negative extra EED
are retained separately; beneficial signs are dropped only upwards.
    """
    kappa = nonnegative(kappa, positive=True)
    lam, lower = map(rational, (cosmological_constant, other_eed_lower))
    field = qsei_coefficients(h0, h1, h2, h3, duration, beta_m=beta_m, hbar=hbar)
    source = max(lam, 0)+kappa*max(-lower, 0)
    return {"Q2": kappa*field["Q2"], "Q0": kappa*field["Q0"]+source,
            "raw_constant_before_nonnegative_coarsening": kappa*field["Q0"]+lam-kappa*lower,
            "separate_source_penalty": source,
            "cosmological_constant": lam, "additional_EED_lower": lower,
            "new_SEE_solution_asserted": False}


def identities():
    rho, pressure, other, kappa, lam = sp.symbols("rho p E_other kappa Lambda", real=True)
    trace = rho-3*pressure
    # Trace of G+Lambda*g=-kappa*T and then its 00 component.
    scalar = kappa*trace+4*lam
    ricci = -kappa*rho-lam+scalar/2
    gamma, curvature = sp.symbols("gamma E_I", real=True)
    radiation, scale = sp.symbols("C a", positive=True)
    return {
        "cosmological_Ricci_sign": sp.expand(ricci-see_ricci((rho+3*pressure)/2, 0, kappa, lam)),
        "additional_matter_EED_not_omitted": sp.expand(see_ricci(rho, other, kappa, lam)
                                                     -see_ricci(rho, 0, kappa, lam)+kappa*other),
        "radiation_EED_is_its_positive_density": sp.simplify((radiation/scale**4+3*radiation/(3*scale**4))/2-radiation/scale**4),
        "explicit_gravitational_curvature_source_if_moved": sp.expand(
            -kappa*(rho+gamma*curvature/kappa)-(-kappa*rho-gamma*curvature)),
        "frozen_A1_gradient_floor": sp.Rational(3, 1)/sp.Rational(3, 4)-4,
    }


def calibration():
    # All four source signs are tested; this is not a chosen cosmology.
    cases = {}
    for lam, lower in ((1, 2), (1, -2), (-1, 2), (-1, -2)):
        cases[f"Lambda={lam},other_lower={lower}"] = geometric_constants(
            1, 1, 1, 1, 1, beta_m=1, kappa=2, hbar=1,
            cosmological_constant=lam, other_eed_lower=lower)
    return {"signed_source_controls": cases,
            "Hmax_times_tau_at_most_one": {"K_trace_times_tau_cap": sp.Integer(3),
                "A1_zeta_zero_gradient_lower": sp.Integer(4), "strict_incompatibility_gap": sp.Integer(1),
                "small_Q2_implies_focusing": False},
            "A1_initial_pointwise_Ricci_premise_removed": False,
            "A1_extended_domain_and_K_threshold_verified": False}
