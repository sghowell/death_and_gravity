"""Exact conditional theorem margin; quantum and source inputs stay distinct."""

import sympy as sp
from p8a_maxwell.domain import nonnegative, rational

from . import envelope, history


def theorem_constants(ratio, contraction, past_caps, future_caps, *, beta_m, delta, sigma):
    """delta bounds kappa*hbar/(8*pi²*tau²); sigma is the signed-source loss.

Nonpositive returned index_upper forces focusing if the normal reaches
tau. This does not verify the conditional geometric hypotheses for a model.
    """
    h = nonnegative(contraction, positive=True)
    delta, sigma = map(nonnegative, (delta, sigma))
    data = envelope.coefficients(past_caps, future_caps, ratio, beta_m=beta_m)
    if min(data["past_caps"][0], data["future_caps"][0]) < h:
        raise ValueError("the contraction history contradicts a stated endpoint Hubble cap")
    gain = history.lower(h, data["ratio"])
    gradient = sp.Rational(18, 5)
    quantum, source = delta*data["total_cost"], sigma*data["source_weight"]
    margin = gain-gradient-quantum-source
    return {"envelope": data, "contraction": h, "delta": delta, "sigma": sigma,
            "history_gain": gain, "future_gradient": gradient,
            "quantum_cost": quantum, "source_cost": source,
            "strict_focusing_margin": margin, "index_plus_K_times_tau_upper": -margin,
            "sufficient_focusing_test": bool(margin >= 0),
            "actual_SEE_or_geometric_hypotheses_verified": False}


def physical_enclosure(*, kappa, hbar, tau, cosmological_constant, other_eed_lower):
    """Exact rational upper delta uses pi>3; no hidden source or unit choice."""
    kappa, tau = (nonnegative(value, positive=True) for value in (kappa, tau))
    hbar = nonnegative(hbar)
    lam, lower = map(rational, (cosmological_constant, other_eed_lower))
    penalty = max(lam, 0)+kappa*max(-lower, 0)
    return {"exact_delta": kappa*hbar/(8*sp.pi**2*tau**2),
            "rational_delta_upper": kappa*hbar/(72*tau**2),
            "sigma": tau**2*penalty, "separate_source_penalty": penalty,
            "Lambda": lam, "other_EED_lower": lower}


def calibration():
    caps = (2, 4, 16, 96)
    kwargs = {"ratio": sp.Rational(1, 100), "contraction": sp.Rational(3, 2),
              "past_caps": caps, "future_caps": caps, "beta_m": 1,
              "delta": sp.Rational(1, 10**8)}
    base = theorem_constants(**kwargs, sigma=0)
    sourced = theorem_constants(**kwargs, sigma=1)
    return {"zero_source": base, "source_sigma_at_most_one": sourced,
            "beta_absolute_cap": sp.Integer(1), "cost_coarsening": sp.Integer(13000000),
            "cost_coarsening_gap": 13000000-base["envelope"]["total_cost"],
            "zero_source_margin_above_three_quarters": base["strict_focusing_margin"]-sp.Rational(3, 4),
            "sourced_margin_above_one_third": sourced["strict_focusing_margin"]-sp.Rational(1, 3),
            "initial_pointwise_SEC_assumed": False,
            "short_contraction_history_is_a_geometric_hypothesis": True,
            "actual_cosmology_or_radiation_witness": False}


def identities():
    kappa, hbar, tau, g, planck = sp.symbols("kappa hbar tau G ellP_squared", positive=True)
    delta = kappa*hbar/(8*sp.pi**2*tau**2)
    return {"Planck_dimensionless_dictionary": sp.simplify(
                delta.subs(kappa, 8*sp.pi*g).subs(g*hbar, planck)-planck/(sp.pi*tau**2)),
            "source_weight_both_sampler_pieces": history.moments()["zeroth"]*(1+sp.Rational(1, 100))-sp.Rational(1313, 3500),
            "history_threshold_numerator": history.lower(sp.Rational(3, 2), sp.Rational(1, 100))-sp.Rational(63351, 14000),
            "future_gradient_cubic": 3*history.moments()["first"]-sp.Rational(18, 5)}
