"""Conditional exact curvature-envelope QSEI constants, beta_M not frozen."""

import sympy as sp

from .domain import nonnegative, rational


def reference_loss_polynomial(h0, h1, h2, h3, beta_abs):
    """Magnitude of E_conf in hbar/(2880pi²) units from absolute jet caps."""
    h0, h1, h2, h3, beta = map(sp.sympify, (h0, h1, h2, h3, beta_abs))
    return sp.expand(186*(h0**4+2*h0**2*h1)
                     +beta*(18*h3+90*h1**2+90*h0*h2+108*h0**2*h1))


def sampler_factor(h0, h1, duration):
    """Use pi>3 in two proper Dirichlet Poincare estimates."""
    h0, h1, duration = map(sp.sympify, (h0, h1, duration))
    return sp.expand(1+2*h0*duration/3+(3*h0**2/4+3*h1/2)*duration**2/9)


def envelope(h0, h1, h2, h3, duration, *, beta_m):
    """Caps are hypotheses, not inferred from an uninspected metric.

The beta_m argument is mandatory even if the caller intends the named
zero-type-D specialization. D is an enclosing proper-time duration.
    """
    h0, h1, h2, h3 = map(nonnegative, (h0, h1, h2, h3))
    duration = nonnegative(duration, positive=True)
    beta = rational(beta_m)
    factor = sampler_factor(h0, h1, duration)
    loss = reference_loss_polynomial(h0, h1, h2, h3, abs(beta))
    return {"hubble_jet_caps": [h0, h1, h2, h3], "proper_duration": duration,
            "beta_M": beta, "sampler_factor": factor,
            "reference_loss_numerator": loss,
            "derivative_coefficient_pi_squared_over_hbar": factor**2/8,
            "zero_coefficient_pi_squared_over_hbar": loss/2880}


def qsei_coefficients(h0, h1, h2, h3, duration, *, beta_m, hbar):
    data = envelope(h0, h1, h2, h3, duration, beta_m=beta_m)
    hbar = nonnegative(hbar)
    return {"Q2": hbar*data["derivative_coefficient_pi_squared_over_hbar"]/sp.pi**2,
            "Q0": hbar*data["zero_coefficient_pi_squared_over_hbar"]/sp.pi**2}


def scaling_identities():
    tau = sp.Symbol("tau", positive=True)
    c0, c1, c2, c3, beta = sp.symbols("c0 c1 c2 c3 beta_abs", nonnegative=True)
    factor = sampler_factor(c0/tau, c1/tau**2, tau)
    loss = reference_loss_polynomial(c0/tau, c1/tau**2, c2/tau**3, c3/tau**4, beta)
    return {"dimensionless_sampler_factor": sp.simplify(factor-sampler_factor(c0, c1, 1)),
            "curvature_loss_fourth_power": sp.simplify(tau**4*loss-reference_loss_polynomial(c0, c1, c2, c3, beta))}


def calibration():
    # An arithmetic control only: these caps and Planck ratios do not give
    # a physical background, field prescription or focusing example.
    positive = envelope(1, 1, 1, 1, 1, beta_m=1)
    negative = envelope(1, 1, 1, 1, 1, beta_m=-1)
    zero = envelope(1, 1, 1, 1, 1, beta_m=0)
    return {"arithmetic_plus_beta": positive, "arithmetic_minus_beta": negative,
            "named_zero_type_D_arithmetic": zero,
            "physical_calibration": False,
            "beta_M_fixed_by_scalar_gamma": False}
