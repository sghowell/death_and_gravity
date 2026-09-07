"""Generic finite-beta affine cost and robust macroscopic sufficient margin."""

import sympy as sp
from p8a_maxwell.domain import nonnegative, rational
from p8a_maxwell_focusing import envelope, history

from . import geometry

COST_UPPER = sp.Integer(13000000)
WEIGHTED_DELTA_MAX = sp.Rational(1, 10**8)
SIGMA_MAX = sp.Integer(5)


def affine_cost(past_caps=geometry.EXPANDED_CAPS, future_caps=geometry.FUTURE_CAPS):
    """C(beta)=C0+Cbeta*abs(beta); no new finite prescription is selected."""
    zero = envelope.coefficients(past_caps, future_caps, geometry.RATIO, beta_m=0)
    one = envelope.coefficients(past_caps, future_caps, geometry.RATIO, beta_m=1)
    c0, cb = zero["total_cost"], one["total_cost"]-zero["total_cost"]
    return {"past_caps": zero["past_caps"], "future_caps": zero["future_caps"], "C0": c0, "Cbeta": cb,
            "C_at_abs_beta_one": one["total_cost"],
            "C0_upper_margin": COST_UPPER-c0, "Cbeta_upper_margin": COST_UPPER-cb,
            "past_V0": zero["past_reference_loss"],
            "past_Vbeta": one["past_reference_loss"]-zero["past_reference_loss"],
            "future_V0": zero["future_reference_loss"],
            "future_Vbeta": one["future_reference_loss"]-zero["future_reference_loss"]}


def cost_at_beta(beta_m):
    beta = rational(beta_m)
    data = affine_cost()
    return data["C0"]+data["Cbeta"]*abs(beta)


def theorem_gate(weighted_delta, sigma):
    """weighted_delta bounds actual delta*(1+abs(beta_M)), not delta alone."""
    weighted_delta, sigma = map(nonnegative, (weighted_delta, sigma))
    if weighted_delta > WEIGHTED_DELTA_MAX or sigma > SIGMA_MAX:
        raise ValueError("the input exceeds the named robust cosmological gate")
    gain = history.lower(geometry.CONTRACTION, geometry.RATIO)
    gradient = sp.Rational(18, 5)
    weight = sp.Rational(13, 35)*(1+geometry.RATIO)
    quantum, source = COST_UPPER*weighted_delta, weight*sigma
    margin = gain-gradient-quantum-source
    return {"weighted_delta_upper": weighted_delta, "sigma_upper": sigma,
            "history_gain": gain, "future_gradient": gradient, "source_weight": weight,
            "quantum_cost_upper": quantum, "source_cost_upper": source,
            "strict_focusing_margin_lower": margin,
            "margin_above_one_eighth": margin-sp.Rational(1, 8),
            "future_relative_caps_or_actual_SEE_verified": False}


def identities():
    b, delta, c0, cb = sp.symbols("b delta C0 Cbeta", nonnegative=True)
    bound = COST_UPPER
    return {"finite_beta_joint_gate": sp.expand(
                bound*delta*(1+b)-delta*(c0+cb*b)
                -delta*((bound-c0)+(bound-cb)*b)),
            "robust_history_gain": history.lower(geometry.CONTRACTION, geometry.RATIO)-sp.Rational(2009079, 350000),
            "universal_final_margin": theorem_gate(WEIGHTED_DELTA_MAX, SIGMA_MAX)["strict_focusing_margin_lower"]-sp.Rational(47079, 350000),
            "universal_margin_excess": theorem_gate(WEIGHTED_DELTA_MAX, SIGMA_MAX)["margin_above_one_eighth"]-sp.Rational(3329, 350000)}


def calibration():
    return {"expanded_affine_cost": affine_cost(),
            "tight_reference_regression": affine_cost(geometry.REFERENCE_CAPS, geometry.REFERENCE_CAPS),
            "tight_regression_is_not_used_as_future_geometry": True,
            "universal_gate": theorem_gate(WEIGHTED_DELTA_MAX, SIGMA_MAX),
            "uniform_bound_on_abs_beta_required": False,
            "generic_finite_beta_remains_explicit": True,
            "COST_UPPER": COST_UPPER}
