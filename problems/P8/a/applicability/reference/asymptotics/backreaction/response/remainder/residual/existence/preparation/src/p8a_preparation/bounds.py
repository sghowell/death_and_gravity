"""Exact full-map bounds for the named, conserved-source preparation.

These are bounds for the complete integrated-trace map, not only its
nonlinear mode term.  The fixed physical amplitude is delta=10^-14 and
the future dimensionless conformal interval has length 10^-10.  The
geometric and state hypotheses are derived in notes/contraction.md.

The lower-level arithmetic interfaces are conditional on their displayed
geometric caps; passing numbers does not establish their physical origin.
The calibration obtains that origin from the immutable A.8 constants.
"""

from functools import cache

import sympy as sp
from p8a_existence.mode_lipschitz import rational
from p8a_residual import bounds as prior_bounds

DELTA = sp.Rational(1, 10**14)
LENGTH = sp.Rational(1, 10**10)
RADIUS = sp.Rational(1, 10**6)
POTENTIAL_CAP = sp.Rational(2, 10**12)
BACKGROUND_DERIVATIVE_CAP = sp.Rational(1, 10**10)
HISTORY_DERIVATIVE_CAP = sp.Rational(1, 10**5)
HISTORY = sp.Integer(3)
SOURCE_CAP = sp.Rational(12, 10**9)
BACKGROUND_Q_CAP = sp.Rational(1, 10**8)
AUXILIARY_Q_CAP = sp.Rational(2, 10**8)
AUXILIARY_P_CAP = sp.Rational(1, 10**6)
LOCAL_CAP = sp.Rational(17, 30)
INVERSE_CAP = sp.Rational(21, 100)
CONTRACTION_CAP = sp.Rational(3, 25)
CENTER_RHS_CAP = sp.Rational(3, 10**7)


def _nonnegative(value, name):
    result = rational(value)
    if result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _positive(value, name):
    result = rational(value)
    if result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def center_bounds(length, potential_cap, source_cap):
    """Bound the zero-perturbation RHS after the source is switched.

Here |c-cbar|<=source_cap, c(0)=cbar(0), a>=2 and
1/3<=h<=1/2. No bound on a derivative of the prescribed cutoff is used.
    """
    length = _positive(length, "length")
    potential_cap = _nonnegative(potential_cap, "potential_cap")
    source_cap = _nonnegative(source_cap, "source_cap")
    p = 8*source_cap*(3+length*(1+9*potential_cap))
    q = length*p/4
    rhs = p*(1+sp.Rational(9, 4)*length)
    return {"P_difference": p, "q_difference": q, "rhs": rhs}


def pair_bounds(delta, length, radius, potential_cap, source_cap,
                background_derivative_cap, q_cap, p_cap):
    """Full pre-inverse Lipschitz bound per ||X1-X2||_infinity.

Requires both metric solutions to share initial a,h, to obey
2<=a<=3, 1/3<=h<=1/2, and |u|<=potential_cap.  The source c is
the same prescribed function in both solutions.  The q and P=a^2 q'
caps are auxiliary-solution caps, not state-independent quantum bounds.
The mode-response term is added separately by calibration().
    """
    delta = _positive(delta, "delta")
    length = _positive(length, "length")
    radius = _positive(radius, "radius")
    potential_cap = _nonnegative(potential_cap, "potential_cap")
    source_cap = _nonnegative(source_cap, "source_cap")
    background_derivative_cap = _nonnegative(
        background_derivative_cap, "background_derivative_cap")
    q_cap = _nonnegative(q_cap, "q_cap")
    p_cap = _nonnegative(p_cap, "p_cap")
    ell, b, c = length, potential_cap, source_cap
    # h1-h2 solves a damped scalar equation, hence the three gains below.
    metric = {"u": ell, "h": ell**2/2, "log_a": ell**3/6,
              "a_squared": 3*ell**3, "inverse_a_squared": ell**3/12,
              "d": ell**3/12, "dprime": ell**2/4}
    einstein = 3*ell**2/(40*delta)+b*ell**4/(80*delta)
    curvature = (b*ell**2/4+ell**2/240
                 +(b+sp.Rational(1, 2))*ell**3/180)
    forcing = 72*c*ell**2+72*c*b*ell**3
    p = einstein+curvature+forcing
    q = ell*p/4+p_cap*ell**4/48
    auxiliary = p+9*q+9*q_cap*ell**2+3*q_cap*ell**3
    local = (LOCAL_CAP+ell/4+b*ell**2/4
             +(radius+background_derivative_cap)*ell**3/12
             +radius*ell**3/4)
    return {"metric_gains": metric, "Einstein_P": einstein,
            "curvature_P": curvature, "source_P": forcing,
            "P": p, "q": q, "auxiliary_Sprime": auxiliary,
            "local_Wick_terms": local,
            "rhs_without_mode_response": auxiliary+local}


def contraction_gate(inverse_cap, rhs_lipschitz, center_rhs, radius):
    """Check the ordinary Banach inequalities for supplied proved bounds."""
    inverse_cap = _nonnegative(inverse_cap, "inverse_cap")
    rhs_lipschitz = _nonnegative(rhs_lipschitz, "rhs_lipschitz")
    center_rhs = _nonnegative(center_rhs, "center_rhs")
    radius = _positive(radius, "radius")
    contraction = inverse_cap*rhs_lipschitz
    center_image = inverse_cap*center_rhs
    self_map = center_image+contraction*radius
    if contraction >= 1:
        raise ValueError("a strict full-map contraction is required")
    if self_map > radius:
        raise ValueError("the proposed ball is not proved invariant")
    return {"contraction": contraction, "center_image": center_image,
            "self_map": self_map, "ball_margin": radius-self_map,
            "fixed_point_distance": center_image/(1-contraction)}


def inverse_bound():
    """Bound J_EulerGamma/(4*pi^2) on the fixed interval LENGTH.

EulerGamma>0, pi>3, log(10)>2 and sqrt(10^-10)=10^-5 give
2 L/(1-L)+(2/9)*sqrt(L)+1/5.  The pole is included.
    """
    return 2*LENGTH/(1-LENGTH)+sp.Rational(2, 9*10**5)+sp.Rational(1, 5)


def elementary_margins():
    """Rational witnesses for the elementary constants used in the proof."""
    pi_lower = 4*sum((sp.Rational((-1)**j, 2*j+1)
                      for j in range(8)), sp.S.Zero)
    y0 = sp.Rational(5, 2)
    a0_fourth = y0**4-DELTA
    h0_fourth = y0**8/a0_fourth**3
    return {
        "pi_integral_lower_minus_three": pi_lower-3,
        "e_squared_upper_below_ten": sp.Integer(10-3**2),
        "exp_half_partial_sum_above_three_halves": sp.Rational(1, 8),
        "exp_three_partial_sum_above_ten": sp.Integer(3),
        "exp_two_partial_sum_above_three": sp.Integer(2),
        "inverse_cap_margin": INVERSE_CAP-inverse_bound(),
        "h_lower_first_exit_margin": sp.Rational(3, 8)-LENGTH/3-sp.Rational(1, 3),
        "h_upper_first_exit_margin": sp.Rational(1, 2)-sp.Rational(5, 12)-LENGTH/3,
        "h_derivative_margin": sp.Rational(1, 3)-sp.Rational(1, 4)-POTENTIAL_CAP,
        "a_upper_first_exit_margin": 3-sp.Rational(13, 5)*sp.Rational(20, 19),
        "initial_a_fourth_lower": a0_fourth-sp.Rational(12, 5)**4,
        "initial_a_fourth_upper": sp.Rational(13, 5)**4-a0_fourth,
        "initial_h_fourth_lower": h0_fourth-sp.Rational(3, 8)**4,
        "initial_h_fourth_upper": sp.Rational(5, 12)**4-h0_fourth,
    }


@cache
def calibration():
    """Reproduce the complete named finite-interval self-map certificate."""
    old = prior_bounds.calibration()
    jets = old["potential_bounds"]
    s0, s1 = (DELTA*old["effective_S_Born_norm_constants"][j]
              +DELTA**2*old["mode_remainder_constants"][j] for j in range(2))
    source = 81*DELTA*(14000+DELTA*8*10**14+sp.Rational(3, 2880*1922))
    input_margins = {
        "source": SOURCE_CAP-source,
        "future_potential": POTENTIAL_CAP-DELTA*jets[0]-LENGTH*RADIUS,
        "background_derivative": BACKGROUND_DERIVATIVE_CAP-DELTA*jets[1],
        "full_history_derivative": HISTORY_DERIVATIVE_CAP-DELTA*jets[1]-RADIUS,
        "background_q": BACKGROUND_Q_CAP-s0/4,
        "background_qprime": BACKGROUND_Q_CAP-(s1+s0)/4,
        "background_observation_stays_on_plateau": 3-sp.Rational(5, 2)-LENGTH,
    }
    center = center_bounds(LENGTH, POTENTIAL_CAP, SOURCE_CAP)
    pairs = pair_bounds(DELTA, LENGTH, RADIUS, POTENTIAL_CAP, SOURCE_CAP,
                        BACKGROUND_DERIVATIVE_CAP, AUXILIARY_Q_CAP, AUXILIARY_P_CAP)
    p_derived = 9*BACKGROUND_Q_CAP+center["P_difference"]+pairs["P"]*RADIUS
    q_derived = BACKGROUND_Q_CAP+LENGTH*AUXILIARY_P_CAP/4
    # A.10 shared-history lemma.  log(3*10^10)<32 and exp(2z)<2.
    # Both estimates are deliberately coarsened to exact rationals.
    z = HISTORY_DERIVATIVE_CAP*HISTORY**3
    amplitude = HISTORY_DERIVATIVE_CAP*HISTORY
    response_quadratic = amplitude*LENGTH**2*(sp.Rational(5, 4)+16)
    response_higher = 36*HISTORY_DERIVATIVE_CAP**2*HISTORY**5*LENGTH
    response = response_quadratic+response_higher
    rhs = pairs["rhs_without_mode_response"]+response
    actual = contraction_gate(inverse_bound(), rhs, center["rhs"], RADIUS)
    rounded_center_image = INVERSE_CAP*CENTER_RHS_CAP
    rounded_self_map = rounded_center_image+CONTRACTION_CAP*RADIUS
    margins = {
        **input_margins,
        "auxiliary_P": AUXILIARY_P_CAP-p_derived,
        "auxiliary_q": AUXILIARY_Q_CAP-q_derived,
        "response_exponential_argument": sp.Rational(1, 2)-2*z,
        "A10_rational_strength": sp.Rational(1, 4)-z,
        "rounded_center_rhs": CENTER_RHS_CAP-center["rhs"],
        "rounded_rhs_lipschitz": sp.Rational(567, 1000)-rhs,
        "rounded_contraction": CONTRACTION_CAP-INVERSE_CAP*rhs,
        "rounded_ball": RADIUS-rounded_self_map,
        **elementary_margins(),
    }
    if any(value <= 0 for value in margins.values()):
        raise ValueError("a full-map input, geometric or strict rounding margin failed")
    return {
        "delta": DELTA, "length": LENGTH, "radius": RADIUS,
        "initial_background_label": sp.Rational(5, 2),
        "geometry": {"a_min": sp.Integer(2), "a_max": sp.Integer(3),
                     "h_min": sp.Rational(1, 3), "h_max": sp.Rational(1, 2),
                     "u_cap": POTENTIAL_CAP, "background_uprime_cap": BACKGROUND_DERIVATIVE_CAP,
                     "history_uprime_cap": HISTORY_DERIVATIVE_CAP, "history": HISTORY},
        "A8_inputs": {"S0": s0, "S1": s1, "source_derived": source,
                      "U0_per_delta": jets[0], "U1_per_delta": jets[1]},
        "source_cap": SOURCE_CAP, "background_q_and_qprime_cap": BACKGROUND_Q_CAP,
        "auxiliary_P_cap": AUXILIARY_P_CAP, "auxiliary_q_cap": AUXILIARY_Q_CAP,
        "auxiliary_P_derived": p_derived, "auxiliary_q_derived": q_derived,
        "center": center, "pairs": pairs,
        "mode_response": {"strength": z, "log_history_ratio_upper": sp.Integer(32),
                          "exponential_upper": sp.Integer(2),
                          "quadratic": response_quadratic, "higher": response_higher,
                          "total": response},
        "local_cap": LOCAL_CAP, "rhs_lipschitz": rhs,
        "inverse_bound": inverse_bound(), "inverse_cap": INVERSE_CAP,
        "actual_gate": actual,
        "rounded": {"center_rhs": CENTER_RHS_CAP,
                    "center_image": rounded_center_image,
                    "rhs_lipschitz": sp.Rational(567, 1000),
                    "contraction": CONTRACTION_CAP,
                    "self_map": rounded_self_map,
                    "ball_margin": RADIUS-rounded_self_map},
        "strict_margins": margins,
    }
