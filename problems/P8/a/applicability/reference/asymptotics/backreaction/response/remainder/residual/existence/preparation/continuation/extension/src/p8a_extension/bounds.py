"""Complete weighted Banach gate with separately verified pointwise barriers."""

from functools import cache

import sympy as sp
from p8a_continuation import bounds as inverse
from p8a_preparation import bounds as original

from .weighted import metric_gains, prefix_response


@cache
def calibration():
    old, frozen = original.calibration(), inverse.calibration()
    delta, source_window = old["delta"], old["length"]
    if delta != frozen["delta"]:
        raise ValueError("the physical amplitude must remain unchanged")
    ell, radius, weight = sp.Rational(1, 10**6), sp.Rational(2, 10**7), sp.Integer(9)
    sigma, k = frozen["weighted_sigma"], frozen["weighted_norm_upper"]
    alpha = 1/sigma
    b, pmax = sp.Rational(3, 10**12), sp.Rational(1, 10**5)
    c, v = old["source_cap"], old["geometry"]["background_uprime_cap"]
    m, history = old["geometry"]["history_uprime_cap"], old["geometry"]["history"]
    qb, qmax = old["background_q_and_qprime_cap"], old["auxiliary_q_cap"]
    ea, df = delta+9*ell, delta+ell/4
    curvature = (b/2+sp.Rational(1, 120))*alpha**2+(b/30+sp.Rational(1, 60))*alpha**3
    source = 144*c*alpha**2+432*c*b*alpha**3
    einstein_full = (9*alpha**2+18*b*alpha**4)/(60*delta)
    p = einstein_full+curvature+source
    q = alpha*p/4+pmax*alpha**4/2
    auxiliary = 9*q+18*qmax*(alpha**2+alpha**3)
    einstein_remainder = (ea*alpha**2+18*b*alpha**4)/(60*delta)
    local = df+(v+weight*radius)*alpha**3/2+alpha/4+b*alpha**2/2
    response = prefix_response(m, history, ell)
    full = einstein_remainder+curvature+source+auxiliary+local+response["total"]
    center_p = 8*c*(3+ell*(1+9*b))
    center_rhs = center_p*(1+9*ell/4)
    p_derived = 9*qb+center_p+p*weight*radius
    q_derived = qb+ell*pmax/4
    contraction, center_image = k*full, k*center_rhs
    self_map = center_image+contraction*radius
    distance = center_image/(1-contraction)
    z = m*history**3
    smooth = k*(df+z+36*z**2)
    af = frozen["scale"]
    a0_fourth = af**4-delta
    h0_fourth = af**8/a0_fourth**3
    margins = {
        "longer_than_original": ell-source_window,
        "future_potential": b-delta*old["A8_inputs"]["U0_per_delta"]-weight*radius*alpha,
        "whole_history_derivative": m-v-weight*radius,
        "background_derivative": v-delta*old["A8_inputs"]["U1_per_delta"],
        "old_plateau_coverage": 3-af-ell,
        "a0_lower": a0_fourth-sp.Rational(12, 5)**4,
        "a0_upper": sp.Rational(13, 5)**4-a0_fourth,
        "h0_lower": h0_fourth-sp.Rational(3, 8)**4,
        "h0_upper": sp.Rational(5, 12)**4-h0_fourth,
        "h_lower_first_exit": sp.Rational(3, 8)-ell/3-sp.Rational(1, 3),
        "h_upper_first_exit": sp.Rational(1, 2)-sp.Rational(5, 12)-ell/3,
        "h_derivative": sp.Rational(1, 3)-sp.Rational(1, 4)-b,
        "a_upper_first_exit": 3-sp.Rational(13, 5)/(1-ell/2),
        "initial_d_difference": 8*a0_fourth-1,
        "pointwise_exponential": weight-sp.Rational(87, 32)**2,
        "P_barrier": pmax-p_derived, "q_barrier": qmax-q_derived,
        "rounded_rhs": sp.Rational(9, 10**6)-full,
        "rounded_contraction": sp.Rational(3, 10**6)-contraction,
        "rounded_center_image": sp.Rational(9, 10**8)-center_image,
        "ball": radius-self_map,
        "weighted_fixed_point": sp.Rational(9, 10**8)-distance,
        "old_uniqueness_ball": original.RADIUS-weight*distance,
        "highest_jet_contraction": sp.Rational(1, 10**4)-smooth,
    }
    if sigma*ell != 2 or any(value <= 0 for value in margins.values()):
        raise ValueError("a complete weighted-extension or pointwise barrier failed")
    return {
        "delta": delta, "source_window": source_window,
        "source_flat_end": source_window/4, "source_off_from": source_window/2,
        "length": ell, "extension_factor": ell/source_window,
        "source_free_span": ell-source_window/2,
        "sigma": sigma, "weight_cap": weight, "weighted_radius": radius,
        "inverse_cap": k, "scale_freeze": af,
        "geometry": {"a_min": sp.Integer(2), "a_max": sp.Integer(3),
                     "h_min": sp.Rational(1, 3), "h_max": sp.Rational(1, 2),
                     "u_cap": b, "background_uprime_cap": v,
                     "history_uprime_cap": m, "history": history,
                     "pointwise_ball_X": weight*radius,
                     "pointwise_ball_W": weight*radius*alpha,
                     "rolling_a_squared_cap": ea, "local_d_difference_cap": df},
        "metric_gains": metric_gains(sigma), "source_cap": c,
        "auxiliary_P_cap": pmax, "auxiliary_q_cap": qmax,
        "auxiliary_P_derived": p_derived, "auxiliary_q_derived": q_derived,
        "pairs": {"Einstein_full_P": einstein_full, "curvature": curvature,
                  "source": source, "P": p, "q": q, "auxiliary_product": auxiliary,
                  "Einstein_remainder": einstein_remainder, "local_Wick": local,
                  "actual_mode_response": response, "full_G": full},
        "center": {"P_difference": center_p, "rhs": center_rhs},
        "gate": {"contraction": contraction, "center_image": center_image,
                 "self_map": self_map, "fixed_point_weighted": distance,
                 "fixed_point_pointwise": weight*distance},
        "regularity": {"strength": z, "highest_jet_contraction": smooth},
        "strict_margins": margins,
    }
