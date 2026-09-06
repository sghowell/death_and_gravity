"""Fraction-only independent replay of all complete-map extension constants."""

from fractions import Fraction as Q


def serialize(value):
    if isinstance(value, dict):
        return {name: serialize(item) for name, item in value.items()}
    return str(value)


def replay(a11, a13):
    old = a11["derived_constants"]["full_map"]
    inv = a13["derived_constants"]["rational_gate"]
    delta, source_window = Q(old["delta"]), Q(old["length"])
    if delta != Q(inv["delta"]):
        raise ValueError("the physical amplitude changed between pinned inputs")
    length, radius, weight = Q(1, 10**6), Q(2, 10**7), Q(9)
    sigma, inverse = Q(inv["weighted_sigma"]), Q(inv["weighted_norm_upper"])
    alpha, af = 1/sigma, Q(inv["scale"])
    b, pmax = Q(3, 10**12), Q(1, 10**5)
    c, v = Q(old["source_cap"]), Q(old["geometry"]["background_uprime_cap"])
    m, history = Q(old["geometry"]["history_uprime_cap"]), Q(old["geometry"]["history"])
    qb, qmax = Q(old["background_q_and_qprime_cap"]), Q(old["auxiliary_q_cap"])
    ea, df = delta+9*length, delta+length/4
    metric = {"u": alpha, "h": alpha**2, "log_a": alpha**3,
              "a_squared": 18*alpha**3, "inverse_a_squared": alpha**3/2,
              "d": alpha**3/2, "dprime": alpha**2/2}
    curvature = alpha*((b/2+Q(1, 120))*metric["u"]+(b/30+Q(1, 60))*metric["h"])
    source = 8*c*(9*metric["h"]+alpha*(9*metric["u"]+54*b*metric["h"]))
    e_full = alpha*(9*metric["u"]+b*metric["a_squared"])/(60*delta)
    p = e_full+curvature+source
    q = alpha*(p/4+pmax*metric["inverse_a_squared"])
    auxiliary = 9*q+18*qmax*metric["h"]+qmax*metric["a_squared"]
    e_rem = alpha*(ea*metric["u"]+b*metric["a_squared"])/(60*delta)
    local = df+(v+weight*radius)*metric["d"]+metric["u"]/4+b*metric["dprime"]
    z = m*history**3
    response = {"strength": z, "log_ratio_upper": Q(20), "exponential_upper": Q(2),
                "quadratic": m*history*length**2*Q(45, 4),
                "higher": 36*m**2*history**5*length}
    response["total"] = response["quadratic"]+response["higher"]
    full = e_rem+curvature+source+auxiliary+local+response["total"]
    center_p = 8*c*(3+length*(1+9*b))
    center_rhs = center_p*(1+Q(9, 4)*length)
    p_derived, q_derived = 9*qb+center_p+p*weight*radius, qb+length*pmax/4
    contraction, center_image = inverse*full, inverse*center_rhs
    self_map, distance = center_image+contraction*radius, center_image/(1-contraction)
    smooth = inverse*(df+z+36*z**2)
    a0_fourth = af**4-delta
    h0_fourth = af**8/a0_fourth**3
    margins = {
        "longer_than_original": length-source_window,
        "future_potential": b-delta*Q(old["A8_inputs"]["U0_per_delta"])-weight*radius*alpha,
        "whole_history_derivative": m-v-weight*radius,
        "background_derivative": v-delta*Q(old["A8_inputs"]["U1_per_delta"]),
        "old_plateau_coverage": 3-af-length,
        "a0_lower": a0_fourth-Q(12, 5)**4, "a0_upper": Q(13, 5)**4-a0_fourth,
        "h0_lower": h0_fourth-Q(3, 8)**4, "h0_upper": Q(5, 12)**4-h0_fourth,
        "h_lower_first_exit": Q(3, 8)-length/3-Q(1, 3),
        "h_upper_first_exit": Q(1, 2)-Q(5, 12)-length/3,
        "h_derivative": Q(1, 3)-Q(1, 4)-b,
        "a_upper_first_exit": 3-Q(13, 5)/(1-length/2),
        "initial_d_difference": 8*a0_fourth-1,
        "pointwise_exponential": weight-Q(87, 32)**2,
        "P_barrier": pmax-p_derived, "q_barrier": qmax-q_derived,
        "rounded_rhs": Q(9, 10**6)-full,
        "rounded_contraction": Q(3, 10**6)-contraction,
        "rounded_center_image": Q(9, 10**8)-center_image,
        "ball": radius-self_map, "weighted_fixed_point": Q(9, 10**8)-distance,
        "old_uniqueness_ball": Q(old["radius"])-weight*distance,
        "highest_jet_contraction": Q(1, 10**4)-smooth,
    }
    if sigma*length != 2 or Q(8, 3)**20 < history/length or z > Q(1, 4):
        raise ValueError("the weighted response/exponential domain changed")
    if any(value <= 0 for value in margins.values()):
        raise ValueError("an independent weighted-extension margin failed")
    return serialize({
        "delta": delta, "source_window": source_window, "source_flat_end": source_window/4,
        "source_off_from": source_window/2, "length": length,
        "extension_factor": length/source_window, "source_free_span": length-source_window/2,
        "sigma": sigma, "weight_cap": weight, "weighted_radius": radius,
        "inverse_cap": inverse, "scale_freeze": af,
        "geometry": {"a_min": Q(2), "a_max": Q(3), "h_min": Q(1, 3), "h_max": Q(1, 2),
                     "u_cap": b, "background_uprime_cap": v, "history_uprime_cap": m,
                     "history": history, "pointwise_ball_X": weight*radius,
                     "pointwise_ball_W": weight*radius*alpha,
                     "rolling_a_squared_cap": ea, "local_d_difference_cap": df},
        "metric_gains": metric, "source_cap": c,
        "auxiliary_P_cap": pmax, "auxiliary_q_cap": qmax,
        "auxiliary_P_derived": p_derived, "auxiliary_q_derived": q_derived,
        "pairs": {"Einstein_full_P": e_full, "curvature": curvature, "source": source,
                  "P": p, "q": q, "auxiliary_product": auxiliary,
                  "Einstein_remainder": e_rem, "local_Wick": local,
                  "actual_mode_response": response, "full_G": full},
        "center": {"P_difference": center_p, "rhs": center_rhs},
        "gate": {"contraction": contraction, "center_image": center_image, "self_map": self_map,
                 "fixed_point_weighted": distance, "fixed_point_pointwise": weight*distance},
        "regularity": {"strength": z, "highest_jet_contraction": smooth},
        "strict_margins": margins,
    })
