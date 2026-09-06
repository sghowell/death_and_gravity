"""Fraction-only fresh replay from A14's pinned complete-map constants."""

from fractions import Fraction as Q


def serialize(value):
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    if isinstance(value, bool):
        return value
    return str(value)


def replay(report):
    old = report["derived_constants"]["complete_weighted_map"]
    geo = old["geometry"]
    ell, delta = Q(old["length"]), Q(old["delta"])
    t, m, b_local = Q(geo["history"]), Q(geo["history_uprime_cap"]), Q(geo["u_cap"])
    a_min, a_max, h0 = Q(geo["a_min"]), Q(geo["a_max"]), Q(geo["h_max"])
    weight, inverse_sigma = Q(old["weight_cap"]), 1/Q(old["sigma"])
    distance = Q(old["gate"]["fixed_point_weighted"])
    distance_cap = Q(9, 10**8)
    # The gains are independently composed from the weight, not read back
    # from the root reference implementation or its reported loss.
    pair = a_max**2*(inverse_sigma+2*h0*inverse_sigma**2)
    pair += (b_local+h0**2)*2*a_max**2*inverse_sigma**3
    loss = weight*distance_cap*pair
    negative = loss/(960*delta*a_min**4)
    ref_margins = {"sharp_weighted_distance": distance_cap-distance,
                   "signed_reference_coarsening": Q(1, 40)-negative}
    credit = {"actual_weighted_distance": distance, "distance_cap": distance_cap,
              "weight_cap": weight, "numerator_weighted_pair": pair,
              "numerator_loss_upper": loss, "baseline_numerator_is_positive": True,
              "positive_reference_credit_asserted": False,
              "negative_EED_magnitude_upper": negative,
              "rounded_negative_EED_magnitude": Q(1, 40),
              "reference_EED_lower_dimensionless": -Q(1, 40),
              "strict_margins": ref_margins}
    b, h1, h2 = m*t, Q(1, 3), Q(1, 3)
    q0, q1 = 2*b, 2*m+2*t*b*b
    a0, a1, a2 = 2*b*t, q0, q1
    c0 = (1+h0)*a0/2
    c1 = ((1+h0)*a1+h1*a0)/2
    c2 = ((1+h0)*a2+2*h1*a1+h2*a0)/2
    back = (1+h0)/4
    p0, p1 = back*q0, back*q1+h1*q0/4
    lm, lp = 2*t*b, b+t*m+t*t*b*b
    remainder = m*lm+b*lp
    infrared = (1+h0)*b*t*t+2*b*t
    mode_margins = {"local_hprime": h1-h0*h0-b_local,
                    "local_hsecond": h2-m-2*h0*(b_local+h0*h0),
                    "strict_modulus_exponent": Q(1, 2)-b*t*t/2}
    modes = {"history": t, "potential_cap": b, "derivative_cap": m,
             "hubble_caps": [h0, h1, h2], "modulus_cap": Q(2),
             "modulus_exponent": b*t*t/2,
             "q_and_qprime": [q0, q1], "A_and_two_derivatives": [a0, a1, a2],
             "forward_two_derivatives": [c0, c1, c2],
             "local_one_derivative": [p0, p1], "backward_coefficient": back,
             "bminus1_inverse_k": lm, "bprime_inverse_k": lp,
             "qprime_minus_uprime_inverse_k": remainder, "infrared_error": infrared,
             "strict_margins": mode_margins}
    f0, f1 = ell*ell/4, ell/2+Q(3, 16)*ell*ell
    f2 = 1+ell/2+Q(117, 448)*ell*ell
    b1 = ell/4+Q(39, 224)*ell*ell
    auxiliary = f2+Q(3, 2)*b1
    root_span = Q(1, 1000)
    roots = {"forward": 2*root_span*(c0*f2+2*c1*f1+c2*f0),
             "local": 2*root_span*(p0*f1+p1*f0),
             "real_history_Parseval": 8*back*m*f0,
             "history_remainder": 2*back*t*remainder*f0,
             "infrared": 2*infrared*f0}
    root = auxiliary+sum(roots.values())
    penalty = Q(2, 5)*ell**4
    coefficient = root*root+penalty
    off = Q(old["source_off_from"])
    sampler_margins = {"source_free_width": ell-off,
                       "proper_Hdot": Q(1, 7)-(b_local+Q(1, 2))/4,
                       "root_coarsening": 1+Q(2, 10**6)-root,
                       "absolute_coefficient": 2-coefficient}
    sampler = {"source_off_from": off, "future_end": ell,
               "source_free_conformal_width": Q(old["source_free_span"]),
               "conformal_envelope": ell, "proper_span_upper": 3*ell,
               "sqrt_conformal_envelope": root_span,
               "proper_H_cap": Q(1, 4), "proper_Hdot_cap": Q(1, 7),
               "sampler_norms": [f0, f1, f2], "weighted_first_norm": b1,
               "auxiliary_root": auxiliary, "spectral_error_roots": roots,
               "total_spectral_root": root, "difference_coefficient": root*root,
               "signed_reference_penalty": penalty,
               "derived_absolute_coefficient": coefficient, "rounded_absolute_coefficient": Q(2),
               "strict_margins": sampler_margins}
    span, ricci = 3*ell, 3*(b_local+h0*h0)/(a_min*a_min)
    lower, expansion, q2 = 3/span-span/8, 3*h0/a_min, 360*delta
    focus_margins = {"Ricci": Q(1, 4)-ricci,
                     "positive_index_energy_factor": 3-span*span/8,
                     "no_comoving_index_trigger": lower-expansion}
    focus = {"proper_span_upper": span, "Ricci_times_T0_squared_cap": Q(1, 4),
             "expansion_times_T0_cap": expansion, "index_times_T0_lower": lower,
             "normal_jacobian_lower": (a_min/a_max)**3,
             "Q2_over_T0_squared": q2,
             "Q2_over_available_duration_squared_lower": q2/(span*span),
             "strict_margins": focus_margins}
    for margins in (ref_margins, mode_margins, sampler_margins, focus_margins):
        if any(value <= 0 for value in margins.values()):
            raise ValueError("an independent extended QSEI/focusing margin failed")
    if root_span**2 != ell or Q(old["source_free_span"]) != ell-off or off != Q(old["source_window"])/2:
        raise ValueError("the original-source or new-sampler domain changed")
    return serialize({"fresh_C1_scattering": modes, "signed_actual_reference": credit,
                      "fresh_proper_H2_sampling": sampler, "actual_comoving_index": focus})
