"""Fraction-only reconstruction from the pinned A.11 input report.

No symbolic engine, primary bound helper, frequency quadrature or rounded
floating-point comparison is used in this independent arithmetic path.
"""

from fractions import Fraction as Q


def serialize(value):
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    return str(value)


def replay(a11):
    full = a11["derived_constants"]["full_map"]
    geometry = full["geometry"]
    t, m = Q(geometry["history"]), Q(geometry["history_uprime_cap"])
    length, delta = Q(full["length"]), Q(full["delta"])
    h0, h1, h2 = Q(geometry["h_max"]), Q(1, 3), Q(1, 3)
    b, ell = m*t, length/2
    q0, q1 = 2*b, 2*m+b*(2*b*t)
    aa0, aa1, aa2 = t*q0, q0, q1
    c0 = (1+h0)*aa0/2
    c1 = (1+h0)*aa1/2+h1*aa0/2
    c2 = (1+h0)*aa2/2+h1*aa1+h2*aa0/2
    back = (1+h0)/4
    p0, p1 = back*q0, back*q1+h1*q0/4
    l0 = 2*t*b
    l1 = (q0+t*q1)/2
    r1 = m*l0+b*l1
    ir = (1+h0)*b*t*t+2*b*t
    # Derive each transformed sampler bound from proper H,Hdot,a^-1.
    hp, hpdot, ainverse = Q(1, 4), Q(1, 7), Q(1, 2)
    f0 = ainverse**2*ell**2
    f1 = ainverse*(ell+Q(3, 2)*hp*ell**2)
    f2 = 1+2*hp*ell+(Q(3, 4)*hp**2+Q(3, 2)*hpdot)*ell**2
    weighted = hp*ell+(hpdot+hp**2/2)*ell**2
    flat = f2+Q(3, 2)*weighted
    forward_root = 2*Q(1, 10**5)*(c0*f2+2*c1*f1+c2*f0)
    local_root = 2*Q(1, 10**5)*(p0*f1+p1*f0)
    history_root = 4*back*m*2*f0
    remainder_root = 2*back*t*r1*f0
    infrared_root = 2*ir*f0
    root = flat+forward_root+local_root+history_root+remainder_root+infrared_root
    distance = Q(full["actual_gate"]["fixed_point_distance"])
    rstar = Q(72, 10**9)
    local_b = Q(geometry["u_cap"])
    loss = (9*length+Q(9, 2)*length**2+(3*local_b+Q(3, 4))*length**3)*rstar
    lower_n = delta/Q(54)
    eed = lower_n/(960*delta*3**4)
    span = 3*ell
    index = (3-span**2/8)/span
    q2 = Q(2880, 16)*2*delta
    margins = {
        "new_root": 1+Q(1, 10**8)-root,
        "new_coefficient": 2-root*root,
        "distance": rstar-distance,
        "new_EED_numerator": delta/54-loss,
        "index_test": index-Q(3, 4),
    }
    if any(value <= 0 for value in margins.values()):
        raise ValueError("an independent SEE/QSEI or nonfocusing margin failed")
    return {"constants": serialize({
        "history": t, "derivative_cap": m, "potential_cap": b,
        "length": length, "delta": delta, "conformal_span": ell,
        "mode_exponent": b*t*t/2,
        "q0": q0, "q1": q1, "A0": aa0, "A1": aa1, "A2": aa2,
        "forward0": c0, "forward1": c1, "forward2": c2,
        "local0": p0, "local1": p1, "backward_coefficient": back,
        "bminus1_inverse_k": l0, "bprime_inverse_k": l1,
        "qprime_remainder": r1, "infrared_error": ir,
        "sampler0": f0, "sampler1": f1, "sampler2": f2,
        "weighted_sampler": weighted, "auxiliary_root": flat,
        "forward_root": forward_root, "local_root": local_root,
        "history_root": history_root, "remainder_root": remainder_root,
        "infrared_root": infrared_root, "total_root": root,
        "derived_coefficient": root*root,
        "actual_distance": distance, "distance_cap": rstar,
        "numerator_loss": loss, "positive_EED": eed,
        "proper_span": span, "index_lower": index,
        "normal_jacobian": Q(2, 3)**3,
        "Q2_over_T0_squared": q2, "Q2_over_duration_squared": q2/span**2,
    }), "strict_margins": serialize(margins)}
