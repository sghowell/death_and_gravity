"""Independent Fraction-only replay of the preparation norm arithmetic.

No imports from the symbolic map, bounds, regularity, or SymPy. The
pointwise polynomial majorants are integrated here independently before
comparison with the main engine. Inputs are the pinned A.7/A.8 reports.
"""

from fractions import Fraction as F


def integral(coefficients, end):
    """Integrate an exact polynomial majorant on [0,end]."""
    return sum((F(coefficient)*end**(degree+1)/(degree+1)
                for degree, coefficient in coefficients.items()), F(0))


def replay(a7_report, a8_report):
    old7 = a7_report["independently_replayed_uniform_bounds"]
    old8 = a8_report["independently_replayed_finite_bounds"]
    delta, ell, radius = F(1, 10**14), F(1, 10**10), F(1, 10**6)
    b, c, v, qcap, pcap = F(2, 10**12), F(12, 10**9), F(1, 10**10), F(2, 10**8), F(1, 10**6)
    qb, m, history = F(1, 10**8), F(1, 10**5), F(3)
    old_s = [delta*F(old8["effective_S_Born_norm_constants"][j])
             +delta**2*F(old7["integrated_mode_remainder_constants"][j]) for j in range(2)]
    source = 3**4*delta*(F(old8["rounded_first_constants"]["density"])
                            +delta*F(old8["rounded_second_constants"]["density"])
                            +F(3, 2880*1922))
    ep = 8*c*(3+integral({0: 1+9*b}, ell))
    eq = ell*ep/4
    ef = ep+9*eq
    einstein = integral({1: F(9)/(60*delta), 3: 3*b/(60*delta)}, ell)
    curvature = integral({1: b/2+F(1, 120), 2: b/60+F(1, 120)}, ell)
    forcing_endpoint = 8*c*F(9, 2)*ell**2
    forcing_bulk = 8*c*integral({1: 9, 2: 27*b}, ell)
    kp = einstein+curvature+forcing_endpoint+forcing_bulk
    kq = ell*kp/4+integral({3: pcap/12}, ell)
    kg = kp+9*kq+2*qcap*(F(9, 2)*ell**2+F(3, 2)*ell**3)
    local_terms = [F(17, 30)+radius*ell**3/12, v*ell**3/12,
                   ell/4+radius*ell**3/4, b*ell**2/4]
    local = sum(local_terms, F(0))
    z = m*history**3
    mode2 = m*history*ell**2*(F(5, 4)+F(32, 2))
    higher = 18*m**2*history**5*ell*2
    rhs = kg+local+mode2+higher
    inverse = 2*ell/(1-ell)+F(2, 9*10**5)+F(1, 5)
    contraction, center_image = inverse*rhs, inverse*ef
    actual_selfmap = center_image+contraction*radius
    rounded_selfmap = F(21, 100)*F(3, 10**7)+F(3, 25)*radius
    smooth = F(21, 100)*(F(17, 30)+z+18*z**2/(1-2*z))
    pderived = 9*qb+ep+kp*radius
    qderived = qb+ell*pcap/4
    jets = list(map(F, old7["U_jet_bounds_per_delta"]))
    y0 = F(5, 2)
    a04 = y0**4-delta
    h04 = y0**8/a04**3
    margins = {
        "source": c-source, "future_potential": b-delta*jets[0]-ell*radius,
        "background_derivative": v-delta*jets[1],
        "full_history_derivative": m-delta*jets[1]-radius,
        "background_q": qb-old_s[0]/4,
        "background_qprime": qb-(old_s[1]+old_s[0])/4,
        "auxiliary_P": pcap-pderived, "auxiliary_q": qcap-qderived,
        "initial_a_lower": a04-F(12, 5)**4,
        "initial_a_upper": F(13, 5)**4-a04,
        "initial_h_lower": h04-F(3, 8)**4,
        "initial_h_upper": F(5, 12)**4-h04,
        "inverse_cap": F(21, 100)-inverse,
        "rhs_cap": F(567, 1000)-rhs,
        "contraction_cap": F(3, 25)-F(21, 100)*rhs,
        "center_rhs_cap": F(3, 10**7)-ef,
        "rounded_ball": radius-rounded_selfmap,
        "smooth_block": F(3, 25)-smooth,
    }
    if any(value <= 0 for value in margins.values()):
        raise ValueError("independent preparation margin failed")
    values = {
        "S0": old_s[0], "S1": old_s[1], "source_derived": source,
        "center_P": ep, "center_q": eq, "center_rhs": ef,
        "Einstein_P": einstein, "curvature_P": curvature,
        "source_P": forcing_endpoint+forcing_bulk, "pair_P": kp, "pair_q": kq,
        "auxiliary_Sprime": kg, "local_Wick_terms": local,
        "response_quadratic": mode2, "response_higher": higher,
        "rhs_lipschitz": rhs, "inverse_bound": inverse,
        "actual_contraction": contraction, "actual_center_image": center_image,
        "actual_self_map": actual_selfmap,
        "actual_distance": center_image/(1-contraction),
        "auxiliary_P_derived": pderived, "auxiliary_q_derived": qderived,
        "rounded_self_map": rounded_selfmap,
        "highest_jet_contraction": smooth,
    }
    return {"constants": {name: str(value) for name, value in values.items()},
            "strict_margins": {name: str(value) for name, value in margins.items()},
            "majorant_integration": "Fraction polynomial coefficients integrated exactly; no binary floating point",
            "actual_state_inputs": "unchanged pinned A7 mode norms and A8 raw-scheme density defect",
            "infinite_frequency_or_smoothness_proved_by_finite_samples": False}
