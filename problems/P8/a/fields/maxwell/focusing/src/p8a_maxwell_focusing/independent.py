"""Fraction-only sampler/envelope replay using the pinned photon EED polynomial."""

from fractions import Fraction as F


def derivative(poly):
    return tuple(F(i)*value for i, value in enumerate(poly) if i)


def multiply(left, right):
    result = [F(0)]*(len(left)+len(right)-1)
    for i, value in enumerate(left):
        for j, other in enumerate(right):
            result[i+j] += value*other
    return tuple(result)


def integral(poly):
    return sum((value/F(i+1) for i, value in enumerate(poly)), F(0))


def moments():
    p = tuple(map(F, (0, 0, 3, -2)))
    first, second = derivative(p), derivative(derivative(p))
    return {"zeroth": integral(multiply(p, p)), "first": integral(multiply(first, first)),
            "second": integral(multiply(second, second)),
            "pole_first": integral(multiply(first[1:], first[1:])),
            "pole_zero": integral(multiply(p[2:], p[2:])),
            "cross": integral(multiply(p, first))}


def loss(parent, caps, beta_abs):
    coefficients = parent["derived_constants"]["stress_jet_coefficients"]["EED"]
    result = F(0)
    for piece, factor in (("universal", F(1)), ("beta_M", beta_abs)):
        for powers, value in coefficients[piece].items():
            term = factor*abs(F(value))
            for cap, power in zip(caps, map(int, powers.split(",")), strict=True):
                term *= cap**power
            result += term
    return result


def countercontrol(parent):
    amplitude, width, delta = F(31, 20), F(1, 1000), F(1, 10**8)
    max_first, max_second, max_third = F(30, 16), F(6), F(60)
    jets = [amplitude, amplitude*max_first, amplitude*max_second, amplitude*max_third]
    source_eed = parent["derived_constants"]["stress_jet_coefficients"]["EED"]
    constant_reference = F(source_eed["universal"]["4,0,0,0"])/360
    if source_eed["beta_M"].get("4,0,0,0", "0") != "0":
        raise ValueError("the pinned finite-beta tensor no longer vanishes at constant H")
    required = 3*amplitude**2+delta*constant_reference*amplitude**4
    magnitude = amplitude*(1-10*width**3)
    return {"amplitude": amplitude, "mollifier_width_over_tau": width,
            "global_Hubble_jet_caps": jets,
            "cap_margins": [cap-value for cap, value in zip(map(F, (2, 4, 16, 96)), jets, strict=True)],
            "past_contraction_magnitude_lower": magnitude, "history_margin": magnitude-F(3, 2),
            "second_switch_derivative_square_margin": max_second**2-225*F(4, 27),
            "future_static_by_time_over_tau": 1+width,
            "future_log_scale_factor_loss_upper": amplitude*(1+width),
            "reference_vacuum_required_sigma_at_Lambda_zero": required,
            "required_sigma_exceeds_one_by": required-1,
            "future_timelike_and_null_complete_geometry": True,
            "initial_pointwise_SEC_holds": False, "actual_allowed_Maxwell_SEE_solution": False,
            "geometric_assumptions_alone_imply_the_conclusion": False}


def serialize(value):
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(item) for item in value]
    if isinstance(value, bool):
        return value
    return str(value)


def replay(parent):
    r, h, delta = F(1, 100), F(3, 2), F(1, 10**8)
    caps = tuple(map(F, (2, 4, 16, 96)))
    m = moments()
    v = loss(parent, caps, F(1))
    q = F(3, 4)*caps[0]**2+F(3, 2)*caps[1]
    # Norms are obtained by exact polynomial integration; the replacements
    # below are explicitly checked upper square-root bounds.
    radical_margins = {"sqrt12": F(7, 2)**2-m["second"],
                       "sqrt_6_over_5": F(11, 10)**2-m["first"],
                       "sqrt_13_over_35": F(2, 3)**2-m["zeroth"],
                       "sqrt_13_over_3": F(25, 12)**2-m["pole_zero"]}
    if not all(value > 0 for value in radical_margins.values()):
        raise ValueError("an independent rational square-root majorant failed")
    past_root = F(7, 2)+2*caps[0]*F(11, 10)*r+q*F(2, 3)*r*r
    future_root = F(7, 2)*(1+2*caps[0])+q*F(25, 12)
    past_cost = past_root**2/r**3+m["zeroth"]*r*v/360
    future_cost = future_root**2+m["pole_zero"]*v/360
    total = past_cost+future_cost
    gain = 6*h*m["cross"]+3*h*h*r*m["zeroth"]
    gradient = 3*m["first"]
    quantum, weight = delta*total, (1+r)*m["zeroth"]
    base, sourced = gain-gradient-quantum, gain-gradient-quantum-weight
    constants = {"ratio": r, "contraction": h, "delta": delta, "caps": list(caps),
                 "beta_absolute_cap": F(1), "past_reference_loss": v, "future_reference_loss": v,
                 "past_root": past_root, "future_root": future_root,
                 "past_cost": past_cost, "future_cost": future_cost, "total_cost": total,
                 "history_gain": gain, "future_gradient": gradient, "quantum_cost": quantum,
                 "source_weight": weight, "zero_source_margin": base, "sigma_one_margin": sourced,
                 "cost_coarsening_gap": F(13000000)-total,
                 "zero_source_margin_above_three_quarters": base-F(3, 4),
                 "sourced_margin_above_one_third": sourced-F(1, 3)}
    if min(constants[key] for key in ("cost_coarsening_gap", "zero_source_margin_above_three_quarters",
                                      "sourced_margin_above_one_third")) <= 0:
        raise ValueError("the independent macroscopic threshold lacks the claimed strict margin")
    return serialize({"constants": constants, "cubic_moments": m,
                      "radical_margins": radical_margins, "geometric_countercontrol": countercontrol(parent),
                      "engine": "stdlib Fraction polynomial integration and pinned Maxwell EED coefficients",
                      "imports_new_primary_formula_modules": False})
