"""Fraction-only all-p tube, separated-envelope and source reconstruction."""

from fractions import Fraction as F
from math import comb, factorial


def serialize(value):
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    if isinstance(value, bool):
        return value
    return str(value)


def derivative(poly):
    return tuple(i*value for i, value in enumerate(poly) if i)


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
    first = derivative(p)
    second = derivative(first)
    return {"zeroth": integral(multiply(p, p)), "first": integral(multiply(first, first)),
            "second": integral(multiply(second, second)),
            "pole_first": integral(multiply(first[1:], first[1:])),
            "pole_zero": integral(multiply(p[2:], p[2:])),
            "cross": integral(multiply(p, first))}


def reference_caps():
    p0 = F(1, 2)
    return [F(factorial(j))*2**(j+1)/p0**j for j in range(4)]


def tube(anchored):
    r = F(1, 100)
    caps = [F(21, 10), F(9), F(70), F(800)]
    errors = [F(1, 100), F(1, 2), F(3), F(16)]
    baseline = F(1, 2)/(F(1, 4)+r)
    actual = baseline-errors[0]
    gaps = [a-b-c for a, b, c in zip(caps, reference_caps(), errors, strict=True)]
    if min(*gaps, actual-F(19, 10)) <= 0:
        raise ValueError("the independent actual-past tube has no strict margin")
    return {"error_caps": errors, "reference_jet_caps": reference_caps(),
            "expanded_past_caps": caps, "strict_cap_margins": gaps,
            "baseline_contraction_lower": baseline, "actual_contraction_lower": actual,
            "strict_contraction_margin": actual-F(19, 10),
            "observed_H0_times_tau_range": [F(2), F(2)] if anchored else [2-errors[0], 2+errors[0]],
            "anchored_at_observer": anchored, "log_scale_factor_error_upper": r*errors[0],
            "scale_factor_comparison_is_normalized_at_observer": True,
            "future_relative_caps_verified_by_tube": False,
            "actual_quantum_SEE_solution_asserted": False}


def geometric_data():
    lower_first = F(2, 3)/(F(1, 3)+F(1, 100))**2
    ricci_upper = 3*(-lower_first+F(1, 2)+(2+F(1, 100))**2)
    if ricci_upper >= 0:
        raise ValueError("the named small tube's initial Ricci control changed")
    return {"power_range": [F(1, 2), F(2, 3)],
            "dimensionless_history_interval": [-F(1, 100), F(0)],
            "ratio": F(1, 100), "contraction": F(19, 10),
            "reference_Hstar_times_tau": F(2),
            "reference_age_over_tau_range": [F(1, 4), F(1, 3)],
            "unanchored_C3_tube": tube(False), "anchored_C3_slice": tube(True),
            "this_tube_Ricci_upper_times_tau_squared": ricci_upper,
            "this_tube_demonstrates_initial_SEC_violation": False,
            "power_law_geometries_are_quantum_SEE_solutions": False}


def loss(photon_report, caps, piece):
    coefficients = photon_report["derived_constants"]["stress_jet_coefficients"]["EED"][piece]
    result = F(0)
    for powers, value in coefficients.items():
        exponents = tuple(map(int, powers.split(",")))
        if len(exponents) != 4 or sum((j+1)*power for j, power in enumerate(exponents)) != 4:
            raise ValueError("the transitive photon EED no longer has proper fourth-order weight")
        term = abs(F(value))
        for cap, power in zip(caps, exponents, strict=True):
            term *= cap**power
        result += term
    return result


def affine(photon_report, past, future):
    r = F(1, 100)
    m = moments()
    v0p, vbp = (loss(photon_report, past, part) for part in ("universal", "beta_M"))
    v0f, vbf = (loss(photon_report, future, part) for part in ("universal", "beta_M"))
    qp = F(3, 4)*past[0]**2+F(3, 2)*past[1]
    qf = F(3, 4)*future[0]**2+F(3, 2)*future[1]
    past_root = F(7, 2)+2*past[0]*F(11, 10)*r+qp*F(2, 3)*r*r
    future_root = F(7, 2)*(1+2*future[0])+qf*F(25, 12)
    c0 = past_root**2/r**3+future_root**2+(m["zeroth"]*r*v0p+m["pole_zero"]*v0f)/360
    cb = (m["zeroth"]*r*vbp+m["pole_zero"]*vbf)/360
    return {"past_caps": list(past), "future_caps": list(future), "C0": c0, "Cbeta": cb,
            "C_at_abs_beta_one": c0+cb,
            "C0_upper_margin": F(13000000)-c0, "Cbeta_upper_margin": F(13000000)-cb,
            "past_V0": v0p, "past_Vbeta": vbp, "future_V0": v0f, "future_Vbeta": vbf}


def costs(photon_report):
    r, h, joint, sigma = F(1, 100), F(19, 10), F(1, 10**8), F(5)
    past = (F(21, 10), F(9), F(70), F(800))
    future = tuple(map(F, (4, 128, 16384, 1048576)))
    data = affine(photon_report, past, future)
    if min(data["C0_upper_margin"], data["Cbeta_upper_margin"]) <= 0:
        raise ValueError("independent finite-beta affine costs exceed the named joint gate")
    m = moments()
    gain = 6*h*m["cross"]+3*h*h*r*m["zeroth"]
    gradient = 3*m["first"]
    weight = (1+r)*m["zeroth"]
    quantum = F(13000000)*joint
    source = weight*sigma
    margin = gain-gradient-quantum-source
    if margin <= F(1, 8):
        raise ValueError("independent sourced margin is not strictly above one eighth")
    gate = {"weighted_delta_upper": joint, "sigma_upper": sigma,
            "history_gain": gain, "future_gradient": gradient, "source_weight": weight,
            "quantum_cost_upper": quantum, "source_cost_upper": source,
            "strict_focusing_margin_lower": margin, "margin_above_one_eighth": margin-F(1, 8),
            "future_relative_caps_or_actual_SEE_verified": False}
    return {"expanded_affine_cost": data,
            "tight_reference_regression": affine(photon_report, reference_caps(), reference_caps()),
            "tight_regression_is_not_used_as_future_geometry": True,
            "universal_gate": gate, "uniform_bound_on_abs_beta_required": False,
            "generic_finite_beta_remains_explicit": True, "COST_UPPER": F(13000000)}


def proper_scales(anchored):
    data = tube(anchored)
    return {"tau": F(1), "actual_short_history_duration": F(1, 100),
            "reference_Hstar": F(2),
            "reference_power_law_age_range": [F(1, 4), F(1, 3)],
            "actual_observed_H0_range": data["observed_H0_times_tau_range"],
            "anchored_at_observer": anchored,
            "reference_power_law_age_is_actual_quantum_age": False}


def physical_data():
    cases = {}
    for lam, lower in ((1, 2), (1, -2), (-1, 2), (-1, -2)):
        penalty = F(max(lam, 0)+2*max(-lower, 0))
        cases[f"Lambda={lam},other_lower={lower}"] = {
            "exact_delta": F(0), "rational_delta_upper": F(0),
            "sigma": penalty, "separate_source_penalty": penalty,
            "Lambda": F(lam), "other_EED_lower": F(lower), "beta_M": F(-7),
            "exact_weighted_delta": F(0), "rational_weighted_delta_upper": F(0),
            "sufficient_quantum_budget": True, "sufficient_source_budget": penalty <= 5,
            "actual_field_state_source_or_metric_verified": False}
    return {"anchored_scales_tau_one": proper_scales(True),
            "unanchored_scales_tau_one": proper_scales(False),
            "Lambda_positive_fraction_of_reference_3Hstar_squared_upper": F(5)/(3*F(2)**2),
            "Lambda_positive_fraction_of_actual_3H0_squared_upper_unanchored": F(5)/(3*F(199, 100)**2),
            "signed_source_controls": cases,
            "ordinary_matter_nonnegative_EED_is_separate_assumption": True,
            "pure_Maxwell_specialization_has_other_EED_zero": True,
            "observed_cosmological_parameters_used": False}


def complete_comparator():
    poly = tuple(map(F, (0, 0, 0, 0, 35, -84, 70, -20)))
    fourth = poly
    for _ in range(4):
        fourth = derivative(fourth)
    unit_caps = [F(1), F(140, 64), F(420, 16), F(840, 4), sum(map(abs, fourth), F(0))]
    switch_caps = [value*8**j for j, value in enumerate(unit_caps)]
    base = [F(factorial(j), 2)/F(1, 8)**(j+1) for j in range(5)]
    raw = [sum((F(comb(j, k))*base[j-k]*switch_caps[k] for k in range(j+1)), F(0)) for j in range(5)]
    width = F(1, 10**12)
    errors = [width*value for value in raw[1:]]
    tube_gaps = [a-b for a, b in zip(tube(False)["error_caps"], errors, strict=True)]
    future_gaps = [a-b for a, b in zip(map(F, (4, 128, 16384, 1048576)), raw[:4], strict=True)]
    if min(tube_gaps) <= 0 or min(future_gaps) < 0:
        raise ValueError("independent complete metric does not satisfy the chosen geometric inputs")
    return {"power_range": [F(1, 2), F(2, 3)],
            "cutoff_ends_at_dimensionless_time": F(1, 8),
            "mollifier_width_over_tau": width, "raw_global_jet_caps_through_four": raw,
            "past_C3_error_upper": errors, "strict_past_tube_margins": tube_gaps,
            "future_cap_margins_from_global_bounds": future_gaps,
            "future_static_by_time_over_tau": F(1, 8)+width,
            "future_log_scale_factor_loss_upper": 4*(F(1, 8)+width),
            "covers_every_reference_power_in_interval": True,
            "future_timelike_and_null_complete": True,
            "actual_quantum_SEE_solution_asserted": False,
            "geometric_assumptions_alone_force_incompleteness": False}


def excluded_tight():
    q = F(1, 20)
    initial = [-F(199, 100), -F(11, 2), -F(33)]
    third = F(800)/(1-q)**4
    upper = initial[0]+initial[1]*q+initial[2]*q*q/2+third*q**3/6
    required = -F(21, 10)/(1-q)
    if required <= upper:
        raise ValueError("the rejected too-tight envelope no longer has the claimed direct contradiction")
    return {"test_time_over_tau": q, "initial_upper_jets": initial,
            "third_jet_upper_on_test_interval": third, "Taylor_Hubble_upper": upper,
            "tight_future_required_Hubble_lower": required, "strict_incompatibility_gap": required-upper,
            "QSEI_used_in_this_exclusion": False, "rejected_as_nontrivial_cosmological_calibration": True}


def replay(theorem_report, photon_report):
    cubic = moments()
    if serialize(cubic) != theorem_report["derived_constants"]["cubic_moments"]:
        raise ValueError("independently integrated cubic moments differ from the pinned theorem")
    return serialize({"geometry": geometric_data(), "calibration": costs(photon_report),
                      "dictionary": physical_data(),
                      "controls": {"all_p_complete_geometric_family": complete_comparator(),
                                   "rejected_past_caps_as_future_caps": excluded_tight()},
                      "cubic_moments": cubic,
                      "engine": "stdlib Fraction, polynomial integration, all-p extrema and transitive pinned Maxwell EED coefficients",
                      "imports_new_primary_formula_modules": False})
