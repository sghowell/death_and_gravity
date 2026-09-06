"""Fraction-only replay from the pinned A.7 certificate, without SymPy.

The independent Laurent-polynomial checks are algebra controls, not physical
flat-past states. Their nonzero curvature and histories test reconstruction
separately from the positive majorant arithmetic.
"""

from fractions import Fraction as F

NAMES = ("density", "pressure", "EED")


def add(*values):
    result = {}
    for value in values:
        for exponent, coefficient in value.items():
            result[exponent] = result.get(exponent, F(0))+coefficient
    return {n: c for n, c in result.items() if c}


def scale(value, coefficient):
    return {n: c*coefficient for n, c in value.items() if c*coefficient}


def multiply(left, right):
    return add(*({i+j: a*b} for i, a in left.items() for j, b in right.items()))


def derivative(value):
    return {n-1: n*c for n, c in value.items() if n*c}


def powerlaw_control(power, s_power):
    """Exact a=eta^power, S=eta^s_power; an algebra-only control."""
    if power < 2 or s_power <= 2:
        raise ValueError("The nonzero-curvature Laurent controls need power>=2 and S power>2")
    h = {-1: F(power)}
    u = {-2: F(-power*(power-1))}
    s = {s_power: F(1)}
    j = {s_power-2: F(2*power*(power-1), s_power-2)}
    ell = {-4: F(-power**3*(power-1)**2, 4)}
    inverse_a2, inverse_a4 = {-2*power: F(1)}, {-4*power: F(1)}
    h2, u2 = multiply(h, h), multiply(u, u)
    h4 = multiply(h2, h2)
    rho_inner = add(scale(h4, F(1, 960)),
        scale(add(scale(multiply(h, derivative(s)), -1), multiply(add(h2, scale(u, -1)), s), j), F(1, 8)),
        scale(ell, F(-1, 32)))
    p_inner = add(scale(add(scale(h4, 5), scale(multiply(h2, u), 4)), F(1, 2880)),
        scale(add(derivative(derivative(s)), scale(multiply(h, derivative(s)), -3),
                  multiply(add(scale(h2, 3), u), s), j), F(1, 24)),
        scale(add(u2, scale(ell, -1)), F(1, 96)))
    rho, pressure = multiply(inverse_a4, rho_inner), multiply(inverse_a4, p_inner)
    box = lambda value: multiply(inverse_a2,
        add(derivative(derivative(value)), scale(multiply(h, derivative(value)), 2)))
    wick = scale(multiply(inverse_a2, add(s, scale(u, F(-1, 10)))), F(1, 4))
    curvature = {-2*power-2: F(6*power*(power-1))}
    proper_h = {-power-1: F(power)}
    proper_h_dot = multiply({-power: F(1)}, derivative(proper_h))
    proper_h_squared = multiply(proper_h, proper_h)
    v1 = add(scale(multiply(curvature, curvature), F(1, 288)),
             scale(multiply(proper_h_squared, add(proper_h_dot, proper_h_squared)), F(-1, 60)),
             scale(box(curvature), F(-1, 120)))
    trace = add(scale(box(wick), F(-1, 2)), scale(v1, F(-1, 4)))
    checks = {"trace": add(rho, scale(pressure, -3), scale(trace, -1)),
              "conservation": add(derivative(rho), scale(multiply(h, add(rho, pressure)), 3)),
              "J_history": add(derivative(j), scale(multiply(derivative(u), s), -1)),
              "L_history": add(derivative(ell), scale(multiply(h, u2), -1))}
    if any(checks.values()):
        raise ValueError("Independent Laurent trace/conservation identity failed")
    return {"a_power": power, "S_power": s_power, "all_exact_Laurent_residuals_zero": True,
            "rho_at_eta_one_in_hbar_over_pi2_units": str(sum(rho.values(), F(0))),
            "pressure_at_eta_one_in_hbar_over_pi2_units": str(sum(pressure.values(), F(0))),
            "control_is_a_physical_flat_past_state": False}


def replay(prior_report):
    if prior_report.get("claim") != "P8-A.7":
        raise ValueError("Need the pinned A.7 finite-amplitude input report")
    geometry = prior_report["derived_finite_geometry"]
    old = prior_report["independently_replayed_uniform_bounds"]
    b = list(map(F, geometry["U_derivative_bounds_through_5_per_delta"]))
    cap = F(geometry["delta_bar"])
    k = [3*b[1]+3*b[0], 3*b[2]+3*b[1]+2*b[0],
         3*b[3]+3*b[2]+4*b[1]+(4+cap*b[0])*b[0]]
    s0, s1, s2 = [k[j]/2+b[j]/10 for j in range(3)]
    h, active, inverse_a4 = F(2), F(3), F(1, 8)
    derivative_parts = {
        "density": inverse_a4*(h*s1+(h*h+cap*b[0]+active*cap*b[1])*s0)/8,
        "pressure": inverse_a4*(s2+3*h*s1+(3*h*h+cap*b[0]+active*cap*b[1])*s0)/24,
        "EED": inverse_a4*(s2+4*h*s1+(4*h*h+2*active*cap*b[1])*s0)/16}
    reference_denominators = dict(zip(NAMES, (960, 576, 320), strict=True))
    euler = {name: F(slope, reference_denominators[name]*2**12)
             for name, slope in zip(NAMES, (8, 12, 12), strict=True)}
    first = {name: derivative_parts[name]+euler[name] for name in NAMES}
    old_second = {name: F(old["stress_error_constants"][name]) for name in NAMES}
    history_norm = active*h*b[0]**2
    r_squared = {"density": inverse_a4*history_norm/32,
                 "pressure": inverse_a4*(b[0]**2+history_norm)/96,
                 "EED": inverse_a4*(b[0]**2+2*history_norm)/64}
    second = {name: old_second[name]+r_squared[name] for name in NAMES}
    first_rounded = dict(zip(NAMES, (14000, 236000, 360000), strict=True))
    second_rounded = dict(zip(NAMES, (8*10**14, 9*10**14, 17*10**14), strict=True))
    first_margins = {name: first_rounded[name]-first[name] for name in NAMES}
    second_margins = {name: second_rounded[name]-second[name] for name in NAMES}
    if any(value <= 0 for margins in (first_margins, second_margins) for value in margins.values()):
        raise ValueError("Independent rounded stress coefficient comparison failed")
    amplitude = F(1, 10**14)
    if not 0 < amplitude <= cap <= F(1, 2):
        raise ValueError("The finite physical-amplitude example is outside its derived domain")
    error = {name: amplitude*first_rounded[name]+amplitude**2*second_rounded[name] for name in NAMES}
    reference_ratios = {name: reference_denominators[name]*3**8*error[name] for name in NAMES}
    c = dict(zip(NAMES, (3, 5, 9), strict=True))
    residual_brackets = {name: F(c[name], 1922)+2880*(first_rounded[name]+amplitude*second_rounded[name]) for name in NAMES}
    classical_multipliers = {"density": 27, "pressure": 81, "EED": 27}
    classical_ratios = {name: classical_multipliers[name]*amplitude**2*residual_brackets[name] for name in NAMES}
    if any(value >= F(1, 100) for value in reference_ratios.values()):
        raise ValueError("Independent one-percent actual/reference control failed")
    if any(value >= F(1, 10**17) for value in classical_ratios.values()):
        raise ValueError("Independent finite SEE residual/reference-classical-scale control failed")
    if F(2048, 961)/2**12 != F(1, 1922):
        raise ValueError("The literal 1922 frozen-defect denominator is wrong")
    return {
        "kernel_norm_constants": list(map(str, k)),
        "effective_S_Born_norm_constants": list(map(str, (s0, s1, s2))),
        "Euler_first_constants": {name: str(value) for name, value in euler.items()},
        "derivative_first_constants": {name: str(value) for name, value in derivative_parts.items()},
        "R_squared_history_second_constants": {name: str(value) for name, value in r_squared.items()},
        "first_constants": {name: str(value) for name, value in first.items()},
        "second_constants": {name: str(value) for name, value in second.items()},
        "strict_first_rounding_margins": {name: str(value) for name, value in first_margins.items()},
        "strict_second_rounding_margins": {name: str(value) for name, value in second_margins.items()},
        "rounded_first_constants": first_rounded, "rounded_second_constants": second_rounded,
        "literal_frozen_denominator": 1922,
        "finite_physical_amplitude_example": {
            "delta_upper": str(amplitude), "delta_bar": str(cap),
            "reference_stress_error_ratios": {name: str(value) for name, value in reference_ratios.items()},
            "all_reference_ratios_below_one_percent": True,
            "actual_reference_state_density_pressure_EED_positive_on_target": True,
            "positivity_claim_for_arbitrary_Hadamard_states": False,
            "actual_SEE_residual_bracket_constants": {name: str(value) for name, value in residual_brackets.items()},
            "actual_SEE_residual_over_reference_classical_scales": {name: str(value) for name, value in classical_ratios.items()},
            "all_classical_scale_ratios_below_1e_minus_17": True,
            "normalization_denominators": "positive A.3 quantum components, or positive zeroth-order classical Einstein/radiation components; never unknown actual values",
            "physical_epsilon_one_regime": "t_star=A*eta_star^2/2 >= 2*10^7*sqrt(d), with d=kappa*hbar/(46080*pi^2)",
            "the_metric_is_an_exact_or_nearby_SEE_solution": False},
        "independent_Laurent_controls": [powerlaw_control(p, m) for p, m in ((2, 5), (2, 7), (3, 5), (4, 6))],
    }
