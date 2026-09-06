"""Fraction-only reconstruction, independent of SymPy and the primary modules.

Polynomials use four proper Hubble jets. Stress is reconstructed from
curvature contractions and the variational tensor, not its expanded API.
"""

from fractions import Fraction as F

ZERO = (0, 0, 0, 0)


def constant(value):
    value = F(value)
    return {ZERO: value} if value else {}


def jet(index):
    powers = [0]*4
    powers[index] = 1
    return {tuple(powers): F(1)}


def add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for powers, value in polynomial.items():
            result[powers] = result.get(powers, F(0))+value
    return {powers: value for powers, value in result.items() if value}


def scale(polynomial, factor):
    return {powers: value*F(factor) for powers, value in polynomial.items()
            if value*F(factor)}


def multiply(left, right):
    result = {}
    for powers, value in left.items():
        for other, coefficient in right.items():
            key = tuple(a+b for a, b in zip(powers, other, strict=True))
            result[key] = result.get(key, F(0))+value*coefficient
    return {powers: value for powers, value in result.items() if value}


def derivative(polynomial):
    """Proper-time jet derivative; never silently discard a missing jet."""
    result = {}
    for powers, value in polynomial.items():
        if powers[-1]:
            raise ValueError("fourth proper Hubble jet lies outside this finite audit")
        for index, exponent in enumerate(powers[:-1]):
            if exponent:
                target = list(powers)
                target[index] -= 1
                target[index+1] += 1
                target = tuple(target)
                result[target] = result.get(target, F(0))+value*exponent
    return {powers: value for powers, value in result.items() if value}


def evaluate(polynomial, values):
    result = F(0)
    for powers, coefficient in polynomial.items():
        term = coefficient
        for value, power in zip(values, powers, strict=True):
            term *= F(value)**power
        result += term
    return result


def coefficients(polynomial):
    return {",".join(map(str, powers)): str(value)
            for powers, value in sorted(polynomial.items(), reverse=True)}


def covariant_polynomials():
    h, hd = jet(0), jet(1)
    h2 = multiply(h, h)
    r = scale(add(hd, scale(h2, 2)), 6)
    r00 = scale(add(hd, h2), 3)
    rii = scale(add(hd, scale(h2, 3)), -1)
    r_squared = multiply(r, r)
    ricci_squared = add(multiply(r00, r00), scale(multiply(rii, rii), 3))
    h3_rho = add(multiply(r00, r00), scale(multiply(r, r00), F(-2, 3)),
                 scale(ricci_squared, F(-1, 2)), scale(r_squared, F(1, 4)))
    h3_trace = add(scale(ricci_squared, -1), scale(r_squared, F(1, 3)))
    i_rho = add(scale(multiply(r, r00), 2), scale(r_squared, F(-1, 2)),
                scale(multiply(h, derivative(r)), -6))
    i_trace = scale(add(derivative(derivative(r)),
                        scale(multiply(h, derivative(r)), 3)), -6)
    universal_rho, universal_trace = scale(h3_rho, 62), scale(h3_trace, 62)
    primary = {"rho": (universal_rho, i_rho), "trace": (universal_trace, i_trace),
               "pressure": (scale(add(universal_rho, scale(universal_trace, -1)), F(1, 3)),
                            scale(add(i_rho, scale(i_trace, -1)), F(1, 3))),
               "EED": (add(universal_rho, scale(universal_trace, F(-1, 2))),
                       add(i_rho, scale(i_trace, F(-1, 2)))),
               "I_rho": (i_rho, {}), "I_trace": (i_trace, {}),
               "I_EED": (add(i_rho, scale(i_trace, F(-1, 2))), {})}
    for component in (0, 1):
        rho, pressure = primary["rho"][component], primary["pressure"][component]
        if add(derivative(rho), scale(multiply(h, add(rho, pressure)), 3)):
            raise ValueError("independent curvature stress is not conserved")
    return primary


def proper_ibp_polynomials():
    h, hd = jet(0), jet(1)
    p = scale(h, -2)
    q = add(scale(multiply(h, h), F(3, 4)), scale(hd, F(-3, 2)))
    first = add(multiply(p, p), scale(derivative(p), -1), scale(q, -2))
    zero = add(multiply(q, q), derivative(derivative(q)),
               scale(derivative(multiply(p, q)), -1))
    return {"first": first, "zero": zero}


def radiation():
    # H^(j)(1)=(-1)^j j!/2; homogeneity supplies the indicated t powers.
    jets = (F(1, 2), F(-1, 2), F(1), F(-3))
    data, ibp = covariant_polynomials(), proper_ibp_polynomials()
    rho, eed = (evaluate(data[key][0], jets)/2880 for key in ("rho", "EED"))
    zero = evaluate(ibp["zero"], jets)
    return {"proper_time_positive": True, "operator_first_times_t": F(-1),
            "operator_zero_times_t_squared": F(3, 4)*jets[0]**2-F(3, 2)*jets[1],
            "squared_first_times_t_squared": evaluate(ibp["first"], jets),
            "squared_zero_times_t_fourth": zero,
            "reference_rho_pi_squared_t_fourth_over_hbar": rho,
            "reference_EED_pi_squared_t_fourth_over_hbar": eed,
            "absolute_zero_times_t_fourth": zero-8*eed,
            "finite_ambiguity_in_density": evaluate(data["rho"][1], jets),
            "finite_ambiguity_in_EED": evaluate(data["EED"][1], jets)}


def envelope_arithmetic(beta):
    jets, duration = [F(1)]*4, F(1)
    h, hd = jets[:2]
    first_norm = duration/3
    zeroth_norm = first_norm**2
    operator = F(1)+2*h*first_norm+(F(3, 4)*h*h+F(3, 2)*hd)*zeroth_norm
    eed = covariant_polynomials()["EED"]
    # At unit absolute caps, the coefficient l1 norm is the envelope.
    loss = sum(map(abs, eed[0].values()), F(0))
    loss += abs(beta)*sum(map(abs, eed[1].values()), F(0))
    return {"hubble_jet_caps": jets, "proper_duration": duration, "beta_M": F(beta),
            "sampler_factor": operator, "reference_loss_numerator": loss,
            "derivative_coefficient_pi_squared_over_hbar": operator**2/8,
            "zero_coefficient_pi_squared_over_hbar": loss/2880}


def source_controls():
    field = envelope_arithmetic(F(1))
    quantum = 2*field["zero_coefficient_pi_squared_over_hbar"]
    derivative_coefficient = 2*field["derivative_coefficient_pi_squared_over_hbar"]
    cases = {}
    for lam, lower in ((1, 2), (1, -2), (-1, 2), (-1, -2)):
        coarse, raw = F(max(lam, 0)+2*max(-lower, 0)), F(lam-2*lower)
        cases[f"Lambda={lam},other_lower={lower}"] = {
            "Q2_pi_squared": derivative_coefficient, "Q0_quantum_pi_squared": quantum,
            "coarsened_source": coarse, "raw_source": raw,
            "upward_coarsening_gap": coarse-raw}
    return cases


def serialize(value):
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    if isinstance(value, bool):
        return value
    return str(value)


def replay():
    stress = {name: {piece: coefficients(poly)
                     for piece, poly in zip(("universal", "beta_M"), pair, strict=True)}
              for name, pair in covariant_polynomials().items()}
    arithmetic = {"arithmetic_plus_beta": envelope_arithmetic(F(1)),
                  "arithmetic_minus_beta": envelope_arithmetic(F(-1)),
                  "named_zero_type_D_arithmetic": envelope_arithmetic(F(0)),
                  "physical_calibration": False, "beta_M_fixed_by_scalar_gamma": False}
    return serialize({"stress_jet_coefficients": stress,
                      "proper_IBP_coefficients": {name: coefficients(poly)
                                                  for name, poly in proper_ibp_polynomials().items()},
                      "radiation_control": radiation(), "envelope_controls": arithmetic,
                      "source_controls": source_controls(),
                      "focusing_incompatibility": {"K_trace_times_tau_cap": F(3),
                          "A1_zeta_zero_gradient_lower": F(3)/F(3, 4),
                          "strict_incompatibility_gap": F(1),
                          "small_Q2_implies_focusing": False},
                      "engine": "stdlib Fraction; curvature contractions and inverse-metric variation",
                      "imports_primary_formula_modules": False})
