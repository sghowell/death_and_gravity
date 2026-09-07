"""Separate Fraction polynomial/causal arithmetic; no SymPy or primary import."""

from fractions import Fraction
from math import comb

F = Fraction
R_MAX = F(1, 100)


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("The independent engine accepts integers and Fractions only")
    return Fraction(value)


def _add(*polynomials):
    result = [F(0)]*max(map(len, polynomials))
    for polynomial in polynomials:
        for index, value in enumerate(polynomial):
            result[index] += value
    return result


def _scale(polynomial, value):
    return [value*coefficient for coefficient in polynomial]


def _multiply(left, right):
    result = [F(0)]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] += a*b
    return result


def _power(polynomial, degree):
    result = [F(1)]
    for _ in range(degree):
        result = _multiply(result, polynomial)
    return result


def _at(polynomial, value):
    answer = F(0)
    for coefficient in reversed(polynomial):
        answer = answer*value+coefficient
    return answer


def positive_envelope_polynomials():
    """Build the complete increasing r^2-polynomial bounds coefficientwise."""
    one, t, cp, dp = [F(1)], [F(0), F(1)], [F(2), F(1, 100)], [F(1), F(1)]
    theta = [F(18, 5)]
    omega = _scale(_multiply([F(48), F(3, 25)], _power(dp, 5)), F(1, 10))
    theta_v = _scale(_add(_scale(_multiply(cp, _power(dp, 11)), F(48, 25)), [F(3, 5)]), 6)
    omega_v = _multiply(omega, _add([F(5)], _scale(_multiply(cp, _power(dp, 11)), F(6, 5))))
    norm = _add(theta, _scale(_multiply(t, theta_v), 2), _multiply(t, _power(theta, 2)))
    sw = _scale(_multiply([F(4), F(1, 100)], _power(dp, 6)), F(1, 10))
    cf_minus1 = _add(_scale(_multiply(_power(cp, 2), _power(dp, 8)), F(1, 4)), [F(-1)])
    if cf_minus1[0] != 0:
        raise ValueError("The E vanishing factor cannot be divided by t")
    cross_over_t = _scale(_multiply(sw, cf_minus1[1:]), 4)
    e_over_t = _add(cross_over_t, _scale(_multiply(omega, theta), 2))
    aa = _add(_scale(_add(_multiply(cp, _power(dp, 12)),
                              _scale(_multiply(_power(cp, 2), _power(dp, 8)), 2)), F(2, 5)), norm)
    cc = _add(_multiply(t, cross_over_t), _scale(omega, 2),
              _scale(_multiply(t, omega_v), 4), _scale(_multiply(t, _multiply(omega, theta)), 2))
    ba = _add(_scale(_add([F(8)], _scale(_multiply(_power(cp, 3), _power(dp, 20)), F(1, 4))), F(2, 5)),
              norm, _scale(_multiply(t, _power(omega, 2)), 4))
    return {"c_upper": cp, "d_upper": dp, "theta_over_u": theta, "omega_over_u": omega,
            "theta_v": theta_v, "omega_v": omega_v,
            "A": aa, "C": cc, "E_over_r_squared": e_over_t,
            "dc_and_fc": _scale(omega, 2), "b_analytic": ba,
            "ag_squared": _scale(_multiply(cp, _power(dp, 6)), F(2, 5)),
            "bg_squared": _scale(one, F(16, 5)),
            "jL_squared": _scale(_multiply(cp, _power(dp, 18)), F(1, 10)),
            "jH_squared": _scale(_power(dp, 6), F(4, 5))}


def coefficient_constants(r=R_MAX):
    r = rational(r)
    if not 0 < r <= F(1, 100):
        raise ValueError("Require 0<r<=1/100")
    result = {key: _at(poly, r*r) for key, poly in positive_envelope_polynomials().items()}
    result.update({"r": r, "t": r*r, "mass_N_delta_derivative": F(32),
                   "mass_N_v_derivative": F(1008), "mass_remainder": F(1211, 8),
                   "total_bounded_B": F(161)})
    return result


def response_constants(r=R_MAX):
    r = rational(r)
    if not 0 < r <= F(1, 100):
        raise ValueError("Require 0<r<=1/100")
    r2, r4 = r*r, r**4
    rr = F(101, 100)
    potential_integral = F(161)*2*r2*rr/3
    exponent = F(1, 2)+potential_integral
    feedback = F(5)*64*r4*(10+60*r2)
    qerror = F(292)*F(101, 100)+F(4, 3)*rr*r2*(10+60*r2)*323
    return {"r": r, "rho_over_r": rr, "heavy_perturbation_exponent": potential_integral,
            "full_heavy_exponent": exponent,
            "full_heavy_exp_upper": 1+exponent+exponent**2/(2*(1-exponent/3)),
            "GH_over_r2": F(4, 3)*rr, "u_GHprime_over_r2": F(14, 3)*rr,
            "Bstar_GH_over_r2": (12*F(4, 3)+10*F(14, 3))*rr,
            "GH_minus_G0_over_r4": F(8, 3)*F(322, 3)*rr**2,
            "GL_over_r2": F(4)/(1-18*r2), "feedback": feedback,
            "l_over_r2": F(5)*(1+64*r2)/(1-feedback),
            "l_error_over_r4": F(320)*(1+60*r2+360*r4),
            "F0_over_S": 1+50*r2+300*r4, "Q_error_over_r4": qerror,
            "g_error_over_r4": F(323+2*296),
            "stiffness_lower_times_r2": F(80)/(8+F(1, 100))-161*r2,
            "stiffness_endpoint_upper_times_r2": 10+161*r2,
            "RMS_over_stiffness_proxy_lower": F(3, 2)/F(13, 4)}


def checks():
    polynomials = positive_envelope_polynomials()
    if any(value < 0 for poly in polynomials.values() for value in poly):
        raise ValueError("The continuous monotone enclosure lost positivity")
    # Expand D=(2+delta)(1+v)^4-2 directly from binomial coefficients.
    d_base = [2*comb(4, j) for j in range(5)]
    d_delta = [comb(4, j) for j in range(5)]
    if d_base != [2, 8, 12, 8, 2] or d_delta != [1, 4, 6, 4, 1]:
        raise ValueError("The exact positive denominator coefficients changed")
    n_v, n_delta = F(-1008), F(-32)
    punctured = n_v/8-F(80)*12/64
    if punctured != -141 or n_delta != -32:
        raise ValueError("The mass remainder path controls changed")
    # Literal source/map center normalization after extracting sqrt(5).
    ag_sqrt5, bg_sqrt5 = F(2), F(-4)
    jl_sqrt5, jh_sqrt5 = F(1), F(-2)
    full = (ag_sqrt5*jl_sqrt5+bg_sqrt5*jh_sqrt5)/5
    projected = ag_sqrt5**2/10
    if full != 2 or projected != F(2, 5) or full-projected != F(8, 5):
        raise ValueError("The physical g source is not normalized")
    # A finite-dimensional Schur sign/source control, not a physical oscillator fit.
    ll = (F(1)-F(2, 5)*(-2))/(3-F(4, 5))
    hh = (F(-2)-2*ll)/5
    if 3*ll+2*hh != 1 or 2*ll+5*hh != -2:
        raise ValueError("The Schur source or feedback sign changed")
    return {"positive_polynomial_coefficients": sum(map(len, polynomials.values())),
            "mass_center_remainder": n_delta, "mass_punctured_remainder": punctured,
            "actual_center_retarded_source": full, "prepared_center_projection": projected,
            "complementary_center_projection": full-projected,
            "source_aware_schur_control": (ll, hh)}
