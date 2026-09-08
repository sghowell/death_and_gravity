"""Independent S4 spectral check with a bounded Euler-Maclaurin remainder."""
from functools import cache

import sympy as sp

from . import heat


@cache
def spectral():
    ell = sp.Symbol("ell", integer=True, positive=True)
    t, x = sp.symbols("proper_time spectral_x", positive=True)
    eigenvalue = (ell+1)*(ell+2)
    multiplicity = ell*(ell+3)*(2*ell+3)/2
    f = (x**3-sp.Rational(9, 4)*x)*sp.exp(-t*x**2)
    half = sp.Rational(1, 2)
    lower_integral = sp.exp(-t/4)*(1/(2*t**2)-1/t)
    endpoint = f.subs(x, half)/2-sum(sp.bernoulli(2*j)/sp.factorial(2*j)
                        *sp.diff(f, x, 2*j-1).subs(x, half) for j in range(1, 4))
    approximation = sp.expand((lower_integral+endpoint)*sp.exp(t/4))
    leading = 1/(2*t**2)-1/t-sp.Rational(11, 30)
    polynomial_tail = sp.Poly(sp.expand(approximation-leading), t)
    if polynomial_tail.nth(0) != 0:
        raise ValueError("The independent spectral constant disagrees")
    polynomial_bound = sum(abs(coefficient) for coefficient in polynomial_tail.all_coeffs())
    y = sp.Symbol("y", nonnegative=True)
    derivative_norms = {}
    for power in (1, 3):
        shape = sp.Poly(sp.expand(sp.diff(y**power*sp.exp(-y**2), y, 6)*sp.exp(y**2)), y)
        derivative_norms[power] = sum(abs(coefficient)*sp.gamma(sp.Rational(exponent[0]+1, 2))/2
                                     for exponent, coefficient in shape.terms())
    # Integral |d^6[x^3 exp(-t x²)]| <= C3*t, and the x term <= C1*t².
    em_constant = (derivative_norms[3]+sp.Rational(9, 4)*derivative_norms[1])/30240
    total_error_constant = sp.factor(polynomial_bound+sp.Rational(4, 3)*em_constant)
    local = heat.coefficients()["proca"]
    sphere_values = {heat.R: 12, heat.Ric2: 36, heat.Riem2: 24, heat.boxR: 0}
    expected = sum(local[order].subs(sphere_values)*t**(order//2-2)/6 for order in (0, 2, 4))
    bernoulli_variable = sp.Symbol("bernoulli_x", real=True)
    b = bernoulli_variable*(1-bernoulli_variable)
    return {"t": t, "eigenvalue": eigenvalue, "multiplicity": multiplicity,
            "approximation": approximation, "leading": leading, "polynomial_tail_bound": polynomial_bound,
            "sixth_derivative_norm_envelopes": derivative_norms, "Euler_Maclaurin_bound_before_exponential": em_constant,
            "trace_remainder_over_t_upper_for_0_lt_t_le_1": total_error_constant,
            "S4_Hodge_eigenvalue_half_shift": sp.expand(eigenvalue-(ell+sp.Rational(3, 2))**2+sp.Rational(1, 4)),
            "S4_vector_multiplicity_half_shift": sp.expand(multiplicity-((ell+sp.Rational(3, 2))**3-sp.Rational(9, 4)*(ell+sp.Rational(3, 2)))),
            "scalar_constant_mode_subtraction": sp.simplify(f.subs(x, half)*sp.exp(t/4)+1),
            "vanishing_extra_half_shift_mode": sp.simplify(f.subs(x, sp.Rational(3, 2))),
            "independent_S4_heat_coefficients": sp.expand(expected-leading),
            "Euler_Maclaurin_spectral_constant": polynomial_tail.nth(0),
            "periodic_Bernoulli_six_bound_polynomial": sp.expand(sp.bernoulli(6, bernoulli_variable)
                                                        -(sp.Rational(1, 42)-b**2/2-b**3))}


@cache
def checks():
    d = spectral()
    names = ("S4_Hodge_eigenvalue_half_shift", "S4_vector_multiplicity_half_shift", "scalar_constant_mode_subtraction",
             "vanishing_extra_half_shift_mode", "independent_S4_heat_coefficients", "Euler_Maclaurin_spectral_constant",
             "periodic_Bernoulli_six_bound_polynomial")
    return {name: d[name] for name in names}


@cache
def proof_checks():
    d = spectral()
    return {"periodic_B6_upper": sp.bernoulli(6) == sp.Rational(1, 42),
            "periodic_B6_lower_above_negative_upper": bool(sp.Rational(1, 42)-sp.Rational(1, 16)/2-sp.Rational(1, 64) > -sp.Rational(1, 42)),
            "sixth_derivative_EM_prefactor": sp.Rational(1, 42)/sp.factorial(6) == sp.Rational(1, 30240),
            "independent_linear_Gaussian_sixth_derivative_envelope": d["sixth_derivative_norm_envelopes"][1] == 2124,
            "independent_cubic_Gaussian_sixth_derivative_envelope": d["sixth_derivative_norm_envelopes"][3] == 14016,
            "Euler_Maclaurin_polynomial_tail_envelope": d["polynomial_tail_bound"] == sp.Rational(59, 576),
            "Euler_Maclaurin_derivative_tail_envelope": d["Euler_Maclaurin_bound_before_exponential"] == sp.Rational(179, 288),
            "unit_proper_time_exponential_envelope": 1/(1-sp.Rational(1, 4)) == sp.Rational(4, 3),
            "complete_S4_trace_error_envelope": d["trace_remainder_over_t_upper_for_0_lt_t_le_1"] == sp.Rational(1609, 1728)}
