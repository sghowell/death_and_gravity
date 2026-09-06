"""Exact Maxwell proper H2 functional and compact-support IBP controls."""

import sympy as sp

from . import stress


def proper_operator(sampler, hubble, clock):
    f, h, t = map(sp.sympify, (sampler, hubble, clock))
    return sp.diff(f, t, 2)-2*h*sp.diff(f, t)+(3*h**2/4-3*sp.diff(h, t)/2)*f


def ibp_coefficients(hubble, clock):
    """Integral Lf² = integral f''² + first*f'² + zero*f²; zero boundary jets."""
    h, t = map(sp.sympify, (hubble, clock))
    p, q = -2*h, 3*h**2/4-3*sp.diff(h, t)/2
    return {"first": sp.expand(p*p-sp.diff(p, t)-2*q),
            "zero": sp.expand(q*q+sp.diff(q, t, 2)-sp.diff(p*q, t))}


def radiation_control():
    t = sp.Symbol("t", positive=True)
    h = 1/(2*t)
    data = stress.reference_jets(h, *(sp.diff(h, t, j) for j in (1, 2, 3)))
    coefficients = ibp_coefficients(h, t)
    return {"proper_time_positive": True,
            "operator_first_times_t": sp.simplify(-2*h*t),
            "operator_zero_times_t_squared": sp.simplify((3*h*h/4-3*sp.diff(h, t)/2)*t*t),
            "squared_first_times_t_squared": sp.simplify(coefficients["first"]*t*t),
            "squared_zero_times_t_fourth": sp.simplify(coefficients["zero"]*t**4),
            "reference_rho_pi_squared_t_fourth_over_hbar": sp.simplify(data["rho"]*t**4/2880),
            "reference_EED_pi_squared_t_fourth_over_hbar": sp.simplify(data["EED"]*t**4/2880),
            "absolute_zero_times_t_fourth": sp.simplify((coefficients["zero"]-data["EED"]/360)*t**4),
            "finite_ambiguity_in_density": sp.simplify(sp.diff(data["rho"], stress.BETA_M)),
            "finite_ambiguity_in_EED": sp.simplify(sp.diff(data["EED"], stress.BETA_M))}


def identities():
    t = sp.Symbol("t", real=True)
    f, h, a = (sp.Function(name)(t) for name in ("f", "H", "a"))
    p, q = -2*h, 3*h*h/4-3*sp.diff(h, t)/2
    coeff = ibp_coefficients(h, t)
    boundary = p*sp.diff(f, t)**2+2*q*f*sp.diff(f, t)+(p*q-sp.diff(q, t))*f**2
    integrand = proper_operator(f, h, t)**2-sp.diff(f, t, 2)**2
    integrand -= coeff["first"]*sp.diff(f, t)**2+coeff["zero"]*f**2
    conformal_second = a*sp.diff(a*sp.diff(a**-sp.Rational(3, 2)*f, t), t)
    scale, flat_energy, tau_measure = sp.symbols("scale flat_energy dt", positive=True)
    frequency = sp.Symbol("v", positive=True)
    upper = sp.Symbol("u", positive=True)
    # FP97: 4/(2pi)^3 times int_0^u v^3 dv; real Parseval adds pi.
    spectral_coefficient = 4/(2*sp.pi)**3*sp.integrate(frequency**3, (frequency, 0, upper))/upper**4*sp.pi
    return {
        "proper_clock_second_derivative": sp.simplify(conformal_second/sp.sqrt(a)
            -proper_operator(f, sp.diff(a, t)/a, t)),
        "exact_compact_support_IBP": sp.simplify(integrand-sp.diff(boundary, t)),
        "Maxwell_two_polarization_coefficient": spectral_coefficient-1/(8*sp.pi**2),
        "conformal_energy_measure_weight": sp.simplify(
            scale**-4*flat_energy*scale*tau_measure-scale**-3*flat_energy*tau_measure),
        "flat_proper_operator": sp.simplify(proper_operator(f, 0, t)-sp.diff(f, t, 2)),
        "radiation_absolute_coefficient": radiation_control()["absolute_zero_times_t_fourth"]-sp.Rational(4601, 1280),
    }
