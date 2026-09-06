"""Exact proper/conformal plateau geometry and literal 1922 denominator."""

import sympy as sp
from p8a_backreaction import dynamics as frozen

Y, DELTA = sp.symbols("y delta", positive=True)


def euler_ratios(z):
    z = sp.sympify(z)
    return {"density": (1-z)**-4,
            "pressure": (1+8*z/5)*(1-z)**-4,
            "EED": (1+4*z/3)*(1-z)**-4}


def euler_error_constants():
    # Derivative caps of the ratios on z<=1/32 are 8,12,12. Since z=delta/y^4,
    # y>=2 contributes y^-12<=1/4096, with the A.3 component normalizations.
    return {"density": sp.Rational(1, 491520), "pressure": sp.Rational(1, 196608),
            "EED": sp.Rational(3, 327680)}


def identities():
    y, delta = Y, DELTA
    f = 1-delta/y**4
    a = y*f**sp.Rational(1, 4)
    derivative = lambda value: f**sp.Rational(1, 4)*sp.diff(value, y)
    h = derivative(a)/a
    u = -derivative(derivative(a))/a
    z = delta/y**4
    normalized = {"density": h**4/(960*a**4),
                  "pressure": (5*h**4+4*h**2*u)/(2880*a**4),
                  "EED": (3*h**4+2*h**2*u)/(960*a**4)}
    denominators = {"density": 960, "pressure": 576, "EED": 320}
    results = {"exact_Hc": sp.simplify(h-1/(y*f**sp.Rational(3, 4))),
               "exact_U": sp.simplify(u-2*delta/(y**6*f**sp.Rational(3, 2)))}
    for name, ratio in euler_ratios(z).items():
        results[f"Euler_{name}_ratio"] = sp.simplify(normalized[name]-ratio/(denominators[name]*y**8))
    normalization, length, epsilon, d, kappa, hbar = sp.symbols("A eta_star epsilon d kappa hbar", positive=True)
    time = normalization*length**2*y**2/2
    for name, value in frozen.frozen_defects(z, time).items():
        c = {"density": 3, "pressure": 5, "EED": 9}[name]
        results[f"frozen_{name}_normalization"] = sp.simplify(
            value-c*delta**2*(2-z)/(normalization**2*length**4*y**12*(1-z)**2))
    amplitude = 16*epsilon*d/(normalization**2*length**4)
    results["exact_quantum_source_coupling_restoration"] = sp.simplify(
        (epsilon*kappa*hbar/sp.pi**2-2880*amplitude*normalization**2*length**4).subs(
            d, kappa*hbar/(46080*sp.pi**2)))
    z0 = sp.Symbol("z", nonnegative=True)
    results["continuous_frozen_ratio_at_zero"] = sp.limit(((1-z0)**-2-1)/z0, z0, 0)-2
    results["literal_1922_denominator"] = sp.Rational(2048, 961)/4096-sp.Rational(1, 1922)
    return results


def elementary_margins():
    floor = sp.Rational(31, 32)
    v = sp.Symbol("v", nonnegative=True)
    z = (1-v)/32
    gap = sp.expand(sp.Rational(2048, 961)*(1-z)**2-(2-z))
    return {"frozen_ratio_positive_polynomial_coefficients": list(reversed(sp.Poly(gap, v).all_coeffs())),
            "Euler_density_slope_margin": 8*floor**5-4,
            "Euler_pressure_slope_margin": 12*floor**5-sp.Rational(23, 4),
            "Euler_EED_slope_margin": 12*floor**5-sp.Rational(131, 24),
            "positive_target_radicand_floor": floor,
            "target_y_twelfth_lower_bound": sp.Integer(4096),
            "literal_denominator": sp.Integer(1922)}
