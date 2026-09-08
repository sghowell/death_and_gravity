"""Dimension-dependent WKB coefficients and exact radial Laurent terms."""
from functools import cache

import sympy as sp

dimension = sp.Symbol("spatial_dimension", real=True)
z = sp.Symbol("physical_momentum_fraction", real=True)
H, Hd, Hdd, Hddd, Hdddd = sp.symbols("H H_dot H_ddot H_third H_fourth", real=True)
alpha, beta, ell = sp.symbols("a_N b_N log_mass_over_scale_squared", real=True)


def derivative(expression):
    jets = ((H, Hd), (Hd, Hdd), (Hdd, Hddd), (Hddd, Hdddd))
    return sp.expand(sum(sp.diff(expression, before)*after for before, after in jets)
                     -2*H*z*(1-z)*sp.diff(expression, z))


@cache
def reference(kind):
    s = (dimension-2)/2
    lam = -H*z
    if kind == "transverse":
        rate = s*H
    elif kind == "longitudinal":
        rate = (s+z)*H
    else:
        raise ValueError("Require transverse or longitudinal mode")
    U = derivative(rate)+rate**2
    p2 = sp.expand(-U/2-derivative(lam)/4+lam**2/8)
    p4 = sp.expand(-p2**2/2-derivative(derivative(p2))/4
                   +sp.Rational(5, 4)*lam*derivative(p2)
                   +(derivative(lam)/2-sp.Rational(3, 2)*lam**2)*p2)
    b2 = derivative(p2)-2*lam*p2
    return {"rate": rate, "U": sp.expand(U), "P2": p2, "P4": p4, "B2": sp.expand(b2), "c1": rate+lam/2}


@cache
def energy_coefficients():
    out = {}
    for kind in ("transverse", "longitudinal"):
        data = reference(kind)
        if kind == "transverse":
            A, B = sp.Integer(1), 1+beta*(1-z)
        else:
            A, B = 1-alpha*z, 1+beta
        p, r, b2, c1 = [data[key] for key in ("P2", "P4", "B2", "c1")]
        out[kind] = {0: sp.expand(A+B),
                     1: sp.expand((A-B)*p+A*c1**2),
                     2: sp.expand((A-B)*r+B*p**2+A*(c1*b2-c1**2*p))}
    out["combined"] = {order: sp.expand((dimension-1)*out["transverse"][order]+out["longitudinal"][order])
                       for order in range(3)}
    return out


@cache
def pressure_coefficients():
    out = {}
    for kind in ("transverse", "longitudinal"):
        data = reference(kind)
        if kind == "transverse":
            A, B = (dimension-2)/dimension, (2*z-dimension+2)/dimension
        else:
            A, B = (dimension-2+2*z)/dimension, -(dimension-2)/dimension
        p, r, b2, c1 = [data[key] for key in ("P2", "P4", "B2", "c1")]
        out[kind] = {0: sp.expand(A+B),
                     1: sp.expand((A-B)*p+A*c1**2),
                     2: sp.expand((A-B)*r+B*p**2+A*(c1*b2-c1**2*p))}
    out["combined"] = {order: sp.expand((dimension-1)*out["transverse"][order]+out["longitudinal"][order])
                       for order in range(3)}
    return out


def radial_polynomial(coefficient, order):
    """Radial angular coefficient relative to m^(4-2n)/(64*pi²).

    Gamma(epsilon+n-2) is left out here and restored in laurent().
    Dimensional regularization uses e^(gamma*epsilon) (mu²/m²)^epsilon.
    """
    if type(order) is not int or order not in (0, 1, 2):
        raise ValueError("Require adiabatic half-order zero, one or two")
    if not isinstance(coefficient, sp.Expr) or coefficient.has(sp.Float, sp.oo, -sp.oo, sp.zoo, sp.nan):
        raise ValueError("Require an exact finite symbolic coefficient")
    try:
        polynomial = sp.Poly(coefficient, z)
    except sp.PolynomialError as error:
        raise ValueError("Require a polynomial in the physical momentum fraction") from error
    result = sp.Integer(0)
    for (power,), value in polynomial.terms():
        result += 2*sp.sqrt(sp.pi)*value*sp.rf(dimension/2, power)/sp.gamma(order+power-sp.Rational(1, 2))
    return sp.expand(result)


def laurent(coefficient, order):
    radial = radial_polynomial(coefficient, order)
    residue_factor = sp.Rational((-1)**(2-order), sp.factorial(2-order))
    pole = sp.factor(residue_factor*radial.subs(dimension, 3))
    evanescent = sp.factor(-2*residue_factor*sp.diff(radial, dimension).subs(dimension, 3))
    finite = sp.factor(pole*(sp.harmonic(2-order)-ell)+evanescent)
    # This is a component-wise radial finite part. Covariant counterterm
    # variation in d dimensions can add a further finite local term.
    return {"pole": pole, "radial_finite_part": finite, "evanescent_finite_contribution": evanescent}


@cache
def integrated():
    data = energy_coefficients()["combined"]
    return {order: laurent(data[order], order) for order in range(3)}


@cache
def integrated_pressure():
    data = pressure_coefficients()["combined"]
    return {order: laurent(data[order], order) for order in range(3)}


@cache
def checks():
    data = integrated()
    expected_pole0 = -3-sp.Rational(3, 2)*alpha-sp.Rational(9, 2)*beta
    expected_finite0 = -expected_pole0*ell-sp.Rational(5, 2)-sp.Rational(5, 4)*alpha-sp.Rational(3, 4)*beta
    frozen_dimension = laurent(energy_coefficients()["combined"][0].subs(dimension, 3), 0)
    return {"flat_actual_mass_direction_pole": sp.factor(data[0]["pole"]-expected_pole0),
            "flat_actual_mass_direction_finite_MSbar": sp.factor(data[0]["radial_finite_part"]-expected_finite0),
            "early_polarization_count_limit_loses_finite_term": sp.factor(data[0]["radial_finite_part"]
                    -frozen_dimension["radial_finite_part"]-2-4*beta),
            "ordinary_Proca_R_pole_lapse_variation": sp.factor(data[1]["pole"].subs({alpha: 0, beta: 0})+6*H**2),
            "ordinary_Proca_curvature_squared_pole_lapse_variation": sp.factor(data[2]["pole"].subs({alpha: 0, beta: 0})
                                                                    -(Hd**2-6*H**2*Hd-2*H*Hdd))}
