"""Compact rational coefficients with one exact rolling-matter parity factor.

l denotes the compact matter velocity and satisfies l^2=(1-x^2)^11/100.
z=1/q compactifies high momentum; it is not the canonical coordinate zeta.
"""

from functools import cache

import sympy as sp
from p8_m1 import nonlinear, series
from p8_m1_physical import quadratic

x = series.x
z = sp.Symbol("z", nonnegative=True)
l = sp.Symbol("l", real=True)
L2 = (1-x*x)**11/100
H = 4*x
THETA = x*(4-(1-x*x)**3)
LAMBDA = 1-sp.Rational(3, 2)*(1-x*x)**3
J = series.compact_even(nonlinear.d*nonlinear.functions()["J"])
J0 = sp.factor(J+L2*LAMBDA**2/2)
FIELD = sp.QQ.frac_field(x, z)
LX, LZ = FIELD.gens
MATTER_SQUARED = FIELD.from_sympy(L2)


def reduce_l(polynomial):
    return sum(coefficient*L2**(powers[0]//2)*l**(powers[0] % 2)
               for powers, coefficient in sp.Poly(polynomial, l).terms())


def multiply(left, right):
    a, b = left
    c, d = right
    return a*c+MATTER_SQUARED*b*d, a*d+b*c


@cache
def parity(value):
    """Evaluate the expression tree in QQ(x,z)[l]/(l^2-L2).

    Polynomial arithmetic avoids the combinatorial multinomial expansion of
    large denominator powers produced by general-purpose expression cancel.
    This is exact arithmetic, not a truncated series or interpolation.
    """
    value = sp.sympify(value)
    zero, one = FIELD.zero, FIELD.one
    if value == l:
        return zero, one
    if not value.has(l):
        return FIELD.from_sympy(value), zero
    if value.is_Add:
        terms = [parity(arg) for arg in value.args]
        return sum((a for a, _ in terms), zero), sum((b for _, b in terms), zero)
    if value.is_Mul:
        result = one, zero
        for arg in value.args:
            result = multiply(result, parity(arg))
        return result
    if value.is_Pow and value.exp.is_Integer:
        a, b = parity(value.base)
        n = int(value.exp)
        if not b:
            return a**n, zero
        if n < 0:
            denominator = a*a-MATTER_SQUARED*b*b
            a, b, n = a/denominator, -b/denominator, -n
        result, base = (one, zero), (a, b)
        while n:
            if n % 2:
                result = multiply(result, base)
            n //= 2
            if n:
                base = multiply(base, base)
        return result
    raise ValueError(f"Unsupported coefficient operation: {value.func}")


def canonical(value):
    a, b = parity(value)
    return FIELD.to_sympy(a)+l*FIELD.to_sympy(b)


def derivative(value, dimension=0):
    """ell^(dimension+1) d_t[value/ell^dimension], fixed comoving momentum."""
    a, b = parity(value)
    def even_derivative(coefficient):
        return (1-LX*LX)*coefficient.diff(LX)+6*LX*LZ*coefficient.diff(LZ)
    a = even_derivative(a)-dimension*LX*a
    b = even_derivative(b)-(dimension+11)*LX*b
    return FIELD.to_sympy(a)+l*FIELD.to_sympy(b)


@cache
def coefficients(chart):
    base = quadratic.symbolic(chart)
    mapping = {quadratic.H: H, quadratic.l: l, quadratic.theta: THETA,
               quadratic.lam: LAMBDA, quadratic.w: -l*LAMBDA, quadratic.J: J,
               quadratic.q: 1/z}
    C = sp.hessian(base["density"], quadratic.Q)
    gamma = (C-base["B"].T*base["alpha"]*base["B"]).applyfunc(sp.factor)
    return {key: matrix.subs(mapping, simultaneous=True).applyfunc(canonical)
            for key, matrix in {"alpha": base["alpha"], "beta": base["beta"],
                                "B": base["B"], "C": C, "gamma": gamma}.items()}


@cache
def factors(chart):
    if chart == "unitary":
        d1, d2, chi, delta = 2*J/THETA**2, sp.Integer(1), -l*LAMBDA/THETA, 0
    elif chart == "gamma":
        R, D = LAMBDA**2-J*z, LAMBDA**2-J0*z
        d1, d2, chi, delta = 2*J/R, R/D, l*LAMBDA**2/R, 1
    else:
        raise ValueError("Use unitary or gamma scalar chart")
    d1, d2, chi = map(canonical, (d1, d2, chi))
    ratio = canonical(d2/d1)
    f1 = canonical(6*x+derivative(d1)/(2*d1)-delta*x)
    f2 = canonical(6*x+derivative(d2)/(2*d2))
    h = derivative(chi, delta)
    return {"d1": d1, "d2": d2, "chi": chi, "delta": delta, "ratio": ratio,
            "f1": f1, "f2": f2, "h": h,
            "log_root_ratio_derivative": canonical(derivative(ratio)/(2*ratio))}


def no_high_frequency_pole(value):
    """Exact z=0 denominator test; no sampled extrapolation to high frequency."""
    for coefficient in parity(value):
        if coefficient.denom.as_expr().subs(z, 0) == 0:
            raise ValueError("Uncancelled high-frequency pole")
    return True
