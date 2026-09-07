"""Exact feasibility algebra for a different, variable-lapse parent.

This is not S6.20's constant-lapse action and has no certificate yet.
Physical g time is u=T/tau; dimensionless checks use M=tau=1.
"""

from functools import cache

import sympy as sp


def profiles(u):
    """The derivative of z=-tau*H_f is explicit; z itself is its integral."""
    if isinstance(u, bool):
        raise TypeError("An exact real input or real symbol is required")
    u = sp.sympify(u)
    if u.is_real is not True or u.is_finite is False or u.has(sp.Float):
        raise ValueError("An exact finite real input or real symbol is required")
    d = 1+u**2
    a, b = d**2, 2/d**2
    y, h = b/a, 4*u/d
    hp = 4*(1-u**2)/d**2
    aa = (1-u**2)*d**10/2
    rad = sp.sqrt(aa**2+2)
    zp = (aa+rad)/2
    w = 1/(10*a**3)
    n = y**3*(rad-aa)
    return {"d": d, "a": a, "b": b, "y": y, "h": h, "h_prime": hp,
            "A": aa, "z_prime": zp, "chi_speed": w,
            "null_stress": n, "clock_squared": n-w**2}


@cache
def reconstruction_checks():
    """All metric and scalar equations for arbitrary regular z(u).

Keeping z and z' as independent jets avoids replacing the canonical clock
by an externally time-dependent coefficient. y'=−2hy and a'=ah are kept.
"""
    a, h, hp, hpp, z, zp, zpp = sp.symbols("a h hp hpp z zp zpp", nonzero=True)
    y, c = 2/a**2, h/z

    def derivative(value):
        return (a*h*sp.diff(value, a)+hp*sp.diff(value, h)
                +hpp*sp.diff(value, hp)+zp*sp.diff(value, z)+zpp*sp.diff(value, zp))

    n = 2*(y**3*zp-hp)
    w = 1/(10*a**3)
    k = n-w**2
    b1 = y**3*zp/(c-y)
    b4 = 3*z**2/2-b1/y**3
    b0 = (3*h**2-n/2)/2-3*b1*y
    p = 2*b1
    result = {
        "g_lapse": 3*h**2-n/2-2*(b0+3*b1*y),
        "f_lapse": 3*z**2-2*(b1+b4*y**3)/y**3,
        "g_null": -2*hp-n-(y-c)*p,
        "f_null": 2*zp/c-(c-y)*p/(c*y**3),
        "free_chi": derivative(w)+3*h*w,
        "sourced_clock_times_speed": derivative(k)/2+3*h*k
        +2*(derivative(b0)+(c+3*y)*derivative(b1)+c*y**3*derivative(b4)),
        "sourced_Bianchi": 3*p*(-y*z-h)-2*(derivative(b1)+y**3*derivative(b4)),
    }
    return {name: sp.factor(value) for name, value in result.items()}


def bernstein_margin_polynomial():
    """B(B−A)>1/2 implies z'<B=(h/y)' on 0<=u²<=1."""
    x = sp.Symbol("x", real=True)
    aa = (1-x)*(1+x)**10/2
    bb = 2*(1+x)**2*(1+7*x)
    return sp.Poly(bb*(bb-aa)-sp.Rational(1, 2), x)


def bernstein_cells():
    """Exact degree-14 Bernstein enclosure, not point sampling."""
    poly = bernstein_margin_polynomial()
    x, degree = poly.gen, poly.degree()
    t = sp.Symbol("t", real=True)
    result = []
    for cell in range(16):
        local = sp.Poly(poly.as_expr().subs(x, (cell+t)/16), t)
        coefficients = tuple(sum(local.nth(i)*sp.binomial(j, i)/sp.binomial(degree, i)
                                 for i in range(j+1)) for j in range(degree+1))
        result.append(coefficients)
    return tuple(result)


def center_jets():
    """Taylor checks for the analytic removable u=0 lapse, no 0/0 verdict."""
    u = sp.Symbol("u", real=True)
    d = profiles(u)
    zp = sp.series(d["z_prime"], u, 0, 5).removeO()
    z = sp.integrate(zp, (u, 0, u))
    c = sp.series(d["h"]/z, u, 0, 3).removeO()
    return {"z_prime_0": zp.subs(u, 0), "z_third_0": sp.diff(z, u, 3).subs(u, 0),
            "c0": c.subs(u, 0), "c_second_0": sp.diff(c, u, 2).subs(u, 0),
            "clock_squared_0": d["clock_squared"].subs(u, 0),
            "clock_squared_second_0": sp.diff(d["clock_squared"], u, 2).subs(u, 0)}


def tail_checks():
    """Rationalized infinity limits avoid catastrophic radical subtraction."""
    r = sp.Symbol("r", positive=True)  # r=1/u on the future tail
    aa = (r**2-1)*(r**2+1)**10/(2*r**22)
    zp = 1/(sp.sqrt(aa**2+2)-aa)
    y, h = 2*r**8/(1+r**2)**4, 4*r/(1+r**2)
    z_inf = sp.Symbol("z_infinity", positive=True)
    c = h/z_inf
    mass = 2*y*zp*(y**3+c)/(c-y)
    n = 2*(y**3*zp-4*(r**2-1)*r**2/(1+r**2)**2)
    return {"u22_zprime": sp.limit(zp/r**22, r, 0),
            "u2_null_stress": sp.limit(n/r**2, r, 0),
            "u30_algebraic_mass_squared": sp.limit(mass/r**30, r, 0),
            "u_c": sp.limit(c/r, r, 0),
            "f_Hubble_tail": -z_inf}
