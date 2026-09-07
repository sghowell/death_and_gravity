"""Independent exact Fraction profile jets and literal coframe variations.

No SymPy, primary module, precomputed verdict or interpolated sign proof is
imported. Finite fixtures verify identities; the interval conclusion uses
the continuous written inequalities and independently rebuilt coefficients.
"""

from fractions import Fraction as Q
from math import comb


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Q)):
        raise TypeError("Only exact Python integers/Fractions are accepted")
    return Q(value)


def profile(time, action_lapse, ratio=2):
    u, c, r = map(rational, (time, action_lapse, ratio))
    if abs(u) > Q(1, 10) or not (c == 1 or 2 < c <= 4) or r <= 0:
        raise ValueError("Outside the specified local coefficient/rational fixture domain")
    d = 1+u*u
    a, y = d*d, 2/d**4
    h, hp = 4*u/d, 4*(1-u*u)/d**2
    yp, hpp = -16*u/d**5, -8*u*(3-u*u)/d**3
    n = 2*(y**3/c-1)*hp
    np = 6*y*y*yp*hp/c+2*(y**3/c-1)*hpp
    numerator, denominator = y**3*hp, c*(c-y)
    numerator_p, denominator_p = 3*y*y*yp*hp+y**3*hpp, -c*yp
    b1 = numerator/denominator
    b1p = (numerator_p*denominator-numerator*denominator_p)/denominator**2
    b4 = 3*h*h/(2*c*c)-b1/y**3
    b4p = 3*h*hp/(c*c)-b1p/y**3+3*b1*yp/y**4
    b0 = 3*h*h/2-n/4-3*b1*y
    b0p = 3*h*hp-np/4-3*b1p*y-3*b1*yp
    return {"a": a, "y": y, "h": h, "h_u": hp, "nbar": n,
            "kbar": n-1/(100*a**6), "b0": b0, "b1": b1, "b4": b4,
            "b0_u": b0p, "b1_u": b1p, "b4_u": b4p,
            "theta": 3*h*h*(c-y)/(2*c*hp), "r_star_cubed": -b1/b4,
            "g_flat_residual": b0+3*b1*r, "f_flat_residual": b1+b4*r**3,
            "clock_profile_residual": b0p+4*r*b1p+r**4*b4p,
            "metric_eliminant": b0**3*b4-27*b1**4}


def fixtures():
    return [(Q(0), Q(1), Q(2)), (Q(1, 10), Q(1), Q(2)),
            (Q(-1, 10), Q(1), Q(3, 2)), (Q(0), Q(4), Q(2)),
            (Q(1, 10), Q(4), Q(9, 4)), (Q(-1, 10), Q(3), Q(2)),
            (Q(0), Q(201, 100), Q(2)), (Q(1, 10), Q(201, 100), Q(2))]


def density(ag, af, ng, nf, betas):
    """The literal elementary generating polynomial in four coframes."""
    result = Q(0)
    for j in range(5):
        if j <= 3:
            result += betas[j]*comb(3, j)*ng*ag**(3-j)*af**j
        if j >= 1:
            result += betas[j]*comb(3, j-1)*nf*ag**(4-j)*af**(j-1)
    return -2*result


def derivative_cubic(function, value):
    step = value/10
    return (function(value-2*step)-8*function(value-step)
            +8*function(value+step)-function(value+2*step))/(12*step)


def coframe_fixture(ag, ng, ratio, betas, slopes):
    ag, ng, ratio = map(rational, (ag, ng, ratio))
    betas, slopes = tuple(map(rational, betas)), tuple(map(rational, slopes))
    if ag <= 0 or ng <= 0 or ratio <= 0 or len(betas) != 5 or len(slopes) != 5:
        raise ValueError("Require regular proportional coframes and five coefficient/slope pairs")
    af, nf = ratio*ag, ratio*ng
    dng = derivative_cubic(lambda value: density(ag, af, value, nf, betas), ng)
    dnf = derivative_cubic(lambda value: density(ag, af, ng, value, betas), nf)
    dag = derivative_cubic(lambda value: density(value, af, ng, nf, betas), ag)
    daf = derivative_cubic(lambda value: density(ag, value, ng, nf, betas), af)
    # A literal linear clock variation is exact, with coframes held fixed.
    plus, minus = tuple(a+b for a, b in zip(betas, slopes)), tuple(a-b for a, b in zip(betas, slopes))
    dphi = (density(ag, af, ng, nf, plus)-density(ag, af, ng, nf, minus))/2
    return {"rho_g": -dng/ag**3, "rho_f": -dnf/af**3,
            "pressure_g": dag/(3*ng*ag**2), "pressure_f": daf/(3*nf*af**2),
            "clock_equation": -dphi/(ng*ag**3), "W": -density(ag, af, ng, nf, betas)/(2*ng*ag**3)}


def coframe_fixtures():
    return [
        (Q(1), Q(1), Q(2), (Q(-26), Q(4), Q(0), Q(0), Q(-1, 2)), (Q(0),)*5),
        (Q(3, 2), Q(2, 3), Q(5, 4), (Q(1, 3), Q(-2), Q(5, 7), Q(4, 3), Q(-9, 5)),
         (Q(1), Q(-3, 2), Q(4), Q(2, 7), Q(5, 9))),
        (Q(7, 3), Q(5, 6), Q(3, 2), (Q(2), Q(3), Q(-4), Q(5), Q(7)),
         (Q(-4), Q(6), Q(0), Q(3), Q(1))),
    ]


def coefficients():
    """Rebuild the polynomial bound from the primitive rational constants."""
    h2, hp, y3, ymax, rmax = Q(4, 25), Q(3), Q(7), Q(2), Q(9, 4)
    T = Q(3, 2)*h2
    polynomial = (T*ymax**2*rmax**2, -hp*y3/2, T+hp/2)
    at2 = sum(value*2**power for power, value in enumerate(polynomial))
    margin_c = (at2-polynomial[0], -polynomial[1], -polynomial[2])
    # c=2+2z, followed by the degree-two power-to-Bernstein transform.
    power_z = [sum(margin_c[j]*comb(j, i)*2**j for j in range(i, 3)) for i in range(3)]
    bernstein = (power_z[0], power_z[0]+power_z[1]/2, sum(power_z))
    return {"polynomial_power_c": polynomial, "upper_at_c2": at2,
            "margin_Bernstein": bernstein, "c1_U_upper": T*(1+ymax**4)-9,
            "positive_branch_U_upper": at2/16}
