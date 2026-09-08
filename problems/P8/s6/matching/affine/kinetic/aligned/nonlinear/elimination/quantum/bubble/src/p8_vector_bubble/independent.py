"""Independent spherical moments and isotropic denominator reduction."""
from functools import cache

import sympy as sp

from . import pole


def angular_average(polynomial, coordinates, radial_squared):
    """Exact degree-four rotational average in four Euclidean dimensions."""
    out = 0
    for powers, coefficient in sp.Poly(polynomial, *coordinates).terms():
        degree = sum(powers)
        if degree > 4:
            raise ValueError("This independent angular control is specified only through degree four")
        if any(power % 2 for power in powers):
            continue
        if degree == 0:
            factor = 1
        elif degree == 2:
            factor = radial_squared/4
        elif degree == 4:
            factor = radial_squared**2/(8 if 4 in powers else 24)
        else:
            raise ValueError("Unexpected even polynomial degree")
        out += coefficient*factor
    return sp.expand(out)


@cache
def direct_shift():
    alpha, beta = sp.symbols("alpha beta", real=True)
    time, space = sp.symbols("external_time external_space", real=True)
    l = sp.Matrix(sp.symbols("loop0:4", real=True))
    k = sp.Matrix([time, space, 0, 0])
    A = sp.diag(alpha, beta, beta, beta)
    x, m2 = pole.x, pole.m2
    p, q = l-x*k, l+(1-x)*k
    original = (m2**2*sp.trace(A*A)+m2*((p.T*A*A*p)[0]+(q.T*A*A*q)[0])+(p.T*A*q)[0]**2)
    radial = sp.Symbol("loop_squared", real=True)
    averaged = angular_average(sp.expand(original), list(l), radial)
    Delta = m2+x*(1-x)*(time**2+space**2)
    radial_residue = {0: 1, 1: -2*Delta, 2: 3*Delta**2}
    actual = sum(coefficient*radial_residue[power[0]] for power, coefficient in sp.Poly(averaged, radial).terms())
    mapping = {pole.TA: alpha+3*beta, pole.TA2: alpha**2+3*beta**2,
               pole.P2: time**2+space**2, pole.PA: alpha*time**2+beta*space**2,
               pole.PA2: alpha**2*time**2+beta**2*space**2}
    expected = pole.parameter_integral()["integrand"].subs(mapping)
    return {"independent_explicit_spherical_pole_integrand": sp.expand(actual-expected),
            "alpha": alpha, "beta": beta, "time": time, "space": space,
            "explicit_integrand": sp.factor(actual)}


@cache
def isotropic():
    A = sp.Symbol("isotropic_A", real=True)
    K, m2 = pole.P2, pole.m2
    dp, dq = sp.symbols("D_p D_q", nonzero=True)
    p_dot_q = (dp+dq-2*m2-K)/2
    numerator = 4*m2**2+m2*(dp+dq-2*m2)+p_dot_q**2
    rational = ((3*m2**2+m2*K+K**2/4)/(dp*dq)
                -K*(1/dp+1/dq)/2+(dp/dq+dq/dp+2)/4)
    # By a loop-momentum shift, integral(Dp/Dq)=k²*tadpole up to
    # scaleless and odd pieces. B0 has residue 1, tadpole residue -m².
    independent = A**2*(3*m2**2+m2*K+K**2/4+K*m2/2)
    data = pole.parameter_integral()
    mapping = {pole.TA: 4*A, pole.TA2: 4*A**2, pole.PA: A*K, pole.PA2: A**2*K}
    expected = (data["constant"]+data["second"]+data["fourth"]).subs(mapping)
    return {"independent_isotropic_denominator_identity": sp.factor(numerator/(dp*dq)-rational),
            "independent_isotropic_tadpole_bubble_pole": sp.factor(independent-expected)}


@cache
def zero_momentum_potential():
    from p8_aligned_quantum import kernel, potential
    r, alpha, beta = sp.symbols("r alpha beta", real=True)
    C = potential.coefficients()["pole_weight"].subs({kernel.a: 1+alpha*r, kernel.b: 1+beta*r})
    actual = sp.diff(C, r, 2).subs(r, 0)/2
    expected = ((alpha+3*beta)**2+2*(alpha**2+3*beta**2))/8
    return {"independent_linear_mass_variation_potential_match": sp.factor(actual-expected)}


@cache
def isotropic_nonlocal():
    """All Taylor orders n>=3 from the independently reduced scalar bubble.

    B_j=integral_0^1 [x(1-x)]^j dx=j!^2/(2j+1)!.
    Divide coefficients by B_(n-2), the common sign and m2**(2-n).
    The denominator reduction leaves a scalar-bubble coefficient
    3*m2**2+m2*K+K**2/4 and a local tadpole. The latter cannot affect
    these orders. Analyticity then identifies the subtracted functions.
    """
    n = sp.Symbol("Taylor_order", integer=True, positive=True)
    b_previous = (n-1)/(2*(2*n-1))
    b_current = n*(n-1)/(4*(2*n-1)*(2*n+1))
    original = (3*b_current/n-(sp.Rational(1, 2)*b_previous+6*b_current)/(n-1)
                +(10*b_current-sp.Rational(1, 2)*b_previous)/(n-2))
    reduced = 3*b_current/n-b_previous/(n-1)+sp.Rational(1, 4)/(n-2)
    return {"independent_isotropic_finite_log_all_Taylor_orders_n_ge_3": sp.factor(original-reduced)}


@cache
def checks():
    out = {"independent_explicit_spherical_pole_integrand": direct_shift()["independent_explicit_spherical_pole_integrand"]}
    out.update(isotropic())
    out.update(zero_momentum_potential())
    out.update(isotropic_nonlocal())
    return out
