"""Universal *retarded memory* kernel; not an instantaneous mass inverse."""

from functools import cache

import sympy as sp

from .exact import MU, R_MAX, number, parameters


def jost(z):
    argument = (1-sp.tanh(z))/2
    return sp.exp(sp.I*MU*z)*sp.hyper((-sp.Rational(1, 2), sp.Rational(3, 2)), (1-sp.I*MU,), argument)


def universal_kernel(output, source, delta, r=R_MAX, momentum_squared=1):
    """Real kernel G0(u,s), zero for u<s, using the actual finite delta clock.

    The K argument enforces the physical family but the universal pole kernel
    is K-independent. No numerical branch continuation is substituted.
    """
    delta, r, _ = parameters(delta, r, momentum_squared)
    uu, ss = number(output, "u"), number(source, "s")
    if any((r-abs(value)).is_nonnegative is not True for value in (uu, ss)):
        raise ValueError("Both retarded-kernel times must lie in [-r,r]")
    if (ss-uu).is_positive is True or uu == ss:
        return sp.S.Zero
    if (uu-ss).is_positive is not True:
        raise ValueError("The retarded ordering must be decidable")
    ru, rs = sp.sqrt(uu**2+delta/8), sp.sqrt(ss**2+delta/8)
    zu, zs = sp.asinh(sp.sqrt(8)*uu/sp.sqrt(delta)), sp.asinh(sp.sqrt(8)*ss/sp.sqrt(delta))
    return sp.sqrt(ru*rs)*sp.im(jost(zu)*sp.conjugate(jost(zs)))/MU


@cache
def checks():
    u, delta = sp.symbols("u delta", positive=True)
    rr = sp.sqrt(u**2+delta/8)
    z = sp.asinh(sp.sqrt(8)*u/sp.sqrt(delta))
    psi = sp.Function("psi")
    b, force = sp.symbols("b F", real=True)
    qq = sp.sqrt(rr)*psi(z)
    transformed = sp.expand(rr**sp.Rational(3, 2)*(sp.diff(qq, u, 2)+(10/rr**2+b)*qq-force))
    target = (sp.Subs(sp.Derivative(psi(sp.Symbol("z")), (sp.Symbol("z"), 2)), sp.Symbol("z"), z)
              +(MU**2+sp.Rational(3, 4)*(1-u**2/rr**2)+rr**2*b)*psi(z)-rr**sp.Rational(3, 2)*force)
    zz = sp.Symbol("z", real=True)
    xx = (1-sp.tanh(zz))/2
    # Hypergeometric substitution by algebraic differential coefficients,
    # avoiding any assumed asymptotic numerical evaluation of 2F1.
    f, fp, fpp = sp.symbols("F Fp Fpp")
    aa, bb, cc = -sp.Rational(1, 2), sp.Rational(3, 2), 1-sp.I*MU
    hyper_second = (aa*bb*f-(cc-(aa+bb+1)*xx)*fp)/(xx*(1-xx))
    direct = ((-MU**2+MU**2+sp.Rational(3, 4)/sp.cosh(zz)**2)*f
              +(2*sp.I*MU*sp.diff(xx, zz)+sp.diff(xx, zz, 2))*fp
              +sp.diff(xx, zz)**2*fpp)
    # For Q=sqrt(rho)psi, the Wronskian is unchanged since z'=1/rho.
    return {"actual_Liouville_chain_rule": sp.simplify(transformed-target),
            "hypergeometric_to_Poschl_Teller": sp.simplify(sp.expand_trig(direct.subs(fpp, hyper_second)).rewrite(sp.exp)),
            "Wronskian_clock_cancellation": sp.simplify(rr*sp.diff(z, u)-1),
            "positive_retarded_jump": sp.simplify((sp.I*MU-(-sp.I*MU))/(2*sp.I*MU)-1),
            "pole_mu": MU**2+sp.Rational(1, 4)-10}
