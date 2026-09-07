"""Continuous real-box coefficient bounds, including the full 0<=K<=4 band.

These are derived afresh from the literal formulas, not the S6.21 K>=1
transfer theorem. v=u^2<=r^2, delta<=r^2/100 and c,d are positive.
"""

from functools import cache

import sympy as sp

from . import operator
from .exact import R_MAX, radius


def envelope(r=R_MAX):
    r = radius(r)
    t, cp, dp = r**2, 2+r**2/100, 1+r**2
    theta = sp.Rational(18, 5)
    # |sqrt(2c)|<=2+delta/2; W=c*d^12+8>=10.
    omega = (48+12*t/100)*dp**5/10
    theta_v = 6*(16*12*cp*dp**11/100+sp.Rational(3, 5))
    omega_v = omega*(5+12*cp*dp**11/10)
    norm = theta+2*t*theta_v+t*theta**2
    sw = (4+t/100)*dp**6/10
    cf_delta_over_t = (cp**2*dp**8/4-1)/t
    cross_over_t = 4*sw*cf_delta_over_t
    ee_over_t = cross_over_t+2*omega*theta
    aa = 4*(cp*dp**12+2*cp**2*dp**8)/10+norm
    cc = t*cross_over_t+2*(omega+2*t*omega_v)+2*t*omega*theta
    ba = 4*(8+cp**3*dp**20/4)/10+norm+4*t*omega**2
    return {"r": r, "t": t, "c_upper": cp, "d_upper": dp,
            "theta_over_u": theta, "omega_over_u": omega,
            "theta_v": theta_v, "omega_v": omega_v,
            "A": aa, "C": cc, "E_over_r_squared": ee_over_t,
            "dc_and_fc": 2*omega, "b_analytic": ba,
            "ag_squared": 4*cp*dp**6/10, "bg_squared": sp.Rational(32, 10),
            "jL_squared": cp*dp**18/10, "jH_squared": 8*dp**6/10,
            "mass_N_delta_derivative": sp.Integer(32),
            "mass_N_v_derivative": sp.Integer(1008),
            "mass_remainder": 126+80*(sp.Rational(41, 320)+sp.Rational(121, 640)),
            "total_bounded_B": sp.Integer(161)}


@cache
def checks():
    c = envelope()
    ceilings = {"A": 9, "C": 12, "E_over_r_squared": 60,
                "dc_and_fc": 10, "b_analytic": 9,
                "ag_squared": 1, "bg_squared": 4,
                "jL_squared": 1, "jH_squared": 1, "mass_remainder": 152}
    margins = {name: sp.factor(upper-c[name]) for name, upper in ceilings.items()}
    margins["theta_negative_on_box"] = 8-c["c_upper"]*c["d_upper"]**12
    margins["D_delta_v_coefficient"] = sp.Rational(41, 10)-(4+sp.Rational(6, 100)+sp.Rational(4, 10000)+sp.Rational(1, 1000000))
    margins["D_v_squared_coefficient"] = sp.Rational(121, 10)-(12+sp.Rational(8, 100)+sp.Rational(2, 10000))
    v, delta = sp.symbols("v delta", nonnegative=True)
    den, base = (2+delta)*(1+v)**4-2, delta+8*v
    nn = 16*(1-v)*(8/((2+delta)*(1+v)**14)+1/(1+v)**2)
    d = operator.derive()
    residuals = {
        "D_positive_polynomial": sp.expand(den-base-delta*v*(4+6*v+4*v**2+v**3)-v**2*(12+8*v+2*v**2)),
        "mixed_denominator_bound": sp.expand(base**2-32*delta*v-(delta-8*v)**2),
        "actual_mass_v_formula": sp.factor(d["mass"].subs({d["u"]: sp.sqrt(v), d["delta"]: delta})-nn/den),
        "N_delta_derivative": sp.factor(sp.diff(nn, delta)+128*(1-v)/((2+delta)**2*(1+v)**14)),
        "N_v_derivative": sp.factor(sp.diff(nn, v)+16*(8/(2+delta)*((1+v)**-14+14*(1-v)*(1+v)**-15)
                                                              +(1+v)**-2+2*(1-v)*(1+v)**-3)),
        "mass_quotient_nonpositive": sp.factor(nn/den-80/base-(nn-80)/den+80*(den-base)/(den*base)),
    }
    if any(value.is_positive is not True for value in margins.values()):
        raise ValueError("A continuous actual-source coefficient margin failed")
    if any(value != 0 for value in residuals.values()):
        raise ValueError("A literal mass/pole identity failed")
    return {"residuals": residuals, "strict_margins": margins}
