"""Continuous analytic coefficient-l1 bounds, not a real-axis sample grid.

For every complex |K|<=4 use sum |f[j,n](K)| v^j R^n. All displayed
rational bounds control this Banach norm uniformly; v=R^2=1/400.
"""

from functools import cache

import sympy as sp

from .exact import DELTA_RADIUS, RADIUS


@cache
def calibration():
    r, v = RADIUS, DELTA_RADIUS
    cp, cm, dp, dm = 2+v, 2-v, 1+v, 1-v
    perturb = cp*dp**12-2
    den = 10-perturb
    derivative = 12*cp*dp**11
    theta = 6*(6+perturb)/(dm*den)
    theta_v = 6*(derivative/(dm*den)+(6+perturb)/(dm**2*den)
                 +(6+perturb)*derivative/(dm*den**2))
    omega = 48*dp**5/((1-v/2)*den)
    omega_v = 48/(1-v/2)*(5*dp**4/den+dp**5*derivative/den**2)
    normalization = theta+2*v*theta_v+v*theta**2
    omega_prime = omega+2*v*omega_v
    w1, w2 = cp*dp**12/den, 8/den
    cf = cp**2*dp**8/4
    sqrtw = 4*dp**6/((1-v/2)*den)
    q2 = 4/dm**4
    cross = sqrtw*(cf-1)
    a_majorant = q2*(w1+cf*w2)+normalization
    c_majorant = q2*cross+2*omega_prime+2*v*omega*theta
    e_majorant = q2*cross+2*v*omega*theta
    b_majorant = q2*(w2+cf*w1)+normalization+4*v*omega**2
    d_norm = cp*dp**4-2
    u_dprime = 8*v*cp*dp**3
    dpp_error = 8*cp*dp**3+48*v*cp*dp**2-16
    n_error = 16*((1+v)*(8/(cm*dm**14)+1/dm**2)-5)
    q_self = ((d_norm-9*v)/(7*v)+(16*cp*dp**3-32)/70
              +(dpp_error+n_error+9*d_norm)/84)
    light_self = 9*v/2
    light_from_q_scaled = (6+sp.Rational(5, 3))*d_norm+5*u_dprime
    q_from_light = (1+5*v)/7
    column_light = light_self+v*q_from_light
    column_q = light_from_q_scaled+q_self
    even_l_error = 2*(light_self+v*5*v/7)
    even_q = (5*v/7+q_from_light*even_l_error)/(1-q_self)
    odd_l_error = 2*column_light*r
    odd_q = q_from_light*(r+odd_l_error)/(1-q_self)
    out = {
        "R": r, "v": v, "denominator_lower": den,
        "theta_over_u": theta, "theta_v": theta_v,
        "omega_over_u": omega, "omega_v": omega_v,
        "normalization": normalization, "omega_prime": omega_prime,
        "A": a_majorant, "C": c_majorant, "E": e_majorant, "b_analytic": b_majorant,
        "D": d_norm, "u_Dprime": u_dprime, "Dpp_minus16": dpp_error,
        "N_minus80": n_error, "q_self": q_self,
        "light_self": light_self, "light_from_q_scaled": light_from_q_scaled,
        "q_from_light": q_from_light, "column_light": column_light, "column_q": column_q,
        "even_l_error": even_l_error, "even_q": even_q,
        "odd_l_error": odd_l_error, "odd_q": odd_q,
        "inverse_norm": sp.Rational(1, 84), "inverse_u_derivative_norm": sp.Rational(1, 70),
        "inverse_u2_second_derivative_norm": sp.Rational(1, 7),
        "inverse_spring": d_norm*(cp*dp**12+8)/(dm**6*(80-n_error)),
        "K_g": dp**6/8, "K_f": 1/(cm*dm**6), "G_g": dp**2/8,
    }
    # sqrt(c)<=3/(2(1-v/2)); sqrt(10)>3. The leading factor is 2.
    out["physical_ag"] = dp**3/((1-v/2)*(1-perturb/10))
    out["physical_bg"] = 2/(dm**3*(1-perturb/10))
    out["physical_bf"] = cp*dp**9/(3*(1-perturb/10))
    return out


def checks():
    c = calibration()
    ceilings = {"theta_over_u": 4, "omega_over_u": 5, "normalization": 4,
                "A": 9, "C": 12, "E": 60*c["v"], "b_analytic": 9,
                "column_light": sp.Rational(1, 2), "column_q": sp.Rational(1, 2),
                "even_l_error": 10*c["v"], "even_q": 3*c["v"],
                "odd_l_error": c["R"]/40, "odd_q": c["R"]/5,
                "inverse_spring": sp.Rational(1, 100),
                "physical_ag": 2, "physical_bg": 4, "physical_bf": 1,
                "N_minus80": 3, "K_g": 1, "K_f": 1, "G_g": 1}
    margins = {name: sp.factor(upper-c[name]) for name, upper in ceilings.items()}
    if any(value.is_positive is not True for value in margins.values()):
        raise ValueError("A continuous coefficient-l1 or contraction bound failed")
    return margins
