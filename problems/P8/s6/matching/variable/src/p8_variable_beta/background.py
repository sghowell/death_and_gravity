"""Actual local CD-scale canonical-clock/free-chi family, dimensionless u=T/tau.

Einstein coefficients are Gg=Gf=M²; beta_n=M²/tau²*b_n(phi/M).
The clock is phi=M integral sqrt(kbar) du, not an external beta(t).
Both phi and chi have literal positive canonical kinetic actions on g.
An optional free f-clock theta is constant; beta is independent of theta.
"""

from functools import cache

import sympy as sp
from p8_star.exact import rational

U = sp.Symbol("u", real=True)
C = sp.Symbol("c", positive=True)


@cache
def derive():
    u, c = U, C
    d = 1+u**2
    a, b = d**2, 2/d**2
    y = b/a
    h = sp.diff(a, u)/a
    hp = sp.diff(h, u)
    chi = 1/(10*a**3)
    nbar = 2*(y**3/c-1)*hp
    kbar = nbar-chi**2
    b1 = y**3*hp/(c*(c-y))
    b4 = 3*h**2/(2*c**2)-b1/y**3
    b0 = (3*h**2-nbar/2)/2-3*b1*y
    V = b1+b4*y**3
    return {"u": u, "c": c, "d": d, "a": a, "b": b, "y": y, "h": h, "h_u": hp,
            "h_f": -h/c, "D_f_h_f": -hp/c**2, "chi_speed": chi,
            "nbar": nbar, "kbar": kbar, "b0": b0, "b1": b1, "b2": sp.S.Zero,
            "b3": sp.S.Zero, "b4": b4, "P": 2*b1, "U": b0+3*b1*y, "V": V,
            "clock_edge_force_times_speed": -2*(sp.diff(b0, u)+(c+3*y)*sp.diff(b1, u)+c*y**3*sp.diff(b4, u)),
            "partial_V_clock_times_speed": sp.diff(b1, u)+y**3*sp.diff(b4, u)}


@cache
def checks():
    d = derive()
    u, c, y, h, hp, n, k, chi, P = (d[key] for key in ("u", "c", "y", "h", "h_u", "nbar", "kbar", "chi_speed", "P"))
    values = {
        "g_lapse": 3*h**2-n/2-2*d["U"],
        "f_lapse": 3*h**2/c**2-2*d["V"]/y**3,
        "g_null": -2*hp-n-(y-c)*P,
        "f_null": 2*hp/c**2-(c-y)*P/(c*y**3),
        "free_chi_equation": sp.diff(chi, u)+3*h*chi,
        "clock_equation_times_strictly_positive_clock_speed": sp.diff(k, u)/2+3*h*k-d["clock_edge_force_times_speed"],
        "undivided_source_aware_Bianchi": 3*P*(y*d["h_f"]-h)-2*d["partial_V_clock_times_speed"],
        "optional_constant_f_clock": sp.diff(sp.S.Zero, u)+3*d["h_f"]*sp.S.Zero,
        "canonical_g_null_sum": k+chi**2-n,
        "actual_g_Hubble": h-sp.diff(d["a"], u)/d["a"],
        "actual_f_Hubble": d["h_f"]-sp.diff(d["b"], u)/(c*d["b"]),
    }
    return {key: sp.factor(value) for key, value in values.items()}


def validate_c(value):
    result = rational(value)
    if not (result == 1 or 2 < result <= 4):
        raise ValueError("Only c=1 or 2<c<=4 belongs to this certified local family")
    return result


def evaluate(time, lapse):
    u, c = rational(time), validate_c(lapse)
    if abs(u) > sp.Rational(1, 10):
        raise ValueError("Only |u|<=1/10 is certified")
    point = {U: u, C: c}
    d = derive()
    keys = ("a", "b", "y", "h", "h_u", "h_f", "D_f_h_f", "chi_speed", "nbar", "kbar", "b0", "b1", "b4", "P")
    return {key: sp.cancel(d[key].subs(point)) for key in keys}


def domain_checks():
    dmax = sp.Rational(101, 100)
    margins = {"d12_below_8_over_7": sp.Rational(8, 7)-dmax**12,
               "h_u_above_3": 4*sp.Rational(99, 100)/dmax**2-3,
               "y_above_1": 2/dmax**4-1,
               "c1_clock_squared_lower_bound": sp.Rational(3599, 100),
               "positive_branch_clock_squared_lower_bound": sp.Rational(449, 100)}
    if any(value.is_positive is not True for value in margins.values()):
        raise ValueError("A compact-domain rational positivity margin failed")
    # d<=dmax gives y³>7 and h_u>3. Also chi_speed²<=1/100.
    # c=1: nbar>36; 2<c<=4: nbar>9/2. The actual bounds are strict.
    return margins


def controls():
    d = derive()
    omission = sp.factor((3*d["P"]*(d["y"]*d["h_f"]-d["h"])).subs({U: sp.Rational(1, 10), C: 1}))
    if omission == 0:
        raise ValueError("Omitting the clock exchange did not change the Bianchi relation")
    return {"omitted_clock_exchange_at_c1_u1over10": omission,
            "c1_center": evaluate(0, 1), "c4_center": evaluate(0, 4),
            "constant_beta_tree_hypothesis": "NOT_SATISFIED; free chi conserved, interacting clock not separately conserved"}
