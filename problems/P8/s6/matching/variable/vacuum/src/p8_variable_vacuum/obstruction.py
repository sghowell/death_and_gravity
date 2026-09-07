"""Exact local profile, independent vacuum ratio, and residual identities."""

from fractions import Fraction
from functools import cache

import sympy as sp
from p8_variable_beta import background


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("Require an exact rational input, not bool or floating point")
    return sp.Rational(value)


def point(time, action_lapse):
    u, c = exact(time), exact(action_lapse)
    if abs(u) > sp.Rational(1, 10):
        raise ValueError("Only the specified local clock interval |u|<=1/10 is in scope")
    if not (c == 1 or 2 < c <= 4):
        raise ValueError("Only the action parameters c=1 or 2<c<=4 are in scope")
    return u, c


@cache
def derive():
    d = background.derive()
    u, c, h, hp, y = (d[key] for key in ("u", "c", "h", "h_u", "y"))
    r = sp.Symbol("r", positive=True)
    b0, b1, b4 = (d[key] for key in ("b0", "b1", "b4"))
    theta = sp.factor(3*h**2*(c-y)/(2*c*hp))
    root_cube = sp.factor(-b1/b4)
    U = b0+3*b1*r
    V = b1+b4*r**3
    clock = sp.diff(b0, u)+4*r*sp.diff(b1, u)+r**4*sp.diff(b4, u)
    return {"background": d, "u": u, "c": c, "r": r, "theta": theta,
            "r_star_cubed": root_cube, "g_flat_residual": U, "f_flat_residual": V,
            "clock_profile_residual": clock,
            "b0_u": sp.diff(b0, u), "b1_u": sp.diff(b1, u), "b4_u": sp.diff(b4, u)}


def checks():
    d = derive()
    old, u, c, r = d["background"], d["u"], d["c"], d["r"]
    h2, hp, y, cc, rr = sp.symbols("h_squared h_prime y c r", positive=True)
    n = 2*(y**3/cc-1)*hp
    b1 = y**3*hp/(cc*(cc-y))
    b4 = 3*h2/(2*cc**2)-b1/y**3
    b0 = 3*h2/2-n/4-3*b1*y
    U, V, S = b0+3*b1*rr, b1+b4*rr**3, rr**2+rr*y+y**2
    rhs = 3*h2/2-n/4+9*h2*y**3*rr**3/(2*cc**2*S)
    M, tau, speed = sp.symbols("M tau kbar", positive=True)
    slope = sp.Symbol("b_u", real=True)
    phi = sp.Symbol("phi", real=True)
    inverse = sp.Function("u_of_phi")(phi)
    profile = sp.Function("b")(inverse)
    actual_beta = M**2*profile/tau**2
    differentiated_beta = sp.diff(actual_beta, phi).subs(
        {sp.diff(inverse, phi): 1/(M*sp.sqrt(speed)), sp.diff(profile, inverse): slope})
    values = {
        "root_cube_from_f_metric_equation": sp.factor(old["b1"]+old["b4"]*d["r_star_cubed"]),
        "root_cube_and_theta_dictionary": sp.factor(d["r_star_cubed"]-old["y"]**3/(1-d["theta"])),
        "full_metric_residual_combination": sp.factor(U+3*y**3*V/S-rhs),
        "positive_root_sum_square_bound": sp.expand(S-3*rr*y-(rr-y)**2),
        "clock_dimensionful_chain_rule": sp.factor(differentiated_beta-M*slope/(tau**2*sp.sqrt(speed))),
        "center_f_flat_root": sp.factor(d["f_flat_residual"].subs({u: 0, r: 2})),
        "center_g_residual_at_f_root": sp.factor(d["g_flat_residual"].subs({u: 0, r: 2})-(2-16/c)),
        "center_clock_equation_at_arbitrary_ratio": sp.factor(d["clock_profile_residual"].subs(u, 0)),
    }
    values.update({f"even_profile_center_derivative_b{j}": sp.factor(sp.diff(old[f"b{j}"], u).subs(u, 0))
                   for j in (0, 1, 4)})
    return values


def evaluate(time, action_lapse, vacuum_ratio=2):
    """Exact coefficient-label evaluation, not a claim of a vacuum solution.

    time is u(phi/M), not the time coordinate of the stationary vacuum.
    The canonical positive clock makes this a bijective label of the field interval.
    """
    u, c = point(time, action_lapse)
    r = exact(vacuum_ratio)
    if r <= 0:
        raise ValueError("The independent proportional vacuum root must be strictly positive")
    d = derive()
    at = {d["u"]: u, d["c"]: c, d["r"]: r}
    names = ("a", "y", "h", "h_u", "nbar", "kbar", "b0", "b1", "b4")
    out = {name: sp.cancel(d["background"][name].subs(at)) for name in names}
    names = ("theta", "r_star_cubed", "g_flat_residual", "f_flat_residual",
             "clock_profile_residual", "b0_u", "b1_u", "b4_u")
    out.update({name: sp.cancel(d[name].subs(at)) for name in names})
    out["metric_eliminant"] = out["b0"]**3*out["b4"]-27*out["b1"]**4
    return out


def controls():
    at_f = evaluate(0, 4, 2)
    at_g = evaluate(0, 4, sp.Rational(13, 6))
    if (at_f["f_flat_residual"], at_f["clock_profile_residual"], at_f["g_flat_residual"]) != (0, 0, -2):
        raise ValueError("Omitting the g metric equation was not detected")
    if (at_g["g_flat_residual"], at_g["clock_profile_residual"], at_g["f_flat_residual"]) != (0, 0, -sp.Rational(469, 432)):
        raise ValueError("Omitting the f metric equation was not detected")
    changed = {"g_flat_residual": at_f["g_flat_residual"]+2,
               "f_flat_residual": at_f["f_flat_residual"],
               "clock_profile_residual": at_f["clock_profile_residual"]}
    if set(changed.values()) != {0}:
        raise ValueError("The explicitly changed-action cosmological-offset control failed")
    d = derive()["background"]
    changed_bounce = sp.factor(3*d["h"]**2-d["nbar"]/2-2*(d["b0"]+2+3*d["b1"]*d["y"]))
    if changed_bounce != -4:
        raise ValueError("The changed action did not expose its failure on the old bounce")
    return {"actual_center_f_and_clock_only": {key: at_f[key] for key in
               ("g_flat_residual", "f_flat_residual", "clock_profile_residual")},
            "actual_center_g_and_clock_only": {key: at_g[key] for key in
               ("g_flat_residual", "f_flat_residual", "clock_profile_residual")},
            "separately_changed_b0_plus2_at_c4_center": changed,
            "changed_action_old_bounce_g_lapse_residual": changed_bounce,
            "changed_action_control_scope": "Adding a constant2 to b0 changes the action on the bounce interval; it is not the frozen candidate",
            "singular_c2": "Not an action in the admitted family; cancelling divergent coefficients alone is not a vacuum"}
