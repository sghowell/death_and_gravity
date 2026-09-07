"""Independent four-coordinate action and Arb whole-box derivations.

Neither calculation imports the primary stationary or Fraction engine.
Their results are compared separately in bridges.py.  The Arb calculation
evaluates the entire box, not sampled points or a fitted numerical root.
"""
from fractions import Fraction
from functools import cache

import sympy as sp
from flint import arb, ctx

T = sp.Symbol("audit_time", real=True)
M2 = sp.Symbol("audit_M_squared", positive=True)


def curvature(metric):
    """Literal time-dependent diagonal metric, all four coordinate indices."""
    inverse = metric.inv()

    def partial(value, index):
        return sp.diff(value, T) if index == 0 else sp.S.Zero

    gamma = [[[sp.simplify(sum(inverse[a, d]*(partial(metric[d, c], b)
                +partial(metric[d, b], c)-partial(metric[b, c], d))
                for d in range(4))/2) for c in range(4)] for b in range(4)]
                for a in range(4)]
    ricci = sp.zeros(4)
    for a in range(4):
        for b in range(4):
            ricci[a, b] = sp.simplify(sum(partial(gamma[c][a][b], c)
                -partial(gamma[c][a][c], b)
                +sum(gamma[c][c][d]*gamma[d][a][b]-gamma[c][b][d]*gamma[d][a][c]
                     for d in range(4)) for c in range(4)))
    scalar = sp.simplify(sum(inverse[a, b]*ricci[a, b]
                             for a in range(4) for b in range(4)))
    return {"Gamma": gamma, "Ricci": ricci, "R": scalar}


@cache
def action():
    b, n, q, h = (sp.Function(name)(T) for name in ("audit_b", "audit_N", "audit_q", "audit_h"))
    p, w = (sp.Function(name)(T) for name in ("audit_beta1", "audit_beta4"))
    g = sp.diag(1, -sp.exp(q), -sp.exp(-q), -1)
    f = sp.diag(n**2, -b**2*sp.exp(h), -b**2*sp.exp(-h), -b**2)
    gc, fc = curvature(g), curvature(f)
    # Positive lapse and scale fix the signed volume and square-root branch.
    roots = (n, b*sp.exp((h-q)/2), b*sp.exp(-(h-q)/2), b)
    potential = -2*p*sum(roots)-2*w*sp.prod(roots)
    boundary = 3*M2*b**2*sp.diff(b, T)/n
    g_eh = -M2*gc["R"]/2
    f_eh = -M2*n*b**3*fc["R"]/2
    f_adm = sp.simplify(f_eh-sp.diff(boundary, T))
    amplitude = sp.Symbol("audit_amplitude", real=True)
    density = g_eh+f_adm+potential
    deformed = density.subs({q: amplitude*q, h: amplitude*h}).doit()
    quadratic = sp.simplify(sp.diff(deformed, amplitude, 2).subs(amplitude, 0)/2)
    isotropic = sp.simplify(density.subs({q: 0, h: 0}).doit())
    constraint = sp.diff(isotropic, n)
    spatial = sp.diff(isotropic, b)-sp.diff(sp.diff(isotropic, sp.diff(b, T)), T)
    euler = {"g": sp.simplify(-2*(sp.diff(quadratic, q)
                     -sp.diff(sp.diff(quadratic, sp.diff(q, T)), T))),
             "f": sp.simplify(-2*(sp.diff(quadratic, h)
                     -sp.diff(sp.diff(quadratic, sp.diff(h, T)), T)))}
    return {"b": b, "N": n, "q": q, "h": h, "p": p, "w": w,
            "g_R": gc["R"], "f_R": fc["R"], "g_EH": g_eh, "f_EH": f_eh,
            "f_ADM": f_adm, "boundary": boundary, "potential": potential,
            "quadratic": quadratic, "isotropic": isotropic,
            "C": constraint, "E": spatial, "tensor_equations": euler}


def arb_box(precision=256):
    """Independent rigorous ball enclosure with separately written derivatives.

    Endpoints and final comparison margins are rational.  Binary floating
    point is never used to construct an input or decide an inequality.
    """
    if isinstance(precision, bool) or not isinstance(precision, int):
        raise TypeError("Arb precision must be an integer")
    if precision < 128:
        raise ValueError("At least 128 bits are required for this audit")
    with ctx.workprec(precision):
        def rational(value):
            value = Fraction(value)
            return arb(value.numerator)/value.denominator

        def interval(lower, upper):
            lo, hi = Fraction(lower), Fraction(upper)
            return arb(rational((lo+hi)/2), rational((hi-lo)/2))

        v, delta, z = interval(0, "1/10000"), interval(0, "1/100"), interval(11, 13)
        c, d = 2+delta, 1+v
        gap = delta+2*(1-d**-4)
        gap_v, gap_vv = 8/d**5, -40/d**6
        j = c*d**4*(1-7*v)+12*v
        jv = c*(4*d**3*(1-7*v)-7*d**4)+12
        jvv = c*(12*d**2*(1-7*v)-56*d**3)
        numerator = 32*(1-v)/(c*d**14)
        ratio = d**8*j/(8*c*(1-v))
        ell = 8/d+jv/j+1/(1-v)
        ell_v = -8/d**2+jvv/j-(jv/j)**2+1/(1-v)**2
        hh = gap*(jv/j-6/d)-gap_v
        hh_v = gap_v*(jv/j-6/d)+gap*(jvv/j-(jv/j)**2+6/d**2)-gap_vv
        one_minus = 1-v*gap*z
        b = (one_minus/ratio).root(3)
        velocity = 2*(v*z*hh-ell)/3
        velocity_v = 2*(z*hh+v*z*hh_v-ell_v)/3
        velocity_z = 2*v*hh/3
        bv_over_b = (-(gap+v*gap_v)*z/one_minus-ell)/3
        bz_over_b = -v*gap/(3*one_minus)
        fixed = 3*b*velocity**2/(2*numerator)
        fixed_v = fixed*(bv_over_b+2*velocity_v/velocity+1/(1-v)+14/d)
        fixed_z = fixed*(bz_over_b+2*velocity_z/velocity)
        z_v = fixed_v/(1-fixed_z)
        b_v = b*(bv_over_b+bz_over_b*z_v)
        lapse = 2*b_v/velocity
        enclosures = {"Q": (numerator, 15, 17), "R": (ratio, "3/25", "13/100"),
                      "b": (b, "19/10", "21/10"), "F": (velocity, "-81/10", "-79/10"),
                      "T": (fixed, 11, 13), "T_zeta": (fixed_z, "-1/100", "1/100"),
                      "N": (lapse, "19/10", "21/10"), "b_v": (b_v, "-83/10", "-77/10")}
        margins = {name: {"lower": str(lo), "upper": str(hi),
                          "strict": bool(value > rational(lo) and value < rational(hi))}
                   for name, (value, lo, hi) in enclosures.items()}
        if not all(item["strict"] for item in margins.values()):
            raise ValueError("An independent Arb whole-box strict enclosure failed")
        return {"precision_bits": precision, "whole_box_not_grid": True,
                "box": {"v": ["0", "1/10000"], "delta": ["0", "1/100"], "zeta": ["11", "13"]},
                "strict_rational_enclosures": margins,
                "endpoint_is_only_analytic_extension": True}
