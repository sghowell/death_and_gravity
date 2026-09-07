"""Conditional fixed-subset null combination and all-point comparison data.

At a time T take the component through edges with P(T)!=0. Only that
fixed vertex subset is differentiated. Boundary P=0 edges have zero null
stress at T; they need not stay zero or have a differentiable branch label.
"""

from functools import cache

import sympy as sp
from p8_star.exact import nonnegative, positive, rational


@cache
def derive(size=3):
    if type(size) is not int or size < 1:
        raise ValueError("A nonempty finite component is required")
    h, hp = sp.symbols("H_r H_r_prime", real=True)
    gs = (sp.Symbol("G_r", positive=True),)+sp.symbols(f"G1:{size}", nonnegative=True)
    ys = (sp.Integer(1),)+sp.symbols(f"y1:{size}", positive=True)
    cs = (sp.Integer(1),)+sp.symbols(f"c1:{size}", positive=True)
    nulls = sp.symbols(f"n0:{size}", nonnegative=True)
    yps = tuple((c-y)*h for c, y in zip(cs, ys, strict=True))
    a = sum(g*y**2 for g, y in zip(gs, ys, strict=True))
    ap = 2*sum(g*y*yp for g, y, yp in zip(gs, ys, yps, strict=True))
    weights = tuple(c*y**3 for c, y in zip(cs, ys, strict=True))
    source = sum(w*n for w, n in zip(weights, nulls, strict=True))
    individual = tuple((hp/y-h*yp/y**2)/c for c, y, yp in zip(cs, ys, yps, strict=True))
    full = -2*sum(g*w*dh for g, w, dh in zip(gs, weights, individual, strict=True))
    lam, eta = sp.factor(ap/(2*a)), sp.factor(source/(2*a))
    return {"size": size, "G": gs, "y": ys, "c": cs, "nulls": nulls,
            "H": h, "Hprime": hp, "yprimes": yps, "weights": weights,
            "A": a, "A_fixed_prime": sp.factor(ap), "weighted_source": source,
            "full_weighted_null_left": full, "residual": sp.factor(-2*a*hp+h*ap-source),
            "lambda": lam, "eta": eta, "reconstructed_Hprime": sp.factor(lam*h-eta)}


@cache
def checks(size=3):
    d = derive(size)
    a, ap, h, hp = (d[key] for key in ("A", "A_fixed_prime", "H", "Hprime"))
    result = {"fixed_component_weighted_null": sp.expand(d["full_weighted_null_left"]+2*a*hp-h*ap),
              "component_source_reconstruction": sp.factor(d["residual"].subs(hp, d["reconstructed_Hprime"])),
              "bounded_coefficient_weighted_average": sp.factor(d["lambda"]-sum(g*y**2*(yp/y) for g, y, yp in zip(d["G"], d["y"], d["yprimes"], strict=True))/a),
              "target_kinetic_lower_bound": sp.expand(a-d["G"][0]-sum(g*y**2 for g, y in zip(d["G"][1:], d["y"][1:], strict=True))),
              "fixed_stratum_normalized_identity": sp.simplify((hp/sp.sqrt(a)-h*ap/(2*a**sp.Rational(3, 2)))+(d["residual"]+d["weighted_source"])/(2*a**sp.Rational(3, 2)))}
    for index, (c, y, yp) in enumerate(zip(d["c"], d["y"], d["yprimes"], strict=True)):
        result[f"internal_clock_ratio_kinematics_{index}"] = sp.expand(yp-y*(c*h/y-h))
    return result


def reconstruct(Gs, ys, cs, nulls, H, component, root=0):
    """Exact conditional dynamic-component data, not a full solution maker.

    G>=0, target G_root>0; arrays are full finite-vertex data. On this
    supplied component, Hi=H/yi and yi'=(ci-yi)H are imposed. Being the
    actual nonzero-P component and satisfying all vertex equations are
    separate premises. No derivative of component membership is taken.
    """
    if not all(isinstance(v, (tuple, list)) for v in (Gs, ys, cs, nulls)):
        raise TypeError("Explicit finite coefficient lists are required")
    count = len(Gs)
    if count < 1 or any(len(v) != count for v in (ys, cs, nulls)):
        raise ValueError("All vertex arrays must have the same nonzero length")
    if type(root) is not int or not 0 <= root < count:
        raise ValueError("A valid physical root is required")
    if not isinstance(component, (tuple, list)) or not component or any(type(i) is not int or not 0 <= i < count for i in component):
        raise ValueError("An explicit nonempty vertex subset is required")
    if root not in component or len(set(component)) != len(component):
        raise ValueError("The component contains the physical root exactly once")
    gs = tuple(nonnegative(g, "Einstein coefficient") for g in Gs)
    if gs[root] <= 0:
        raise ValueError("The distinguished physical Einstein coefficient must be positive")
    yy = tuple(positive(y, "physical scale ratio") for y in ys)
    cc = tuple(positive(c, "physical lapse ratio") for c in cs)
    nn = tuple(nonnegative(n, "separate NEC null density") for n in nulls)
    h = rational(H)
    if yy[root] != 1 or cc[root] != 1:
        raise ValueError("The distinguished root's own ratios must equal one")
    subset = tuple(sorted(component))
    yp = {i: (cc[i]-yy[i])*h for i in subset}
    a = sum(gs[i]*yy[i]**2 for i in subset)
    ap = 2*sum(gs[i]*yy[i]*yp[i] for i in subset)
    source = sum(cc[i]*yy[i]**3*nn[i] for i in subset)
    lam, eta = ap/(2*a), source/(2*a)
    bound = max(abs(yp[i]/yy[i]) for i in subset)
    if not (a >= gs[root] > 0 and eta >= 0 and abs(lam) <= bound):
        raise ValueError("Positive component denominator or logarithmic-rate bound failed")
    return {"component": subset, "A": a, "A_fixed_prime": ap, "weighted_source": source,
            "lambda": lam, "eta": eta, "Hprime": lam*h-eta,
            "component_log_rate_bound": bound, "yprimes": yp}


def cd_endpoint_test(length, normalized_error):
    length, error = positive(length, "L"), nonnegative(normalized_error, "normalized endpoint error")
    threshold = 4*length/(1+length**2)
    return {"necessary_normalized_endpoint_error": threshold,
            "status": "EXCLUDED_WITH_ACTUAL_VERTEX_METRIC_AND_PROPER_CLOCK" if error < threshold else "INCONCLUSIVE"}
