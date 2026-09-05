"""Tensor pullbacks and finite-slice tests of NEC Einstein-scalar parents.

Global concavity is proved in the accompanying note, not by finite sampling.
The finite-slice inequalities are sufficient obstruction tests under that lemma.
"""

from functools import cache

import sympy as sp

from .conformal import exact_rational


def tensor_map(C, W):
    """Normalized physical (G_T,F_T) for an Einstein seed; C,W must be positive."""
    C, W = map(sp.sympify, (C, W))
    return C**sp.Rational(3, 2)/sp.sqrt(W), sp.sqrt(C*W)


def three_slice_bound(L, epsilon):
    """Uniform tensor-value errors at t=(-L*tau,0,L*tau), exact target geometry."""
    L, epsilon = map(exact_rational, (L, epsilon))
    if not (L > 0 and 0 <= epsilon < 1):
        raise ValueError("Require L>0 and 0<=epsilon<1")
    growth2 = (1+L**2)**2
    threshold = (growth2-1)/(growth2+1)
    margin = sp.factor(growth2*(1-epsilon)-(1+epsilon))
    return {"L": str(L), "epsilon": str(epsilon), "squared_scale_margin": str(margin),
            "necessary_error_lower_bound": str(threshold),
            "status": "EXCLUDED_BY_EINSTEIN_CONCAVITY" if margin > 0 else "BOUND_INCONCLUSIVE",
            "warning": "Necessary only; equality does not construct a parent"}


@cache
def identities():
    C, W, a, M = sp.symbols("C W a M", positive=True)
    velocity, gradient = sp.symbols("tensor_velocity tensor_gradient", real=True)
    # Pull back the complete two-derivative tensor principal action, including
    # proper-time measure, physical volume and each derivative's scale factor.
    aE, nE = sp.sqrt(C)*a, sp.sqrt(W)
    lag = M**2*nE*aE**3*(velocity**2/nE**2-gradient**2/aE**2)/8
    g = sp.simplify(4*sp.diff(lag, velocity, 2)/(M**2*a**3))
    f = sp.simplify(-4*sp.diff(lag, gradient, 2)/(M**2*a))
    expected_g, expected_f = tensor_map(C, W)
    gt, ft = sp.symbols("g_T f_T", positive=True)
    inverse_C, inverse_W = sp.sqrt(gt*ft), ft**sp.Rational(3, 2)/sp.sqrt(gt)
    round_g, round_f = tensor_map(inverse_C, inverse_W)
    t = sp.Symbol("t", real=True)
    scale, spatial, lapse = (sp.Function(name)(t) for name in ("a", "C", "n"))
    H = sp.diff(scale, t)/scale
    a_aux = sp.sqrt(spatial)*scale
    HE = sp.simplify(sp.diff(a_aux, t)/(a_aux*lapse))
    HEd = sp.simplify(sp.diff(HE, t)/lapse)
    HE_expected = (H+sp.diff(spatial, t)/(2*spatial))/lapse
    HEd_expected = (sp.diff(H, t)+sp.diff(spatial, t, 2)/(2*spatial)
                    -sp.diff(spatial, t)**2/(2*spatial**2)
                    -(H+sp.diff(spatial, t)/(2*spatial))*sp.diff(lapse, t)/lapse)/lapse**2
    # Jensen/secant identity for arbitrary positive proper-time intervals.
    left, right = sp.symbols("left_duration right_duration", positive=True)
    ym, y0, yp = sp.symbols("log_a_minus log_a_zero log_a_plus", real=True)
    weighted_excess = right*ym+left*yp-(left+right)*y0
    slopes = (y0-ym)/left-(yp-y0)/right
    L, epsilon = sp.symbols("L epsilon", positive=True)
    threshold = ((1+L**2)**2-1)/((1+L**2)**2+1)
    margin = (1+L**2)**2*(1-epsilon)-(1+epsilon)
    residuals = {"tensor_kinetic_pullback": sp.simplify(g-expected_g),
                 "tensor_gradient_pullback": sp.simplify(f-expected_f),
                 "tensor_product": sp.simplify(g*f-C**2),
                 "inverse_G": sp.simplify(round_g-gt), "inverse_F": sp.simplify(round_f-ft),
                 "H_E": sp.simplify(HE-HE_expected), "H_E_dot": sp.simplify(HEd-HEd_expected),
                 "concavity_secants": sp.expand(left*right*slopes+weighted_excess),
                 "tensor_error_threshold": sp.factor(margin.subs(epsilon, threshold))}
    return {"normalized_G_T": g, "normalized_F_T": f,
            "inverse_C": inverse_C, "inverse_W": inverse_W,
            "three_slice_error_threshold": sp.factor(threshold),
            "all_time_threshold_limit": sp.limit(threshold, L, sp.oo),
            "residuals": residuals}


def derive():
    result = identities()
    if any(sp.simplify(value) != 0 for value in result["residuals"].values()):
        raise ValueError("Einstein tensor/frame identity failed")
    return {**{key: str(value) for key, value in result.items() if key != "residuals"},
            "exact_residuals": {key: "0" for key in result["residuals"]},
            "three_slice_examples": {"ten_percent_L1": three_slice_bound(1, sp.Rational(1, 10)),
                                     "L1_threshold_not_sufficient": three_slice_bound(1, sp.Rational(3, 5)),
                                     "L2_threshold_not_sufficient": three_slice_bound(2, sp.Rational(12, 13))},
            "global_result": "No regular disformal pullback of a classical flat Einstein/NEC-scalar solution can have physical a growing in both tails and F_T,G_T tending to positive constants",
            "scope": "Exact Einstein seed with positive scalar field metric (or separately proved NEC); excludes no general Horndeski/DHOST/quantum parent"}
