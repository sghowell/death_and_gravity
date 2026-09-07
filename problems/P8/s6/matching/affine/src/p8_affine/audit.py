"""Independent source-normal-form and CD-domain audit.

No primary affine dictionary or connection solver is imported.  The
principal source formulas are checked by independent substitution and
literal tensor contractions.  Lower-order matching is conditional on
the separately audited correction of the printed Q2 bracket.
"""

from fractions import Fraction
from functools import cache

import sympy as sp

U = sp.Symbol("clock_u", real=True)
X = sp.Symbol("paper_x", negative=True)
Y = sp.Symbol("repository_X", positive=True)
H = sp.Symbol("h", positive=True)


def _zero(expression):
    return sp.factor(sp.cancel(sp.simplify(expression)))


@cache
def parameters():
    """Paper variables, with h independent of x during x differentiation."""
    w = H - 1 - X
    f = w / (2 * H)
    p = sp.sqrt(w / H) / 2
    c = 2 * (p - f) / X
    quartic = (sp.Rational(1, 2) - f) / X**2
    return {
        "u": U, "x": X, "h": H, "w": w, "f": f, "p": p,
        "c": c, "F4": quartic,
        "p_x": sp.diff(p, X), "c_x": sp.diff(c, X),
        "denominator": _zero(2 * p - c * X + 2 * quartic * X**2),
        "time_substitution": {H: (1 + U**2)**3},
    }


@cache
def principal():
    """Literal source Eqs.4.5,4.9--4.12, not the primary lift code."""
    d = parameters()
    p, c, quartic, px, cx, den = (
        d[name] for name in ("p", "c", "F4", "p_x", "c_x", "denominator")
    )
    a1 = -c / 2 - p * (c - 2 * quartic * X) / den
    a3 = 2 * cx + (
        4 * p * quartic + (4 * px - c) * (c - 2 * quartic * X)
    ) / den
    a4 = (
        -2 * cx + 2 * px * (3 * px - c) / p
        + px * X * (px * c - 4 * p * cx) / p**2
        + (c**2 - 4 * p * quartic - 2 * c * quartic * X) / den
    )
    a5 = (
        -px * (px * c - 4 * p * cx) / p**2
        + 2 * px * (4 * p * quartic + (3 * px - c) * (c - 2 * quartic * X))
        / (p * den)
    )
    raw = {"f": p - c * X / 2, "alpha1": a1, "alpha2": -a1,
           "alpha3": a3, "alpha4": a4, "alpha5": a5}
    return {"raw": raw, "simplified": {key: _zero(value) for key, value in raw.items()}}


@cache
def target():
    """Independently specialize the frozen A1=0 Ia formulas to CD.

    This is a source-normal-form calculation, not a fresh proof of Ia
    constraint rank.  No affine primary module or p8.matter is imported.
    """
    curvature = -sp.Rational(1, 2) + (1 - Y) / (2 * H)
    derivative = sp.diff(curvature, Y)
    a3 = 1 / (H * Y)
    a4 = (
        -Y**2 * curvature * a3**2
        + 8 * curvature * (Y * derivative - curvature) * a3
        + 48 * curvature * derivative**2
    ) / (8 * curvature**2)
    a5 = (4 * derivative + Y * a3) * (4 * curvature * a3) / (8 * curvature**2)
    repo = {"F2": curvature, "A1": sp.S.Zero, "A2": sp.S.Zero,
            "A3": a3, "A4": _zero(a4), "A5": _zero(a5)}
    paper = {
        "f": -curvature.subs(Y, -X),
        "alpha1": sp.S.Zero, "alpha2": sp.S.Zero,
        "alpha3": -a3.subs(Y, -X), "alpha4": -a4.subs(Y, -X),
        "alpha5": a5.subs(Y, -X),
    }
    return {"repository": repo, "paper": {key: _zero(value) for key, value in paper.items()}}


def principal_checks():
    actual, expected = principal()["simplified"], target()["paper"]
    return {f"target_{name}": _zero(actual[name] - expected[name])
            for name in ("f", "alpha1", "alpha3", "alpha4", "alpha5")}


def _exact(value):
    if isinstance(value, bool) or value is sp.true or value is sp.false:
        raise TypeError("exact rational data exclude booleans")
    if type(value) is int or isinstance(value, Fraction):
        return Fraction(value)
    if isinstance(value, sp.Rational):
        return Fraction(int(value.p), int(value.q))
    raise TypeError("expected exact int, Fraction, or SymPy Rational")


def _vector(values):
    values = tuple(values)
    if len(values) != 4:
        raise ValueError("four components are required")
    return tuple(_exact(value) for value in values)


def _matrix(values):
    values = tuple(tuple(row) for row in values)
    if len(values) != 4 or any(len(row) != 4 for row in values):
        raise ValueError("a four-by-four matrix is required")
    result = tuple(tuple(_exact(value) for value in row) for row in values)
    if any(result[i][j] != result[j][i] for i in range(4) for j in range(4)):
        raise ValueError("a symmetric covariant matrix is required")
    return result


def contractions(gradient, hessian, *, paper=False):
    """Independent Fraction contractions for a Lorentz frame and its sign flip."""
    if type(paper) is not bool:
        raise TypeError("paper must be an actual bool")
    v, tensor = _vector(gradient), _matrix(hessian)
    metric = (-1, 1, 1, 1) if paper else (1, -1, -1, -1)
    raised = tuple(metric[i] * v[i] for i in range(4))
    x = sum(metric[i] * v[i]**2 for i in range(4))
    trace = sum(metric[i] * tensor[i][i] for i in range(4))
    vHv = sum(raised[i] * tensor[i][j] * raised[j] for i in range(4) for j in range(4))
    hv = tuple(sum(tensor[i][j] * raised[j] for j in range(4)) for i in range(4))
    return {
        "X": x, "Box": trace, "vHv": vHv,
        "L1": sum(metric[i] * metric[j] * tensor[i][j]**2
                  for i in range(4) for j in range(4)),
        "L2": trace**2, "L3": trace * vHv,
        "L4": sum(metric[i] * hv[i]**2 for i in range(4)),
        "L5": vHv**2,
    }


def signature_checks():
    gradient = (2, Fraction(1, 3), -1, Fraction(2, 5))
    hessian = ((1, 2, -1, 3), (2, -2, 4, 0), (-1, 4, 5, -2), (3, 0, -2, 6))
    repo = contractions(gradient, hessian)
    paper = contractions(gradient, hessian, paper=True)
    signs = {"X": -1, "Box": -1, "vHv": 1, "L1": 1,
             "L2": 1, "L3": -1, "L4": -1, "L5": 1}
    ricci_diagonal = (Fraction(3, 2), -2, Fraction(5, 3), 7)
    scalar_repo = sum(g * r for g, r in zip((1, -1, -1, -1), ricci_diagonal, strict=True))
    scalar_paper = sum(g * r for g, r in zip((-1, 1, 1, 1), ricci_diagonal, strict=True))
    values = {f"signature_{name}": paper[name] - sign * repo[name]
              for name, sign in signs.items()}
    values["signature_R"] = scalar_repo + scalar_paper
    values["minimal_free_matter"] = -paper["X"] / 2 - repo["X"] / 2
    return values


@cache
def lower_order():
    d = parameters()
    actual = {name: d[name].subs(d["time_substitution"]) for name in ("f", "p")}
    f, p = actual["f"], actual["p"]
    fu, pu, px = sp.diff(f, U), sp.diff(p, U), sp.diff(p, X)
    q = sp.Function("q")(U, X)
    source = sp.Function("F_repository")(U, -X)
    cubic = -(q + fu) / (4 * p * X)
    affine_scalar = source + X * sp.diff(q, U) - 3 * X * (q + 2 * fu)**2 / (16 * p**2)
    numerator = pu - cubic * X
    q1 = -2 * fu + 4 * p * numerator
    corrected_q2 = 2 * fu / X - 4 * (p - 3 * X * px) * numerator / X
    printed_q2 = 2 * fu / X - 4 * (p - 3 * px) * numerator / X
    scalar = affine_scalar + 3 * X * numerator**2
    coefficient = 1 / (2 * X) - 3 * px / (2 * p)
    forcing = 3 * px * fu / p
    ode = sp.diff(q, X) + coefficient * q - forcing
    integrating_factor = sp.sqrt(-X) * p**sp.Rational(-3, 2)
    return {
        "u": U, "x": X, "f": f, "p": p, "f_u": fu, "q": q,
        "F3": cubic, "F2_parent": affine_scalar, "F_repository": source,
        "Q1": q1, "Q2_corrected": corrected_q2, "Q2_printed": printed_q2,
        "P": scalar, "ode_coefficient": coefficient, "ode_forcing": forcing,
        "ode_residual": ode, "integrating_factor": integrating_factor,
        "printed_error": _zero(printed_q2 - corrected_q2),
        "lower_order_claim_conditional_on_literal_Q2": True,
    }


def lower_order_checks():
    d = lower_order()
    q = d["q"]
    qx, qu = sp.diff(q, X), sp.diff(q, U)
    # The generic divergence is also checked directly in Fraction fixtures.
    box, vHv = sp.symbols("Box vHv", real=True)
    divergence = q * box + X * qu + 2 * qx * vHv
    combined = d["P"] + q * box + 2 * qx * vHv
    return {
        "Q1_is_q": _zero(d["Q1"] - q),
        "corrected_Q2_minus_2qx": _zero(d["Q2_corrected"] - 2 * qx + 2 * d["ode_residual"]),
        "P_scalar_counterterm": _zero(d["P"] - d["F_repository"] - X * qu),
        "full_IBP_scalar_match": _zero(combined - divergence - d["F_repository"]),
        "integrating_factor_log_derivative": _zero(
            sp.diff(d["integrating_factor"], X) / d["integrating_factor"] - d["ode_coefficient"]
        ),
        "initial_forcing_zero": _zero(d["ode_forcing"].subs(X, -1)),
        "initial_f_u_zero": _zero(d["f_u"].subs(X, -1)),
        "p_f_clock_derivative": _zero(d["f_u"] - 4 * d["p"] * sp.diff(d["p"], U)),
        "printed_error_formula": _zero(
            d["printed_error"] - 3 * (1 - X) * sp.diff(d["p"], X)
            * (q + 2 * d["f_u"]) / (d["p"] * X)
        ),
    }


def divergence_fixture(gradient, hessian, q, q_u, q_x, *, paper=True):
    """Product-rule divergence from individual lower components."""
    values = contractions(gradient, hessian, paper=paper)
    q, q_u, q_x = (_exact(item) for item in (q, q_u, q_x))
    v, tensor = _vector(gradient), _matrix(hessian)
    metric = (-1, 1, 1, 1) if paper else (1, -1, -1, -1)
    raised = tuple(metric[i] * v[i] for i in range(4))
    x_gradient = tuple(2 * sum(tensor[i][j] * raised[j] for j in range(4)) for i in range(4))
    direct = sum((q_u * v[i] + q_x * x_gradient[i]) * raised[i]
                 + q * metric[i] * tensor[i][i] for i in range(4))
    expanded = q_u * values["X"] + 2 * q_x * values["vHv"] + q * values["Box"]
    return {"direct": direct, "expanded": expanded, "residual": direct - expanded}


@cache
def palatini_control():
    """Independent pure-pR conformal completion, not the full affine inverse.

    For gbar=p*g, p*g^{-1}Ric[gbar] differs from p*R[g] by
    3(grad p)^2/(2p)-3 Box p.  Expanding grad p keeps the clock/X
    cross term and isolates Q2 without the general Eq.4.8 shortcut.
    """
    p = sp.Symbol("positive_p", positive=True)
    pu, px, mixed, fourth = sp.symbols("p_u p_x vHv L4", real=True)
    norm = pu**2 * X + 4 * pu * px * mixed + 4 * px**2 * fourth
    completion = 3 * norm / (2 * p)
    corrected = 2 * pu / X - 4 * (p - 3 * X * px) * pu / (2 * p * X)
    printed = 2 * pu / X - 4 * (p - 3 * px) * pu / (2 * p * X)
    return {
        "p": p, "p_u": pu, "p_x": px, "vHv": mixed, "L4": fourth,
        "completion": completion, "cross_coefficient": sp.diff(completion, mixed),
        "Q2_corrected": corrected, "Q2_printed": printed,
        "boundary": "-3 Box(p)",
    }


def palatini_checks():
    d = palatini_control()
    return {
        "pure_pR_Q2_cross_coefficient": _zero(
            d["cross_coefficient"] - 6 * d["p_u"] * d["p_x"] / d["p"]
        ),
        "pure_pR_corrected_Q2": _zero(d["Q2_corrected"] - d["cross_coefficient"]),
        "pure_pR_L4_coefficient": _zero(
            sp.diff(d["completion"], d["L4"]) - 6 * d["p_x"]**2 / d["p"]
        ),
    }


def domain_at(u, x):
    """Exact closed CD tube data; no connection-health verdict is returned."""
    u, x = _exact(u), _exact(x)
    if not Fraction(-11, 10) <= x <= Fraction(-9, 10):
        raise ValueError("x must be in the closed source-signature CD tube")
    h = (1 + u**2)**3
    w = h - 1 - x
    p_squared = w / (4 * h)
    return {
        "u": u, "x": x, "h": h, "w": w, "f": w / (2 * h),
        "p_squared": p_squared, "quotient_factor": 8 * p_squared - 1,
        "denominator": Fraction(1), "F4": (1 + x) / (2 * h * x**2),
        "q_abs_bound": Fraction(9, 310) / h**2,
    }


def calibration():
    return {
        "closed_x_min": Fraction(-11, 10), "closed_x_max": Fraction(-9, 10),
        "open_x_min": Fraction(-6, 5), "open_x_max": Fraction(-4, 5),
        "ratio_min": Fraction(9, 10), "ratio_max": Fraction(11, 10),
        "f_min": Fraction(9, 20), "f_max": Fraction(11, 20),
        "p_squared_min": Fraction(9, 40), "p_squared_max": Fraction(11, 40),
        "p_lower": Fraction(9, 20), "p_upper": Fraction(11, 20),
        "quotient_factor_min": Fraction(4, 5),
        "open_quotient_factor_lower": Fraction(3, 5),
        "ode_abs_coefficient": Fraction(25, 18),
        "ode_forcing_bound_times_h_squared": Fraction(1, 4),
        "initial_distance_max": Fraction(1, 10),
        "exponent_max": Fraction(5, 36),
        "exp_upper": Fraction(36, 31),
        "q_bound_times_h_squared": Fraction(9, 310),
        "initial_x": Fraction(-1),
        "domain": "all real u, closed -11/10<=x<=-9/10 inside open -6/5<x<-4/5",
        "units": "frozen CD M=tau=1 coefficient units; u is its scalar clock label",
        "source_Q2_correction_required": True,
        "connection_rank_or_health_proved_here": False,
        "UV_or_cutoff_claim": False,
    }


def domain_checks():
    d = parameters()
    h_of_u = (1 + U**2)**3
    return {
        "denominator_one": d["denominator"] - 1,
        "p_squared_relation": _zero(d["p"]**2 - d["f"] / 2),
        "ratio_lower_cleared": _zero(
            (d["w"] / H - sp.Rational(9, 10)) * H
            - ((H - 1) / 10 + (-X - sp.Rational(9, 10)))
        ),
        "ratio_upper_cleared": _zero(
            (sp.Rational(11, 10) - d["w"] / H) * H
            - ((H - 1) / 10 + (sp.Rational(11, 10) + X))
        ),
        "quotient_factor_relation": _zero(8 * d["p"]**2 - 1 - (2 * d["w"] / H - 1)),
        "p_log_x_derivative": _zero(d["p_x"] / d["p"] + 1 / (2 * d["w"])),
        "f_clock_derivative": _zero(
            sp.diff(d["f"].subs(H, h_of_u), U)
            - (1 + X) * sp.diff(h_of_u, U) / (2 * h_of_u**2)
        ),
        "time_derivative_majorant": _zero(
            9 * h_of_u**2 - sp.diff(h_of_u, U)**2
            - 9 * (1 + U**2)**4 * (U**2 - 1)**2
        ),
    }


def identities():
    values = {}
    for block in (principal_checks(), signature_checks(), lower_order_checks(),
                  palatini_checks(), domain_checks()):
        values.update(block)
    return values


def checks():
    d = calibration()
    endpoint = domain_at(0, Fraction(-9, 10))
    return {
        "closed_tube_negative": d["closed_x_max"] < 0,
        "closed_inside_open_left": d["open_x_min"] < d["closed_x_min"],
        "closed_inside_open_right": d["closed_x_max"] < d["open_x_max"],
        "initial_inside_closed": d["closed_x_min"] < d["initial_x"] < d["closed_x_max"],
        "positive_p": d["p_lower"] > 0,
        "p_lower_squared_margin": d["p_lower"]**2 < d["p_squared_min"],
        "p_upper_squared_margin": d["p_upper"]**2 > d["p_squared_max"],
        "positive_closed_quotient": d["quotient_factor_min"] > 0,
        "closed_quotient_endpoint_equality": endpoint["quotient_factor"] == d["quotient_factor_min"],
        "positive_open_quotient": d["open_quotient_factor_lower"] > 0,
        "exponential_series_domain": 0 < d["exponent_max"] < 1,
        "ode_exponent_product": d["ode_abs_coefficient"] * d["initial_distance_max"] == d["exponent_max"],
        "exponential_geometric_majorant": 1 / (1 - d["exponent_max"]) == d["exp_upper"],
        "q_bound_product": (d["ode_forcing_bound_times_h_squared"]
                            * d["initial_distance_max"] * d["exp_upper"]
                            == d["q_bound_times_h_squared"]),
        "q_bound_below_one_over_32": d["q_bound_times_h_squared"] < Fraction(1, 32),
        "printed_Q2_is_not_corrected_generically": lower_order()["printed_error"] != 0,
    }
