"""General HR-potential identities in the unchanged physical-g frame.

MG2 and MF2 are squared Planck masses: the curvature-square coefficient is
MF2**2*r**3/(4P), not MF2*r**3/(4P).  The polynomial APIs admit P=0; only
the stationary inverse requires a positive simple proportional root.
No functional inverse, heavy-state prescription or EFT error is supplied.
"""

from fractions import Fraction
from functools import cache

import sympy as sp

ETA = sp.diag(1, -1, -1, -1)
MG2, MF2, RATIO = sp.symbols("M_G_squared M_F_squared r", positive=True)
BETA = sp.symbols("beta0:5", real=True)
P = sp.Symbol("P", real=True, nonzero=True)
Z = sp.Symbol("Z", real=True)
S, SPHI, SPHIPHI = sp.symbols("s s_phi s_phiphi", real=True)
KAPPA, KP, KPP, KPPP = sp.symbols("kappa kappa_phi kappa_phiphi kappa_phiphiphi", real=True)


def _symmetric(prefix):
    entries = sp.symbols(f"{prefix}0:10", real=True)
    result = sp.zeros(4)
    pairs = tuple((i, j) for i in range(4) for j in range(i, 4))
    for index, (i, j) in enumerate(pairs):
        result[i, j] = result[j, i] = entries[index]
    return result, entries


def _binomial(n, k):
    return sp.binomial(n, k) if 0 <= k <= n else sp.S.Zero


def _rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("use a finite exact rational, not bool, float, text or a symbolic value")
    return sp.Rational(value)


def polynomial_at(beta, ratio):
    """Evaluate the polynomials on rational positive-r data; no inverse used.

    Neither F_root=0 nor P!=0 is imposed here.  This is an exact rational
    evaluation interface, not a restriction of the symbolic theorem to Q.
    """
    if not isinstance(beta, (tuple, list)) or len(beta) != 5:
        raise ValueError("beta must contain exactly beta0,...,beta4")
    coefficients = tuple(_rational(value) for value in beta)
    ratio = _rational(ratio)
    if ratio <= 0:
        raise ValueError("the specified proportional chart has r>0")
    b0, b1, b2, b3, b4 = coefficients
    p = b1+2*ratio*b2+ratio**2*b3
    return {"beta": coefficients, "r": ratio,
            "F_root": b1+3*ratio*b2+3*ratio**2*b3+ratio**3*b4,
            "P": p, "Z": b1+3*ratio*b2+2*ratio**2*b3,
            "stationary_potential_polynomial": b0+3*ratio*b1+3*ratio**2*b2+ratio**3*b3,
            "hessian_determinant_on_root": 3*(ratio*p)**10/16}


def inverse_at(beta, ratio, mg2=1, mf2=1):
    """Rational simple-root inverse data, with no spectral-health assertion."""
    result = polynomial_at(beta, ratio)
    mg2, mf2 = _rational(mg2), _rational(mf2)
    if mg2 <= 0 or mf2 <= 0:
        raise ValueError("both constant Einstein coefficients must be positive")
    if result["F_root"] != 0:
        raise ValueError("the supplied ratio is not a stationary proportional root")
    if result["P"] == 0:
        raise ValueError("P=0 is a nonsimple root with zero quadratic potential Hessian")
    r, p = result["r"], result["P"]
    return {**result, "MG2": mg2, "MF2": mf2,
            "root_derivative": -3*p/r, "relative_inverse_coefficient": mf2*r/p,
            "additive_inverse_coefficient": mf2*r**3/p,
            "kappa": mf2**2*r**3/(4*p), "leading_Planck_squared": mg2+mf2*r*r,
            "scope": "algebraic inverse only; P sign is not a parent-health verdict"}


@cache
def potential():
    """Literal elementary-symmetric expansion, all ten metric directions.

    e_n(r sqrt(I+epsilon H)) is expanded through epsilon^3 before imposing
    the stationary root.  The full invariant Hessian is not a TT-only fit.
    """
    h, entries = _symmetric("h")
    mixed = ETA*h
    tr1, tr2, tr3 = (sp.trace(mixed**n) for n in (1, 2, 3))
    r = RATIO
    series = {}
    for n in range(5):
        c1, c2, c3 = (_binomial(3, n-1), _binomial(2, n-2), _binomial(1, n-3))
        series[n] = (r**n*_binomial(4, n), r**n*c1*tr1/2,
                     r**n*(c2*tr1**2-(c1+c2)*tr2)/8,
                     r**n*(c3*tr1**3-3*(c2+c3)*tr1*tr2
                           +(3*(c1+c2)+2*c3)*tr3)/48)
    literal = tuple(-2*sum(BETA[n]*series[n][degree] for n in range(5))
                    for degree in range(4))
    b0, b1, b2, b3, b4 = BETA
    f_root = b1+3*r*b2+3*r*r*b3+r**3*b4
    p, z = b1+2*r*b2+r*r*b3, b1+3*r*b2+2*r*r*b3
    root_rule = {b4: -(b1+3*r*b2+3*r*r*b3)/r**3}
    quadratic = -r*P*(tr2-tr1**2)/4
    cubic = r*(Z*tr1**3-3*(P+Z)*tr1*tr2+(3*P+2*Z)*tr3)/24
    hessian = sp.hessian(quadratic, entries)
    return {"h": h, "entries": entries, "H": mixed, "traces": (tr1, tr2, tr3),
            "elementary_series": series, "literal_potential_coefficients": literal,
            "F_root": f_root, "P_polynomial": p, "Z_polynomial": z,
            "root_rule": root_rule, "root_derivative": sp.diff(f_root, r),
            "stationary_potential_polynomial": b0+3*r*b1+3*r*r*b2+r**3*b3,
            "quadratic": quadratic, "cubic": cubic, "hessian": hessian,
            "hessian_determinant": sp.factor(hessian.det())}


@cache
def stationary():
    """Own-f Einstein source, inverse and stationary fourth-order value."""
    p = potential()
    b, entries = _symmetric("B")
    trb = sp.trace(ETA*b)
    linear = MF2*RATIO**2*sp.trace(ETA*b*ETA*p["h"])/2
    solution = MF2*RATIO/P*(b-ETA*trb/3)
    substitution = dict(zip(p["entries"],
                            (solution[i, j] for i in range(4) for j in range(i, 4)),
                            strict=True))
    action = linear+p["quadratic"]
    kappa = MF2**2*RATIO**3/(4*P)
    expected = kappa*(sp.trace(ETA*b*ETA*b)-trb**2/3)
    ric, _ = _symmetric("Ric")
    ric_trace = sp.trace(ETA*ric)
    einstein = ric-ETA*ric_trace/2
    b_rule = dict(zip(entries, (einstein[i, j] for i in range(4) for j in range(i, 4)),
                      strict=True))
    return {"B": b, "entries_B": entries, "trace_g_B": trb,
            "linear_Einstein": linear, "quadratic_action": action,
            "relative_h2": solution, "additive_f2": RATIO**2*solution,
            "substitution": substitution, "kappa": kappa,
            "stationary_value": sp.expand(action.subs(substitution, simultaneous=True)),
            "expected_value": expected, "Ricci_f0": ric,
            "relative_Schouten": solution.subs(b_rule, simultaneous=True),
            "expected_relative_Schouten": MF2*RATIO/P*(ric-ETA*ric_trace/6)}


@cache
def leading_action():
    """Complete scalar normal form through formal degree four.

    The independent kappa*(Ricci^2-R_B^2/3) is retained, not discarded or
    absorbed into the displayed quadratic-DHOST coefficient names.
    """
    x, y, ric_scalar, rphi = sp.symbols("X Y R_B r_phi", real=True)
    p = potential()
    q = MG2+MF2*RATIO**2
    f4 = 2*(KP*(2*S**3+2*S*SPHI+SPHIPHI)
            +KPP*(S*S+2*SPHI)+S*KPPP)
    leading = (-q*ric_scalar/2+(sp.Rational(1, 2)-3*MF2*rphi*rphi)*x+y/2
               -2*p["stationary_potential_polynomial"])
    return {"X": x, "Y": y, "R_B": ric_scalar, "r_phi": rphi,
            "leading_Planck_squared": q, "leading_density": leading,
            "F2": -q/2-2*S*KP*x, "F2_X": -2*S*KP,
            "A1": -4*S*KP, "A2": 4*S*KP,
            "A3": sp.S.Zero, "A4": sp.S.Zero, "A5": sp.S.Zero,
            "K": 6*((S*S+SPHI)*KP+S*KPP)*x,
            "F": -2*p["stationary_potential_polynomial"]
                  +(sp.Rational(1, 2)-3*MF2*rphi*rphi)*x+f4*x*x,
            "f4": f4, "retained_curvature_square": "kappa*(Ricci^2-R_B^2/3)",
            "matter": "physical g and canonical free chi unchanged"}


@cache
def cubic_freedom():
    """Change Z while fixing r,P and every degree<=4 stationary coefficient."""
    p, eta, r = potential(), sp.Symbol("eta_deformation", real=True), RATIO
    increments = (3*r*eta, -2*eta, eta/r, sp.S.Zero, -eta/r**3)
    rule = {old: old+delta for old, delta in zip(BETA, increments, strict=True)}
    tr1, tr2, tr3 = p["traces"]
    cubic_shift = r*eta*(tr1**3-3*tr1*tr2+2*tr3)/24
    return {"eta": eta, "increments": increments, "rule": rule,
            "cubic_shift": cubic_shift,
            "scope": "formal degree-six action difference after the same h2 substitution"}


@cache
def nonboundary_control():
    """A nonzero Euler trace of the Z-deformation on an off-shell Einstein g.

    r,P,eta and MF2 are constant in this fixture.  The metric need not solve
    the parent equations.  A compact conformal variation is used, so this
    test does not infer non-boundary character merely from a density value.
    """
    lam, sigma, epsilon = sp.symbols("Lambda_Einstein sigma epsilon", real=True)
    hessian_sigma, _ = _symmetric("H_sigma")
    box_sigma = sp.trace(ETA*hessian_sigma)
    alpha = MF2*RATIO/P
    eta = cubic_freedom()["eta"]
    p = lam/3
    factor = RATIO*eta*alpha**3/4
    # Schouten_mixed=p I and delta Schouten_mixed=-2 sigma p I-2 Hess(sigma).
    mixed = p*sp.eye(4)+epsilon*(-2*sigma*p*sp.eye(4)-2*ETA*hessian_sigma)
    tr1, tr2, tr3 = (sp.trace(mixed**n) for n in (1, 2, 3))
    e3 = (tr1**3-3*tr1*tr2+2*tr3)/6
    density = 4*factor*p**3
    variation = factor*sp.diff((1+4*epsilon*sigma)*e3, epsilon).subs(epsilon, 0)
    expected = -2*sigma*density-6*factor*p*p*box_sigma
    hubble = sp.Symbol("H_deSitter", positive=True)
    return {"Lambda": lam, "sigma": sigma, "Box_sigma": box_sigma,
            "eta": eta, "alpha": alpha, "delta_L6": density,
            "literal_conformal_variation": variation, "expected_variation": expected,
            "covariant_metric_Euler_trace": -density,
            "H_deSitter": hubble,
            "deSitter_delta_L6": sp.factor(density.subs(lam, -3*hubble**2)),
            "scope": "nonzero for eta*Lambda!=0; off-shell physical g, not an actual parent solution"}


@cache
def curvature_linear_cubic():
    """Exact potential-only six-derivative curvature insertion and complement.

    This is V3[h2], not the full S6 action.  The named complement is kept
    explicitly; a coefficient projection is not an invariant of its full
    higher-derivative action under changes of that complement by IBP.
    """
    ric, _ = _symmetric("C_Ric")
    hess, _ = _symmetric("C_H")
    v = sp.Matrix(sp.symbols("v0:4", real=True))
    scalar = sp.trace(ETA*ric)
    g_tensor = ric-ETA*scalar/2
    x, trace_h = (v.T*ETA*v)[0], sp.trace(ETA*hess)
    h2, vhv = sp.trace(ETA*hess*ETA*hess), (v.T*ETA*hess*ETA*v)[0]
    a = ric-ETA*scalar/6
    d_tensor = -2*S*hess+2*(S*S-SPHI)*v*v.T-ETA*S*S*x
    da, d = sp.trace(ETA*a), sp.trace(ETA*d_tensor)
    d2 = sp.trace(ETA*d_tensor*ETA*d_tensor)
    # The exact coefficient linear in A of the full Newton cubic.
    raw = (Z*da*d*d-(P+Z)*(da*d2+2*d*sp.trace(ETA*a*ETA*d_tensor))
           +(3*P+2*Z)*sp.trace(ETA*a*ETA*d_tensor*ETA*d_tensor))
    gh = sp.trace(ETA*g_tensor*ETA*hess)
    gvv = (v.T*ETA*g_tensor*ETA*v)[0]
    operators = {
        "G_H_squared": sp.trace(ETA*g_tensor*ETA*hess*ETA*hess),
        "T_G_H": trace_h*gh, "X_G_H": x*gh,
        "G_v_Hv": (v.T*ETA*g_tensor*ETA*hess*ETA*v)[0],
        "T_G_vv": trace_h*gvv, "X_G_vv": x*gvv,
        "R_H_squared_minus_T_squared": scalar*(h2-trace_h**2),
        "R_vHv": scalar*vhv, "R_X_T": scalar*x*trace_h,
        "R_X_squared": scalar*x*x,
    }
    coefficients = {
        "G_H_squared": 4*S*S*(3*P+2*Z),
        "T_G_H": -8*S*S*(P+Z),
        "X_G_H": 4*S*(P*S*S-2*(P+Z)*SPHI),
        "G_v_Hv": -8*S*(3*P+2*Z)*(S*S-SPHI),
        "T_G_vv": 8*S*(P+Z)*(S*S-SPHI),
        "X_G_vv": 4*(S*S-SPHI)*(2*(P+Z)*S*S-P*SPHI),
        "R_H_squared_minus_T_squared": sp.Rational(4, 3)*(2*P+Z)*S*S,
        "R_vHv": -sp.Rational(8, 3)*(2*P+Z)*S*(S*S-SPHI),
        "R_X_T": sp.Rational(4, 3)*S*((P+2*Z)*S*S-2*(2*P+Z)*SPHI),
        "R_X_squared": (P+2*Z)*S**4-4*P*S*S*SPHI,
    }
    n = RATIO*(MF2*RATIO/P)**3/8
    complement = tuple(key for key in operators if key not in ("X_G_vv", "R_X_squared"))
    return {"Ricci": ric, "H_phi": hess, "v": v, "G": g_tensor,
            "R_B": scalar, "X": x, "T": trace_h, "vHv": vhv,
            "A_curvature": a, "D_clock": d_tensor, "N": n,
            "raw_over_N": raw, "operators": operators, "coefficients": coefficients,
            "table_over_N": sum(coefficients[key]*value for key, value in operators.items()),
            "retained_complement": complement,
            "s_zero_coefficient_X_G_vv": 4*n*P*SPHI**2,
            "H_free_projection_Xi": 4*n*(P*(S**4+S*S*SPHI+SPHI**2)-2*Z*S*S*SPHI)*x,
            "scope": "potential-only identity; complement retained; not full S6 or a matching bound"}


@cache
def checks():
    p, h, leading = potential(), stationary(), leading_action()
    rule = {P: p["P_polynomial"], Z: p["Z_polynomial"]}
    literal, tr1 = p["literal_potential_coefficients"], p["traces"][0]
    out = {
        "core01_literal_linear": sp.factor(literal[1]+RATIO*p["F_root"]*tr1),
        "core02_literal_quadratic": sp.factor(literal[2].subs(p["root_rule"])
                                                -p["quadratic"].subs(rule)),
        "core03_simple_root_derivative": sp.factor((p["root_derivative"]
                                                     +3*p["P_polynomial"]/RATIO).subs(p["root_rule"])),
        "core04_stationary_constant": sp.factor(literal[0].subs(p["root_rule"])
                                                 +2*p["stationary_potential_polynomial"]),
    }
    for index, variable in enumerate(p["entries"], start=5):
        out[f"core{index:02}_full_stationary_equation"] = sp.factor(
            sp.diff(h["quadratic_action"], variable).subs(h["substitution"], simultaneous=True))
    out["core15_stationary_value"] = sp.factor(h["stationary_value"]-h["expected_value"])
    out["core16_hessian_determinant"] = sp.factor(p["hessian_determinant"]-3*(RATIO*P)**10/16)
    out["core17_literal_cubic"] = sp.factor(literal[3].subs(p["root_rule"])-p["cubic"].subs(rule))
    out["core18_Horndeski_locus"] = leading["A1"]-2*leading["F2_X"]
    phi = sp.Symbol("phi", real=True)
    kappa, q = sp.Function("kap")(phi), sp.Function("q")(phi)
    j = 4*sp.diff(kappa, phi)*q
    m = sp.diff(j, phi)+4*sp.diff(kappa, phi)*q*q
    replacements = {q: S, sp.diff(q, phi): SPHI, sp.diff(q, phi, 2): SPHIPHI,
                    sp.diff(kappa, phi): KP, sp.diff(kappa, phi, 2): KPP,
                    sp.diff(kappa, phi, 3): KPPP}
    ibp_f4 = (sp.diff(m, phi)/2+4*sp.diff(kappa, phi)*q**3).xreplace(replacements)
    out["core19_weighted_IBP_F4"] = sp.expand(ibp_f4-leading["f4"])
    for index, entry in enumerate(h["relative_Schouten"]-h["expected_relative_Schouten"]):
        out[f"Schouten_component_{index:02}"] = sp.expand(entry)
    freedom = cubic_freedom()
    for name in ("F_root", "P_polynomial", "stationary_potential_polynomial"):
        out[f"deformation_preserves_{name}"] = sp.expand(p[name].subs(freedom["rule"], simultaneous=True)-p[name])
    out["deformation_shifts_Z"] = sp.expand(p["Z_polynomial"].subs(freedom["rule"], simultaneous=True)
                                            -p["Z_polynomial"]-freedom["eta"])
    out["deformation_cubic"] = sp.expand(p["cubic"].subs(Z, Z+freedom["eta"])
                                         -p["cubic"]-freedom["cubic_shift"])
    table = curvature_linear_cubic()
    out["full_Lorentzian_curvature_table"] = sp.expand(table["raw_over_N"]-table["table_over_N"])
    out["raw_s_zero_Z_cancellation"] = sp.expand(table["raw_over_N"].subs(S, 0)
                                                  -4*P*SPHI**2*table["operators"]["X_G_vv"])
    out["declared_H_free_projection"] = sp.expand(
        table["N"]*(table["coefficients"]["X_G_vv"]-4*table["coefficients"]["R_X_squared"])*table["X"]
        -table["H_free_projection_Xi"])
    bulk = nonboundary_control()
    out["deformation_compact_conformal_variation"] = sp.expand(
        bulk["literal_conformal_variation"]-bulk["expected_variation"])
    out["deformation_deSitter_value"] = sp.factor(
        bulk["deSitter_delta_L6"]+RATIO*bulk["eta"]*bulk["alpha"]**3*bulk["H_deSitter"]**6)
    return out


def calibration():
    positive = inverse_at((7, 3, Fraction(1, 2), Fraction(1, 4), Fraction(-9, 8)), 2, 2, 5)
    negative = inverse_at((7, -3, Fraction(-1, 2), Fraction(-1, 4), Fraction(9, 8)), 2, 2, 5)
    double = polynomial_at((0, 1, Fraction(-2, 3), Fraction(1, 3), 0), 1)
    return {"positive_P_fixture": positive, "negative_P_fixture": negative,
            "double_root_fixture": double, "core_identity_count": 19,
            "quartic_Xi": sp.S.Zero, "retained_action_order": 4,
            "first_omitted_action_order": 6,
            "root_scope": "positive simple stationary proportional root; not necessarily a vacuum or a rolling solution",
            "matching_scope": "retained-order no-CD only; quantitative omitted-symbol floor still required",
            "sixth_order_scope": "potential-only table and formal Z-deformation, not full S6 or a controlled EFT"}
