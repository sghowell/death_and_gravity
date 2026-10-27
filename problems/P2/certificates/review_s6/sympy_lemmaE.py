"""S6 review, check B5/B6: the algebra of Lemma E re-derived in sympy (no p2 code).

E1  The phase-factored equations: R = e^{-i s w (x - c)} Z, x' = r^2/Delta, turns (Delta R')' + (w^2 r^4/Delta - L) R = 0
    into Delta Z'' + (Delta' - 2 i s w r^2) Z' - (2 i s w r + L) Z = 0 (s = +1 at r_+, -1 at r_c).
E2  Frobenius at r_+: the coefficient of t^(m-1) in L[t^m] is m(Delta_1 (m - 1) + B_0) = Delta_1 m (m - i w/kappa_+);
    the first step D_1 a_1 + C_0 = 0, D_1(omega0) = 0, C_0(omega0) = -(2 kappa_+ r_+ + L).
E3  The log term at omega0: Z = 1 + c Z2 log t + d_1 t + ...; the log-free t^0 coefficient is Delta_1 c + C_0, so
    c = (2 kappa_+ r_+ + L)/Delta_1 != 0.  Checked by expanding L[1 + c Z2 log t + d1 t] with a generic Z2.
E4  The residue: D_1 = -2 i r_+^2 (w - omega0), and (w - omega0) a_1 -> i (2 kappa_+ r_+ + L)/(2 r_+^2).
E5  The Wronskian form: Delta (R_+ R_c' - R_+' R_c) = e^{i w (c_+ - c_c)} [Delta (Z_+ Z_c' - Z_+' Z_c) + 2 i w r^2 Z_+ Z_c].
E6  The closed form (E-int), section 2.4(f): W[rR1, rR2] (d/dx) = Delta (R1 R2' - R1' R2); V_l Phi2 Psi dx =
    e^{kappa_+ c_+} (L + r F') Z2 dr/r; and the two relations of (e)-(f) give Wt = -(1/r_+) int (L + r F') Z2 dr/r.
E7  The limits at r_c: W[Phi~_c, Phi_c] -> 2 kappa_+ r_c^2 and W[Phi~_c, Psi] -> 2 kappa_+ r_c e^{kappa_+ c_c}
    (with H(0) = H~(0) = r_c), and W[Phi_c, Psi] = 2 kappa_c sigma H' e^{kappa_+ (2x - c_c)}.
E8  r F' = (2/r)(M - Q^2/r) - (2 Lambda/3) r^2 for F = 1 - 2M/r + Q^2/r^2 - Lambda r^2/3.
E9  W6' (regime L): in the tau-scaled equation D^2 R_tt + s D E R_t + (w^2 r^4 - (1 - h^2/4) D) R = 0 (r = r_+ + h^2 tau,
    Delta = h^2 D), the pole-free series R~ = (1 + 2 rho) R has, at omega0 (rho = -1/2),
    a~_1 = -(A_3 rho (rho - 1) + B_2 rho + C_1)/A_2 = (1 + 2 kappa_+ r_+ h^2 - h^2/4)/Delta'(r_+) = h^2 (2 kappa_+ r_+ + L)/Delta'(r_+),
    nonzero for every h (the scaled log coefficient); at h = 0 it is 1/Delta'(r_+).

Run:  python sympy_lemmaE.py   -> sympy_lemmaE_out.json
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent


def e1():
    r, w, L, c = sp.symbols("r w L c")
    out = {}
    for s in (1, -1):
        D = sp.Function("Delta")(r)
        Z = sp.Function("Z")(r)
        x = sp.Function("x")(r)
        R = sp.exp(-sp.I * s * w * (x - c)) * Z
        expr = sp.diff(D * sp.diff(R, r), r) + (w ** 2 * r ** 4 / D - L) * R
        expr = expr.subs(sp.Derivative(x, (r, 2)), sp.diff(r ** 2 / D, r)).subs(sp.Derivative(x, r), r ** 2 / D)
        expr = sp.simplify(expr * sp.exp(sp.I * s * w * (x - c)))
        target = D * sp.diff(Z, r, 2) + (sp.diff(D, r) - 2 * sp.I * s * w * r ** 2) * sp.diff(Z, r) \
            - (2 * sp.I * s * w * r + L) * Z
        out[f"s={s}"] = bool(sp.simplify(sp.expand(expr - target)) == 0)
    return out


def frob_setup():
    t, w, L, rp, kp = sp.symbols("t w L r_p kappa_p")
    d2, d3, d4 = sp.symbols("d2 d3 d4")
    d1 = 2 * kp * rp ** 2                                    # Delta'(r_+) = 2 kappa_+ r_+^2
    Delta = d1 * t + d2 * t ** 2 + d3 * t ** 3 + d4 * t ** 4
    r = rp + t
    B = sp.diff(Delta, t) - 2 * sp.I * w * r ** 2
    C = -(2 * sp.I * w * r + L)
    Lop = lambda Z: Delta * sp.diff(Z, t, 2) + B * sp.diff(Z, t) + C * Z
    return t, w, L, rp, kp, d1, d2, d3, d4, Delta, B, C, Lop


def e2():
    t, w, L, rp, kp, d1, d2, d3, d4, Delta, B, C, Lop = frob_setup()
    m = sp.symbols("m", positive=True, integer=True)
    out = {}
    B0 = B.subs(t, 0)
    # lowest coefficient of L[t^m] is t^(m-1) with coefficient m(d1 (m-1) + B0)
    ok = True
    for mm in range(1, 6):
        e = sp.expand(Lop(t ** mm))
        c = e.coeff(t, mm - 1)
        ok = ok and sp.simplify(c - mm * (d1 * (mm - 1) + B0)) == 0
        ok = ok and all(sp.simplify(e.coeff(t, k)) == 0 for k in range(0, mm - 1))
    out["lowest_coefficient_m(Delta1(m-1)+B0)_m=1..5"] = bool(ok)
    Dm = m * (d1 * (m - 1) + B0)
    out["D_m = Delta1 m (m - i w/kappa_+)"] = bool(sp.simplify(Dm - d1 * m * (m - sp.I * w / kp)) == 0)
    a1 = sp.symbols("a1")
    e = sp.expand(Lop(1 + a1 * t))
    c0 = e.coeff(t, 0)
    D1 = B0
    C0 = C.subs(t, 0)
    out["t0_coefficient_is_D1_a1_plus_C0"] = bool(sp.simplify(c0 - (D1 * a1 + C0)) == 0)
    w0 = -sp.I * kp
    out["D1(omega0) = 0"] = bool(sp.simplify(D1.subs(w, w0)) == 0)
    out["C0(omega0) = -(2 kappa_+ r_+ + L)"] = bool(sp.simplify(C0.subs(w, w0) + 2 * kp * rp + L) == 0)
    out["D_m(omega0) = Delta1 m (m-1) != 0 for m >= 2"] = bool(sp.simplify(Dm.subs(w, w0) - d1 * m * (m - 1)) == 0)
    out["exponents at omega0"] = str(sp.solve(sp.Symbol("rho") * (d1 * (sp.Symbol("rho") - 1) + B0.subs(w, w0)),
                                              sp.Symbol("rho")))
    return out


def e3():
    """log coefficient: L[1 + c Z2 log t + d1 t] at omega0, with Z2 = t + b2 t^2 + b3 t^3 solving L[Z2] = O(t^3)"""
    t, w, L, rp, kp, d1, d2, d3, d4, Delta, B, C, Lop = frob_setup()
    w0 = -sp.I * kp
    b2, b3, cc, dd1 = sp.symbols("b2 b3 c dd1")
    LopW = lambda Z: sp.expand((Delta * sp.diff(Z, t, 2) + B * sp.diff(Z, t) + C * Z).subs(w, w0))
    Z2 = t + b2 * t ** 2 + b3 * t ** 3
    e = LopW(Z2)
    sol = sp.solve([e.coeff(t, 1), e.coeff(t, 2)], [b2, b3], dict=True)[0]
    Z2s = Z2.subs(sol)
    Y = 1 + cc * Z2s * sp.log(t) + dd1 * t
    e = LopW(Y)
    e = sp.expand(sp.simplify(e))
    # split log and non-log parts; the t^0 non-log coefficient
    e0 = e.subs(sp.log(t), 0)
    e_log = sp.expand((e - e0) / sp.log(t))
    nonlog_t0 = sp.simplify(sp.series(e0, t, 0, 1).removeO().subs(t, 0))
    out = {"nonlog_t0_coefficient": str(sp.factor(nonlog_t0)),
           "equals_Delta1_c_plus_C0": bool(sp.simplify(nonlog_t0 - (d1 * cc - (2 * kp * rp + L))) == 0),
           "c = (2 kappa_+ r_+ + L)/Delta_1": str(sp.solve(nonlog_t0, cc)),
           "log_part_O(t^2)": bool(sp.simplify(sp.series(e_log, t, 0, 2).removeO()) == 0)}
    return out


def e4():
    w, L, rp, kp = sp.symbols("w L r_p kappa_p")
    d1 = 2 * kp * rp ** 2
    w0 = -sp.I * kp
    D1 = d1 - 2 * sp.I * w * rp ** 2
    C0 = -(2 * sp.I * w * rp + L)
    a1 = -C0 / D1
    lim = sp.limit((w - w0) * a1, w, w0)
    return {"D1 = -2 i r_+^2 (w - omega0)": bool(sp.simplify(D1 + 2 * sp.I * rp ** 2 * (w - w0)) == 0),
            "lim (w - omega0) a1": str(sp.simplify(lim)),
            "equals i(2 kappa_+ r_+ + L)/(2 r_+^2)": bool(sp.simplify(lim - sp.I * (2 * kp * rp + L) / (2 * rp ** 2)) == 0)}


def e5():
    r, w, cp, cc = sp.symbols("r w c_p c_c")
    D = sp.Function("Delta")(r)
    x = sp.Function("x")(r)
    Zp, Zc = sp.Function("Zp")(r), sp.Function("Zc")(r)
    Rp = sp.exp(-sp.I * w * (x - cp)) * Zp
    Rc = sp.exp(sp.I * w * (x - cc)) * Zc
    lhs = D * (Rp * sp.diff(Rc, r) - sp.diff(Rp, r) * Rc)
    lhs = lhs.subs(sp.Derivative(x, r), r ** 2 / D)
    rhs = sp.exp(sp.I * w * (cp - cc)) * (D * (Zp * sp.diff(Zc, r) - sp.diff(Zp, r) * Zc) + 2 * sp.I * w * r ** 2 * Zp * Zc)
    return {"Wronskian_form": bool(sp.simplify(sp.expand(lhs - rhs)) == 0)}


def e6():
    r, L, M, Q, lam, kp, cp, cc, rp, rc = sp.symbols("r L M Q Lambda kappa_p c_p c_c r_p r_c", positive=True)
    F = 1 - 2 * M / r + Q ** 2 / r ** 2 - lam * r ** 2 / 3
    D = r ** 2 * F
    R1, R2 = sp.Function("R1")(r), sp.Function("R2")(r)
    # W[rR1, rR2] with d/dx = F d/dr
    ddx = lambda f: F * sp.diff(f, r)
    W = (r * R1) * ddx(r * R2) - ddx(r * R1) * (r * R2)
    ok1 = sp.simplify(W - D * (R1 * sp.diff(R2, r) - sp.diff(R1, r) * R2)) == 0
    # V_l Phi2 Psi dx with Phi2 = r e^{-kp (x - c_p)} Z2, Psi = e^{kp x}, dx = dr/F
    x = sp.Symbol("x")
    Z2 = sp.Function("Z2")(r)
    V = F * (L / r ** 2 + sp.diff(F, r) / r)
    integrand = V * (r * sp.exp(-kp * (x - cp)) * Z2) * sp.exp(kp * x) / F
    ok2 = sp.simplify(integrand - sp.exp(kp * cp) * (L + r * sp.diff(F, r)) * Z2 / r) == 0
    # the final assembly: e^{kp (cp - cc)} r_p r_c Wt = 2 kp r_c^2 B and 2 kp r_c e^{kp cc} B = -e^{kp cp} I
    Wt, Bc, Iv = sp.symbols("Wt B I")
    sol = sp.solve([sp.exp(kp * (cp - cc)) * rp * rc * Wt - 2 * kp * rc ** 2 * Bc,
                    2 * kp * rc * sp.exp(kp * cc) * Bc + sp.exp(kp * cp) * Iv], [Wt, Bc], dict=True)[0]
    ok3 = sp.simplify(sol[Wt] + Iv / rp) == 0
    # E8: r F'
    ok4 = sp.simplify(r * sp.diff(F, r) - ((2 / r) * (M - Q ** 2 / r) - 2 * lam * r ** 2 / 3)) == 0
    return {"W[rR1,rR2] = Delta (R1 R2' - R1' R2)": bool(ok1),
            "V Phi2 Psi dx = e^{kp c_p} (L + r F') Z2 dr/r": bool(ok2),
            "Wt = -I/r_+": bool(ok3), "Wt": str(sol[Wt]),
            "E8 r F' identity": bool(ok4)}


def e7():
    x, kp, kc, cc, rc = sp.symbols("x kappa_p kappa_c c_c r_c", positive=True)
    s = sp.Symbol("sigma")
    h0, h1, g0, g1 = sp.symbols("h0 h1 g0 g1")
    H = h0 + h1 * s                     # H(sigma) to first order (H(0) = h0)
    Ht = g0 + g1 * s
    sig = sp.exp(-2 * kc * (x - cc))
    Phic = sp.exp(kp * (x - cc)) * H.subs(s, sig)
    Phit = sp.exp(-kp * (x - cc)) * Ht.subs(s, sig)
    Psi = sp.exp(kp * x)
    Wr = lambda f, g: sp.simplify(f * sp.diff(g, x) - sp.diff(f, x) * g)
    W1 = sp.simplify(Wr(Phit, Phic))
    W2 = sp.simplify(Wr(Phit, Psi))
    W3 = sp.simplify(Wr(Phic, Psi))
    lim1 = sp.limit(W1.subs({h0: rc, g0: rc}), x, sp.oo) if False else sp.simplify(W1.subs(sig, 0))
    # the sigma-independent parts (sigma -> 0 as x -> oo)
    W1_0 = sp.simplify(sp.expand(W1).subs({h1: 0, g1: 0}))
    W2_0 = sp.simplify(sp.expand(W2).subs({g1: 0}))
    W3_expected = 2 * kc * sig * h1 * sp.exp(kp * (2 * x - cc))
    return {"W[Phi~_c, Phi_c] (sigma^0 part)": str(W1_0),
            "equals 2 kp h0 g0 (-> 2 kappa_+ r_c^2)": bool(sp.simplify(W1_0 - 2 * kp * h0 * g0) == 0),
            "W[Phi~_c, Psi] (sigma^0 part)": str(sp.simplify(W2_0)),
            "equals 2 kp g0 e^{kp c_c} (-> 2 kappa_+ r_c e^{kappa_+ c_c})": bool(sp.simplify(W2_0 - 2 * kp * g0 * sp.exp(kp * cc)) == 0),
            "W[Phi_c, Psi] = 2 kc sigma H' e^{kp(2x - c_c)}": bool(sp.simplify(W3 - W3_expected) == 0)}


def e9():
    tau, h2, kp, rp = sp.symbols("tau h2 kappa_p r_p", positive=True)
    d2, d3, d4 = sp.symbols("d2 d3 d4")
    d1 = 2 * kp * rp ** 2
    s = 1
    w = -sp.I * kp
    D = d1 * s * tau + d2 * h2 * tau ** 2 + d3 * s * h2 ** 2 * tau ** 3 + d4 * h2 ** 3 * tau ** 4
    r = rp + s * h2 * tau
    E = d1 + 2 * d2 * (r - rp) + 3 * d3 * (r - rp) ** 2 + 4 * d4 * (r - rp) ** 3
    A = sp.Poly(sp.expand(D * D), tau)
    B = sp.Poly(sp.expand(s * D * E), tau)
    C = sp.Poly(sp.expand(w ** 2 * r ** 4 - (1 - h2 / 4) * D), tau)
    co = lambda P, k: P.coeff_monomial(tau ** k)
    rho = sp.Rational(-1, 2)
    indicial_ok = sp.simplify(co(A, 2) * rho * (rho - 1) + co(B, 1) * rho + co(C, 0)) == 0
    a1 = -(co(A, 3) * rho * (rho - 1) + co(B, 2) * rho + co(C, 1)) / co(A, 2)
    target = (1 + 2 * kp * rp * h2 - h2 / 4) / d1
    return {"rho = -1/2 is an indicial root at omega0": bool(indicial_ok),
            "a~_1 at omega0": str(sp.simplify(a1)),
            "a~_1 = (1 + 2 kappa_+ r_+ h^2 - h^2/4)/Delta'(r_+)": bool(sp.simplify(a1 - target) == 0),
            "a~_1 at h = 0": str(sp.simplify(a1.subs(h2, 0)))}


def main():
    out = {"E1": e1(), "E2": e2(), "E3": e3(), "E4": e4(), "E5": e5(), "E6": e6(), "E7": e7(), "E9": e9()}
    flat = []
    for k, v in out.items():
        for kk, vv in v.items():
            if isinstance(vv, bool):
                flat.append(vv)
    out["pass"] = all(flat)
    (HERE / "sympy_lemmaE_out.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
