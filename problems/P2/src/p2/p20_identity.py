"""P2-0 (mode stability on the closed upper half-plane) and D4 (simple zero mode):
the exact sympy identity chain.  FORMULATION.md section 2 (P2-0), D4; s0-gate.md section 2.

Every step below is an *exact* symbolic identity (sympy ``expand``/``cancel`` to 0) or an exact
substitution; nothing is quoted.  The written steps (integration of a total derivative over
[r_+, r_c]; sign bookkeeping) are recorded verbatim in the certificate and in
``notes/s1-p20-d1prime.md``.

Chain (u = e^{-i w T} R, T' = eta/F, w = a + i b, L = l(l+1)):
  C1  (Delta R')' + (w^2 r^4/Delta - L) R = e^{i w T} E[u],
      E[u] = (r^2 F u')' + i w (2 r^2 eta u' + (r^2 eta)' u) + w^2 r^2 (1 - eta^2) F^{-1} u - L u
      (generic F, eta, u; Delta = r^2 F).
  C2  conj(u) E[u] = J' - r^2 F |u'|^2 - 2 w r^2 eta Im(u' conj u) + w^2 r^2(1-eta^2)/F |u|^2 - L |u|^2,
      J = r^2 F u' conj(u) + i w r^2 eta |u|^2     (pointwise, u = p + i q, p, q real).
  C3  J(r_c) - J(r_+) = i w B,  B = r_c^2 |u(r_c)|^2 + r_+^2 |u(r_+)|^2   (F = 0, eta = -+1 there).
  C4  0 = i w B - D - 2 w X + w^2 W - L N  splits into
        Re: -D - b B - 2 a X + (a^2 - b^2) W - L N = 0,   Im: a B - 2 b X + 2 a b W = 0,
      D = int r^2 F |u'|^2, N = int |u|^2, X = int r^2 eta Im(u' conj u), W = int r^2(1-eta^2)/F |u|^2.
  C5  b Re - a Im = -[ b (D + L N + |w|^2 W) + |w|^2 B ]       ("(E)" without dividing by b).
  C6  for eta = (2r - r_+ - r_c)/(r_c - r_+):  r^2 (1 - eta^2)/F = 4 r^4 / ((Lambda/3)(r_c - r_+)^2 (r - r_n)(r - r_-))
      (> 0 and bounded on [r_+, r_c]).
  C7  T' = eta/F has residue -1/(2 kappa_+) at r_+ and -1/(2 kappa_c) at r_c; the E[u] indicial
      exponents are {0, i w/kappa_+} at r_+ and {0, i w/kappa_c} at r_c.  Hence
      u_+ := e^{-i w T} R_+ = e^{-i w T_reg} g_+ is analytic at r_+ with u_+(r_+) != 0.
  C8  (D4)  d/dw [Delta (R_+ R_c' - R_+' R_c)] at (l, w) = (0, 0) equals i (r_+^2 + r_c^2);
      so W_0'(0) = i (r_+^2 + r_c^2)/(r_+ r_c) != 0 and the zero of W_0 at 0 is simple.

Replay:  .venv/bin/python -m p2.p20_identity [write|verify]   (cwd problems/P2/src; ~10 s)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp

from .arbutil import Timer, git_commit

CERT = Path(__file__).resolve().parents[2] / "certificates" / "p2_0_identity.json"

r = sp.Symbol("r", positive=True)
a, b, Lam = sp.symbols("a b Lambda", real=True)
L = sp.Symbol("L", nonnegative=True)
w = sp.Symbol("omega")                      # complex frequency (C1)
Ff, eta, p, qf, uf, h = (sp.Function(n) for n in ("F", "eta", "p", "q", "u", "h"))


def _zero(expr) -> bool:
    return sp.simplify(sp.expand(expr)) == 0


def c1_transformed_ode(eps=0):
    """C1: generic F, eta, u (complex), h with h' = eta/F.  eps != 0: negative control."""
    F_, et, u = Ff(r), eta(r), uf(r)
    hh = h(r)
    R = sp.exp(sp.I * w * hh) * u
    lhs = sp.diff(r ** 2 * F_ * sp.diff(R, r), r) + (w ** 2 * r ** 2 / F_ - L) * R
    lhs = lhs.subs(sp.Derivative(hh, (r, 2)), sp.diff(et / F_, r)).subs(sp.Derivative(hh, r), et / F_)
    E = (sp.diff(r ** 2 * F_ * sp.diff(u, r), r)
         + sp.I * w * (2 * r ** 2 * et * sp.diff(u, r) + sp.diff(r ** 2 * et, r) * u)
         + w ** 2 * r ** 2 * (1 - et ** 2) / F_ * u - L * u + eps * sp.I * w * u)
    diff = sp.expand(sp.simplify(lhs * sp.exp(-sp.I * w * hh) - E))
    return diff == 0, "e^{-i w T}[(r^2 F R')' + (w^2 r^2/F - L) R] - E[u] == 0 for R = e^{i w T} u, T' = eta/F"


def _E_of(u, F_, et, ww):
    return (sp.diff(r ** 2 * F_ * sp.diff(u, r), r)
            + sp.I * ww * (2 * r ** 2 * et * sp.diff(u, r) + sp.diff(r ** 2 * et, r) * u)
            + ww ** 2 * r ** 2 * (1 - et ** 2) / F_ * u - L * u)


def c2_divergence_identity(eps=0):
    """C2: u = p + i q with p, q real; w = a + i b; F, eta real."""
    F_, et = Ff(r), eta(r)
    P, Qv = p(r), qf(r)
    u, ub = P + sp.I * Qv, P - sp.I * Qv
    ww = a + sp.I * b
    E = _E_of(u, F_, et, ww)
    up, ubp = sp.diff(u, r), sp.diff(ub, r)
    mod2 = P ** 2 + Qv ** 2
    im_upub = P * sp.diff(Qv, r) - Qv * sp.diff(P, r)          # Im(u' conj u)
    J = r ** 2 * F_ * up * ub + sp.I * ww * r ** 2 * et * mod2
    rhs = (sp.diff(J, r) - r ** 2 * F_ * (up * ubp) - 2 * ww * r ** 2 * et * im_upub
           + ww ** 2 * r ** 2 * (1 - et ** 2) / F_ * mod2 - L * mod2 + eps * ww * mod2)
    ok = _zero(ub * E - rhs)
    ok2 = _zero(sp.expand(up * ub) - (sp.diff(P, r) * P + sp.diff(Qv, r) * Qv + sp.I * im_upub))
    return ok and ok2, "conj(u) E[u] = J' - r^2F|u'|^2 - 2w r^2 eta Im(u' conj u) + w^2 r^2(1-eta^2)/F |u|^2 - L|u|^2"


def c3_boundary_flux(eps=0):
    """C3: J at the horizons with F = 0, eta(r_+) = -1, eta(r_c) = +1 (u, u' finite there)."""
    ww = a + sp.I * b
    rp, rc = sp.symbols("r_p r_c", positive=True)
    up_, uc_, dup, duc = sp.symbols("u_p u_c du_p du_c")        # values of u, u' (complex)
    upb, ucb = sp.conjugate(up_), sp.conjugate(uc_)

    def J(rr, Fv, etv, uval, ubar, du):
        return rr ** 2 * Fv * du * ubar + sp.I * ww * rr ** 2 * etv * uval * ubar

    flux = J(rc, 0, 1, uc_, ucb, duc) - J(rp, 0, -1, up_, upb, dup)
    B = rc ** 2 * uc_ * ucb + (1 + eps) * rp ** 2 * up_ * upb
    return _zero(flux - sp.I * ww * B), "J(r_c) - J(r_+) = i w (r_c^2|u(r_c)|^2 + r_+^2|u(r_+)|^2) = i w B"


def c4_c5_real_imag_and_elimination(eps=0):
    D, B, X, Wt, N = sp.symbols("D B X W N", real=True)
    ww = a + sp.I * b
    total = sp.I * ww * B - D - 2 * ww * X + ww ** 2 * Wt - L * N
    re, im = (sp.expand(t) for t in total.as_real_imag())
    ok_re = _zero(re - (-D - b * B - 2 * a * X + (a ** 2 - b ** 2) * Wt - L * N))
    ok_im = _zero(im - (a * B - 2 * b * X + 2 * a * b * Wt))
    comb = sp.expand(b * re - a * im)
    ok_el = _zero(comb + (b * (D + L * N + (a ** 2 + b ** 2 + eps) * Wt) + (a ** 2 + b ** 2) * B))
    if eps:
        ok_re = ok_re and _zero(re - (-D - b * B - 2 * a * X + (a ** 2 - b ** 2 + eps) * Wt - L * N))
    return (ok_re and ok_im), ok_el


def c6_weight_positive(eps=0):
    rn, rm, rp, rc = sp.symbols("r_n r_m r_p r_c", real=True)
    Fv = -(Lam / 3) * (r - rn) * (r - rm) * (r - rp) * (r - rc) / r ** 2
    et = (2 * r - rp - rc) / (rc - rp)
    lhs = r ** 2 * (1 - et ** 2) / Fv
    rhs = (4 + eps) * r ** 4 / ((Lam / 3) * (rc - rp) ** 2 * (r - rn) * (r - rm))
    ok = sp.cancel(sp.together(lhs - rhs)) == 0
    ok_ends = sp.simplify(et.subs(r, rp) + 1) == 0 and sp.simplify(et.subs(r, rc) - 1) == 0
    return ok and ok_ends


def c7_exponents(eps=0):
    """Residues of T' = eta r^2/Delta and the indicial exponents of E[u] at both horizons."""
    rp, rc, kp, kc = sp.symbols("r_p r_c kappa_p kappa_c", positive=True)
    s, ww = sp.symbols("s omega")
    D1 = sp.Function("Delta1")          # Delta = (r - r_p) Delta1(r), Delta1(r_p) = Delta'(r_p) = 2 kappa_p r_p^2
    D2 = sp.Function("Delta2")          # Delta = (r - r_c) Delta2(r), Delta2(r_c) = Delta'(r_c) = -2 kappa_c r_c^2
    et = (2 * r - rp - rc) / (rc - rp)
    res_p = sp.simplify((et * r ** 2 / D1(r)).subs(r, rp).subs(D1(rp), 2 * kp * rp ** 2))
    res_c = sp.simplify((et * r ** 2 / D2(r)).subs(r, rc).subs(D2(rc), -2 * kc * rc ** 2))
    ok_res = _zero(res_p + (1 + eps) / (2 * kp)) and _zero(res_c + 1 / (2 * kc))
    # indicial polynomial of p2 u'' + p1 u' + ... at a simple zero of p2 = Delta:
    # p2 ~ Delta'(r_h) (r - r_h), p1 -> Delta'(r_h) + 2 i w r_h^2 eta(r_h)
    def indicial(dprime, rh, eta_h):
        p2c, p1c = dprime, dprime + 2 * sp.I * ww * rh ** 2 * eta_h
        return sp.expand(s * (s - 1) * p2c + s * p1c)
    ind_p = indicial(2 * kp * rp ** 2, rp, -1)
    ind_c = indicial(-2 * kc * rc ** 2, rc, 1)
    ok_ind = (_zero(ind_p - 2 * kp * rp ** 2 * s * (s - (1 + eps) * sp.I * ww / kp))
              and _zero(ind_c + 2 * kc * rc ** 2 * s * (s - sp.I * ww / kc)))
    return ok_res, ok_ind


def c8_zero_mode_derivative(eps=0):
    """D4: with Delta = (r - r_p)(r - r_c) G(r), G(r_p), G(r_c) from kappa's.

    y = R_+(.; w), z = R_c(.; w); at w = 0 (l = 0) y = z = 1; Delta y_w' and Delta z_w' are constant
    (differentiate (Delta R')' + (w^2 r^4/Delta) R = 0 in w at w = 0) and equal to their horizon
    limits; y_w = -(i/(2 kappa_p)) ln(r - r_p) + (analytic), z_w = -(i/(2 kappa_c)) ln(r_c - r) + (analytic)."""
    rp, rc, kp, kc = sp.symbols("r_p r_c kappa_p kappa_c", positive=True)
    G = sp.Function("G")
    ww = sp.Symbol("omega")
    Delta = (r - rp) * (r - rc) * G(r)
    # Delta'(r_p) = (r_p - r_c) G(r_p) = 2 kappa_p r_p^2 ;  Delta'(r_c) = (r_c - r_p) G(r_c) = -2 kappa_c r_c^2
    Gp = 2 * kp * rp ** 2 / (rp - rc)
    Gc = -2 * kc * rc ** 2 / (rc - rp)
    ok_dp = _zero(sp.diff(Delta, r).subs(r, rp).subs(G(rp), Gp) - 2 * kp * rp ** 2)
    ok_dc = _zero(sp.diff(Delta, r).subs(r, rc).subs(G(rc), Gc) + 2 * kc * rc ** 2)
    # singular parts of the w-derivatives
    yw_sing = -(sp.I / (2 * kp)) * sp.log(r - rp)
    zw_sing = -(sp.I / (2 * kc)) * sp.log(rc - r)
    cp = sp.limit(sp.simplify(Delta * sp.diff(yw_sing, r)).subs(G(r), Gp), r, rp)
    cc = sp.limit(sp.simplify(Delta * sp.diff(zw_sing, r)).subs(G(r), Gc), r, rc)
    ok_cp = _zero(cp + sp.I * rp ** 2)
    ok_cc = _zero(cc - sp.I * rc ** 2)
    # W'(0) = Delta (z_w' - y_w') at y = z = 1: = cc - cp
    ok_w = _zero((cc - cp) - sp.I * (rp ** 2 + (1 + eps) * rc ** 2))
    # differentiated ODE at w = 0, l = 0: (Delta y_w')' = -(2 w r^4/Delta) y - (w^2 r^4/Delta) y_w -> 0
    yw = sp.Function("yw")(r)
    dode = (sp.diff(Delta * sp.diff(yw, r), r) + 2 * ww * r ** 4 / Delta * 1 + ww ** 2 * r ** 4 / Delta * yw)
    ok_ode = _zero(dode.subs(ww, 0) - sp.diff(Delta * sp.diff(yw, r), r))
    # the constant function 1 solves the l = 0, w = 0 radial equation
    ok_one = _zero(sp.diff(Delta * sp.diff(sp.Integer(1), r), r))
    return all((ok_dp, ok_dc, ok_cp, ok_cc, ok_w, ok_ode, ok_one))


WRITTEN = {
    "integration": "u = e^{-i w T} R_+ is analytic on [r_+, r_c] when W_l(w) = 0 (C7 at r_+; R_+ = const*R_c "
                   "and C7 at r_c), so D, N, X, W, B are finite (C6 bounds the W-weight) and integrating C2 "
                   "over [r_+, r_c] with C3 gives 0 = i w B - D - 2 w X + w^2 W - L N (C4).",
    "b_pos": "b > 0: by C5, b (D + L N + |w|^2 W) + |w|^2 B = 0 with every term >= 0, so W = 0 (w != 0); the "
             "weight of W is > 0 on (r_+, r_c) (C6), hence u = 0: no QNM with Im w > 0.",
    "b_zero": "b = 0, a != 0: Im part of C4 gives a B = 0, so u(r_+) = 0, contradicting u_+(r_+) != 0 (C7).  "
              "w = 0: C4 gives D + L N = 0, so u' = 0 and (l >= 1) u = 0; for l = 0, u = const: the zero mode.",
    "D4": "C8: W_0(0) = 0 (R = 1) and W_0'(0) = i (r_+^2 + r_c^2)/(r_+ r_c) != 0 (the exponential Jost factor is "
          "entire and does not change the order of the zero), so m_0(0) = 1.",
}


def run_all() -> dict:
    steps = []
    with Timer() as tm:
        ok, txt = c1_transformed_ode()
        steps.append({"id": "C1", "statement": txt, "method": "sympy simplify/expand == 0, generic F, eta, u, T", "ok": ok})
        ok, txt = c2_divergence_identity()
        steps.append({"id": "C2", "statement": txt, "method": "sympy expand == 0, u = p + i q, w = a + i b", "ok": ok})
        ok, txt = c3_boundary_flux()
        steps.append({"id": "C3", "statement": txt, "method": "exact substitution", "ok": ok})
        ok45, ok_el = c4_c5_real_imag_and_elimination()
        steps.append({"id": "C4", "statement": "Re/Im of 0 = i w B - D - 2wX + w^2 W - L N", "method": "sympy as_real_imag", "ok": ok45})
        steps.append({"id": "C5", "statement": "b Re - a Im = -[b(D + L N + |w|^2 W) + |w|^2 B]", "method": "sympy expand", "ok": ok_el})
        steps.append({"id": "C6", "statement": "r^2(1-eta^2)/F = 4 r^4/((Lambda/3)(r_c-r_+)^2 (r-r_n)(r-r_-)), eta(r_+)=-1, eta(r_c)=1",
                      "method": "sympy cancel", "ok": c6_weight_positive()})
        ok_res, ok_ind = c7_exponents()
        steps.append({"id": "C7a", "statement": "Res_{r_+} T' = -1/(2 kappa_+), Res_{r_c} T' = -1/(2 kappa_c)", "method": "exact substitution", "ok": ok_res})
        steps.append({"id": "C7b", "statement": "indicial exponents of E[u]: {0, i w/kappa_+} at r_+, {0, i w/kappa_c} at r_c", "method": "sympy expand", "ok": ok_ind})
        steps.append({"id": "C8", "statement": "D4: d/dw Delta(R_+R_c' - R_+'R_c) at (l,w)=(0,0) is i(r_+^2 + r_c^2)", "method": "sympy limit/expand", "ok": c8_zero_mode_derivative()})
    with Timer() as tc:
        controls = {
            "C1": not c1_transformed_ode(1)[0], "C2": not c2_divergence_identity(1)[0],
            "C3": not c3_boundary_flux(1)[0], "C4": not c4_c5_real_imag_and_elimination(1)[0],
            "C5": not c4_c5_real_imag_and_elimination(1)[1], "C6": not c6_weight_positive(1),
            "C7a": not c7_exponents(1)[0], "C7b": not c7_exponents(1)[1], "C8": not c8_zero_mode_derivative(1)}
    return {"claim": "P2-0 (closed upper half-plane: W_l(w) = 0 with Im w >= 0 iff (l, w) = (0, 0)) and D4 (m_0(0) = 1)",
            "level": "CERTIFIED (exact identities) + written steps",
            "steps": steps, "written_steps": WRITTEN,
            "negative_controls": {"meaning": "each identity with a deliberately wrong term must FAIL (True = it failed)",
                                  "results": controls, "all_discriminating": all(controls.values())},
            "all_ok": all(s["ok"] for s in steps) and all(controls.values()),
            "sympy": sp.__version__, "commit": git_commit(),
            "timing_s": {"wall": round(tm.wall + tc.wall, 2), "cpu": round(tm.cpu + tc.cpu, 2)}}


def write():
    cert = run_all()
    CERT.write_text(json.dumps(cert, indent=1) + "\n")
    return cert


def verify() -> list[str]:
    stored = json.loads(CERT.read_text())
    fresh = run_all()
    probs = []
    if not fresh["all_ok"]:
        probs.append("replay: some identity failed")
    for s_old, s_new in zip(stored["steps"], fresh["steps"], strict=True):
        if s_old["id"] != s_new["id"] or s_old["ok"] != s_new["ok"] or s_old["statement"] != s_new["statement"]:
            probs.append(f"step {s_old['id']}: stored {s_old['ok']} vs replay {s_new['ok']}")
    return probs


if __name__ == "__main__":
    if (sys.argv[1] if len(sys.argv) > 1 else "verify") == "write":
        c = write()
        for s in c["steps"]:
            print(s["id"], s["ok"])
        print("all_ok", c["all_ok"], c["timing_s"])
    else:
        pr = verify()
        print("p2_0_identity.json replay:", "OK" if not pr else pr)
        sys.exit(1 if pr else 0)
