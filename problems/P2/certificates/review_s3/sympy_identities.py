"""S3 review: the reviewer's own sympy derivation of the identities behind zsol (no p2 imports).

S1  For a quartic Delta = -(Lam/3) prod (r - r_i) with distinct roots, r^2/Delta = sum A_i/(r - r_i)
    with A_i = r_i^2/Delta'(r_i)  (so r*' = r^2/Delta for r* = sum A_i log|r - r_i|).
S2  With R = exp(-i s w r*) Z (s = +-1), the radial operator factors:
    (Delta R')' + (w^2 r^4/Delta - L) R = exp(-i s w r*) [Delta Z'' + (Delta' - 2 i s w r^2) Z' - (2 i s w r + L) Z].
    Derived here by expanding R' and (Delta R')' by hand (chain rule with r*' = r^2/Delta), with
    Delta a generic function.
S3  Phase factor = 1.  With Z_+ = exp(i w (r* - c_+)) R_+, Z_c = exp(-i w (r* - c_c)) R_c:
    [Delta (Z_+ Z_c' - Z_+' Z_c) + 2 i w r^2 Z_+ Z_c] / (r_+ r_c)
      = Delta (R_+ R_c' - R_+' R_c) exp(-i w c_+ + i w c_c) / (r_+ r_c)   (the frozen W_l), exactly.
S4  Normalisation: exp(i w (r* - c_+)) (r - r_+)^(-i w A_+) -> 1 as r -> r_+, i.e. the factor
    exp(i w sum_{i != +} A_i log|(r - r_i)/(r_+ - r_i)|) is analytic at r_+ with value 1, and
    A_+ = 1/(2 kappa_+) (Delta'(r_+) = 2 kappa_+ r_+^2 > 0), A_c = -1/(2 kappa_c).  So Z_+ =
    g_+ * (analytic, = 1 at r_+): the zsol solution normalised by Z_+(r_+) = 1 is the frozen one.
S5  Indicial roots of the Z-equation at a simple root r_h: {0, 2 i s w r_h^2/Delta'(r_h)}
    = {0, i w/kappa_+} at r_+ (s = +1) and {0, i w/kappa_c} at r_c (s = -1): the exponent-0
    recursion breaks only on E_+ = {-i n kappa_+} resp. E_c.
S6  D4 from the Z-form: at l = 0, Z_h = 1 + w zeta_h + O(w^2) with Delta zeta_+' = i (r^2 - r_+^2),
    Delta zeta_c' = -i (r^2 - r_c^2); then dW_0/dw (0) = i (r_+^2 + r_c^2)/(r_+ r_c).
S7  Conjugation: if Z_h(r; -conj w) = conj Z_h(r; w) for real r, then W(-conj w) = conj W(w)
    for the factored formula (coefficient-level statement checked symbolically).
Writes sympy_identities.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent


def main():
    out = {}
    r = sp.symbols("r", real=True)
    w, L = sp.symbols("w L")
    # S1
    r1, r2, r3, r4, lam = sp.symbols("r1 r2 r3 r4 Lam", positive=True)
    rts = [r1, r2, r3, r4]
    Dq = -(lam / 3) * sp.Mul(*[(r - x) for x in rts])
    dDq = sp.diff(Dq, r)
    A = [x ** 2 / dDq.subs(r, x) for x in rts]
    out["S1_partial_fractions"] = sp.simplify(sum(A[i] / (r - rts[i]) for i in range(4)) - r ** 2 / Dq) == 0
    # S2 (hand expansion with a generic Delta)
    D = sp.Function("D")(r)
    Z = sp.Function("Z")(r)
    s = sp.symbols("s")
    E = sp.Function("E")(r)          # E = exp(-i s w r*), E' = -i s w (r^2/D) E
    dE = -sp.I * s * w * r ** 2 / D * E
    R = E * Z
    Rp = dE * Z + E * sp.diff(Z, r)
    DRp = D * Rp
    DRp_p = sp.diff(DRp, r).subs(sp.Derivative(E, r), dE)
    lhs = DRp_p + (w ** 2 * r ** 4 / D - L) * R
    rhs = E * (D * sp.diff(Z, r, 2) + (sp.diff(D, r) - 2 * sp.I * s * w * r ** 2) * sp.diff(Z, r)
               - (2 * sp.I * s * w * r + L) * Z)
    diff = sp.expand(lhs - rhs)
    out["S2_factored_equation_s=+1"] = sp.simplify(diff.subs(s, 1)) == 0
    out["S2_factored_equation_s=-1"] = sp.simplify(diff.subs(s, -1)) == 0
    # S3
    cp, cc = sp.symbols("c_p c_c", real=True)
    rp_, rc_ = sp.symbols("r_p r_c", positive=True)
    Ep = sp.Function("Ep")(r)        # exp(i w r*), Ep' = i w r^2/D Ep
    Rp_ = sp.Function("Rp")(r)
    Rc_ = sp.Function("Rc")(r)
    dEp = sp.I * w * r ** 2 / D * Ep
    Zp = Ep * sp.exp(-sp.I * w * cp) * Rp_
    Zc = sp.exp(sp.I * w * cc) * Rc_ / Ep
    dZp = sp.diff(Zp, r).subs(sp.Derivative(Ep, r), dEp)
    dZc = sp.diff(Zc, r).subs(sp.Derivative(Ep, r), dEp)
    Wz = (D * (Zp * dZc - dZp * Zc) + 2 * sp.I * w * r ** 2 * Zp * Zc) / (rp_ * rc_)
    Wf = D * (Rp_ * sp.diff(Rc_, r) - sp.diff(Rp_, r) * Rc_) * sp.exp(-sp.I * w * cp + sp.I * w * cc) / (rp_ * rc_)
    out["S3_factored_W_equals_frozen_W"] = sp.simplify(sp.expand(Wz - Wf)) == 0
    # S4: A_+ = 1/(2 kappa_+), A_c = -1/(2 kappa_c) with kappa = |Delta'|/(2 r^2) and the signs of
    # Delta' at r_+ (> 0) and r_c (< 0); the remaining factor -> 1 at r_+.
    kp, kc = sp.symbols("kappa_p kappa_c", positive=True)
    Ap = rp_ ** 2 / (2 * kp * rp_ ** 2)
    Ac = rc_ ** 2 / (-2 * kc * rc_ ** 2)
    out["S4_A_plus"] = sp.simplify(Ap - 1 / (2 * kp)) == 0
    out["S4_A_c"] = sp.simplify(Ac + 1 / (2 * kc)) == 0
    x = sp.symbols("x", positive=True)           # x = r - r_+
    ri = sp.symbols("r_i", real=True)
    Ai = sp.symbols("A_i", real=True)
    fac = sp.exp(sp.I * w * Ai * (sp.log(rp_ + x - ri) - sp.log(rp_ - ri)))
    out["S4_factor_to_1"] = sp.simplify(sp.limit(fac.subs({ri: 0, rp_: 1}), x, 0) - 1) == 0
    # S5 indicial roots
    m, d1, d2, rh = sp.symbols("m d1 d2 r_h")
    t = sp.symbols("t", positive=True)
    Dt = d1 * t + d2 * t ** 2
    Bt = sp.diff(Dt, t) - 2 * sp.I * s * w * (rh + t) ** 2
    lead = sp.expand(sp.powsimp(sp.expand((Dt * sp.diff(t ** m, t, 2) + Bt * sp.diff(t ** m, t)) / t ** (m - 1)), force=True))
    ind = sp.factor(sp.simplify(lead.subs(t, 0)))
    roots = sp.solve(sp.Eq(ind, 0), m)
    out["S5_indicial_roots"] = [str(x_) for x_ in roots]
    z_plus = [sp.simplify(x_.subs({s: 1, d1: 2 * kp * rh ** 2})) for x_ in roots]
    z_c = [sp.simplify(x_.subs({s: -1, d1: -2 * kc * rh ** 2})) for x_ in roots]
    out["S5_roots_at_r_plus"] = [str(x_) for x_ in z_plus]
    out["S5_roots_at_r_c"] = [str(x_) for x_ in z_c]
    out["S5_ok"] = (set(map(sp.simplify, z_plus)) == {0, sp.I * w / kp}
                    and set(map(sp.simplify, z_c)) == {0, sp.I * w / kc})
    # S6 D4 from the Z-form (O(w) terms)
    zp_d = sp.I * (r ** 2 - rp_ ** 2) / D       # zeta_+'
    zc_d = -sp.I * (r ** 2 - rc_ ** 2) / D      # zeta_c'
    # check the O(w) equation (D zeta')' = 2 i s r with these choices (and regularity: D zeta' -> 0 at r_h)
    out["S6_eq_plus"] = sp.simplify(sp.diff(D * zp_d, r) - 2 * sp.I * r) == 0
    out["S6_eq_c"] = sp.simplify(sp.diff(D * zc_d, r) + 2 * sp.I * r) == 0
    W1 = (D * (zc_d - zp_d) + 2 * sp.I * r ** 2) / (rp_ * rc_)
    out["S6_D4"] = sp.simplify(W1 - sp.I * (rp_ ** 2 + rc_ ** 2) / (rp_ * rc_)) == 0
    # S7 conjugation of the factored formula
    a_, b_ = sp.symbols("a b", real=True)
    zr = sp.symbols("zp0 zp1 zc0 zc1", complex=True)    # Z_+, Z_+', Z_c, Z_c' at w
    Dr, rr, rpr, rcr = sp.symbols("Dr rr rpr rcr", real=True)
    wv = a_ + sp.I * b_
    def Wform(wx, Zp0, Zp1, Zc0, Zc1):
        return (Dr * (Zp0 * Zc1 - Zp1 * Zc0) + 2 * sp.I * wx * rr ** 2 * Zp0 * Zc0) / (rpr * rcr)
    lhs7 = Wform(-sp.conjugate(wv), *[sp.conjugate(v) for v in zr])
    rhs7 = sp.conjugate(Wform(wv, *zr))
    out["S7_conjugation"] = sp.simplify(sp.expand(lhs7 - rhs7)) == 0
    out["pass"] = all(v for k, v in out.items() if isinstance(v, bool))
    (HERE / "sympy_identities.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
