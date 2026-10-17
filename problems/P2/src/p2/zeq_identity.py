"""Exact (sympy) checks of the phase-factored formulation used by zsol (S3).

Z1. With r*' = r^2/Delta and R = exp(-i s w r*) Z (s = +1 or -1):
      (Delta R')' + (w^2 r^4/Delta - L) R
        = exp(-i s w r*) [Delta Z'' + (Delta' - 2 i s w r^2) Z' - (2 i s w r + L) Z].
Z2. With R_+ = exp(-i w (r* - c_+)) Z_+ and R_c = exp(+i w (r* - c_c)) Z_c:
      Delta (R_+ R_c' - R_+' R_c) exp(-i w c_+ + i w c_c)
        = Delta (Z_+ Z_c' - Z_+' Z_c) + 2 i w r^2 Z_+ Z_c.
Z3. Indicial factor at a simple root r_h of Delta (Delta(r_h) = 0, Delta'(r_h) = 2 k_h r_h^2 sgn):
      the coefficient of a_m is m (Delta'(r_h)(m - 1) + Delta'(r_h) - 2 i s w r_h^2)
        = Delta'(r_h) m (m - 2 i s w r_h^2/Delta'(r_h)).
Z4. Conjugation: the coefficients of the s-equation are mapped to their complex conjugates by
    w -> -conj(w) for real r (so W(-conj w) = conj W(w), as Z_h(r_h) = 1 is real).
Run: python -m p2.zeq_identity  (writes certificates/part3/zeq_identity.json).
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]


def checks() -> dict:
    r, w, L = sp.symbols("r w L")
    s = sp.symbols("s")
    cp, cc = sp.symbols("c_p c_c", real=True)
    D = sp.Function("Delta")(r)
    rs = sp.Function("rstar")(r)
    Z = sp.Function("Z")(r)
    Zp = sp.Function("Zp")(r)
    Zc = sp.Function("Zc")(r)
    drs = r ** 2 / D

    def d(expr):
        return sp.diff(expr, r).subs(sp.Derivative(rs, r), drs)

    out = {}
    # Z1
    R = sp.exp(-sp.I * s * w * rs) * Z
    lhs = d(D * d(R)) + (w ** 2 * r ** 4 / D - L) * R
    rhs = sp.exp(-sp.I * s * w * rs) * (D * sp.diff(Z, r, 2) + (sp.diff(D, r) - 2 * sp.I * s * w * r ** 2) * sp.diff(Z, r)
                                       - (2 * sp.I * s * w * r + L) * Z)
    # s^2 = 1 enters through w^2 r^4/Delta: check both signs explicitly
    z1p = sp.simplify(sp.expand((lhs - rhs).subs(s, 1)))
    z1m = sp.simplify(sp.expand((lhs - rhs).subs(s, -1)))
    out["Z1_s=+1"] = z1p == 0
    out["Z1_s=-1"] = z1m == 0
    # Z2
    Rp = sp.exp(-sp.I * w * (rs - cp)) * Zp
    Rc = sp.exp(sp.I * w * (rs - cc)) * Zc
    lhs2 = D * (Rp * d(Rc) - d(Rp) * Rc) * sp.exp(-sp.I * w * cp + sp.I * w * cc)
    rhs2 = D * (Zp * sp.diff(Zc, r) - sp.diff(Zp, r) * Zc) + 2 * sp.I * w * r ** 2 * Zp * Zc
    out["Z2"] = sp.simplify(sp.expand(lhs2 - rhs2)) == 0
    # Z3: Delta = d1 t + d2 t^2 + ..., Z = t^m: lowest power t^(m-1)
    t, m, d1, d2, rh = sp.symbols("t m d1 d2 r_h")
    Dt = d1 * t + d2 * t ** 2
    Bt = sp.diff(Dt, t) - 2 * sp.I * s * w * (rh + t) ** 2
    expr = sp.expand(sp.powsimp(sp.expand((Dt * sp.diff(t ** m, t, 2) + Bt * sp.diff(t ** m, t)) * t ** (1 - m)),
                                force=True))
    low = sp.simplify(expr.subs(t, 0))
    out["Z3"] = sp.simplify(low - d1 * m * (m - 2 * sp.I * s * w * rh ** 2 / d1)) == 0
    # Z4: B(r) and C(r) with w -> -conj(w) equal conj of B, C (r, Delta real)
    a, b = sp.symbols("a b", real=True)
    rr = sp.symbols("rr", real=True)
    Dp = sp.symbols("Dp", real=True)
    Lr = sp.symbols("Lr", real=True)
    ww = a + sp.I * b
    ok = True
    for sg in (1, -1):
        B = Dp - 2 * sp.I * sg * ww * rr ** 2
        C = -(2 * sp.I * sg * ww * rr + Lr)
        Bm = B.subs(a, -a)                      # -conj(a + i b) = -a + i b
        Cm = C.subs(a, -a)
        ok &= sp.simplify(sp.expand(Bm - sp.conjugate(B))) == 0
        ok &= sp.simplify(sp.expand(Cm - sp.conjugate(C))) == 0
    out["Z4"] = bool(ok)
    out["pass"] = all(v for k, v in out.items())
    return out


if __name__ == "__main__":
    res = checks()
    print(res)
    p = ROOT / "certificates" / "part3"
    p.mkdir(parents=True, exist_ok=True)
    (p / "zeq_identity.json").write_text(json.dumps(res, indent=1) + "\n")
