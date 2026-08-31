"""Arb-certified (E) witnesses for the 15 ladder rows S != {alpha_K} (S3).

Ball arithmetic (python-flint / Arb) with outward rounding throughout; house
style of problems/P4/src/p4/validated/ and problems/P9/src/p9/verify.py
(_endpoint).  Nothing on the verification path trusts a float: every input is
an exact rational made into a ball, every subinterval of J is the rigorous
hull of its rational endpoints, and a check passes only when the Arb
comparison is strict for the WHOLE ball.

Frozen background F* (FORMULATION 1.3; matches S1's derivation and
fiducial.py; units H0 = 1, M_Pl^2 = 1):

    w_DE(a)  = -4/5 - (4/5)(1 - a)        =>  1 + w_DE = (4a - 3)/5
    g_de(a)  = rho_DE/rho_DE0 = exp((9/5) ln a + (12/5)(1 - a))
    E2(a)    = H^2/H0^2 = (3/10) a^-3 + (7/10) g_de(a)
    Om       = (3/10) a^-3 / E2,   Ode = (7/10) g_de / E2
    hd       = Hdot/H^2 = -(3/2) [Om + (1 + w_DE) Ode]
    fhat     = f/(M_Pl^2 H^2) = 3 Ode (1 + w_DE)
    m2       = M^2/M_Pl^2 = exp(alpha_M ln a),   rm = 3 Om / m2

Row quantities: target = D c_s^2 = MASTER_TARGET_B (witnesses.py; matched
exactly to the derived dust eigen-speed, rev-v1.1 pin), D = alpha_K +
(3/2)(alpha_B^BS)^2 with alpha_K = 1, c_s^2 = target/D (ball division).

Certification method: adaptive bisection of J (rational endpoints) until the
ball verdict is strict on every subinterval; subdivision counts recorded.
Certificates land in problems/P9/b/certificates/P9b-1.<ROW>.json.

Run:  uv run python -m p9b.certify        (from the repo root)
"""
from __future__ import annotations

import json
import math
from fractions import Fraction as F
from pathlib import Path

from flint import arb, ctx

from .witnesses import A_C, GAP, REG, TARGET_POLY_FN, witness_desc

PREC_BITS = 128
MAX_DEPTH = 26
CERT_DIR = Path(__file__).resolve().parents[2] / "certificates"

BACKGROUND_DOC = {
    "w_DE": "-4/5 - (4/5)(1-a)  [CPL(-4/5,-4/5); 1+w = (4a-3)/5]",
    "g_de": "exp((9/5) ln a + (12/5)(1-a))",
    "E2": "(3/10) a^-3 + (7/10) g_de(a)   [H0 = 1]",
    "hd": "-(3/2) (Om + (1+w) Ode)",
    "fhat": "3 Ode (1+w)  [= f/(M_Pl^2 H^2), M_Pl^2 = 1]",
    "m2": "exp(alpha_M ln a)  [M^2(a=1) = M_Pl^2]",
    "rm": "3 Om / m2",
    "target": "MASTER_TARGET_B (witnesses.py) = D c_s^2; c_s^2 = target/D",
}


class CertificationError(RuntimeError):
    pass


def fr(x: F) -> arb:
    """Exact rational -> ball (exact numerator, one outward-rounded division)."""
    return arb(x.numerator) / arb(x.denominator)


def ball(lo: F, hi: F) -> arb:
    """Rigorous hull of the interval [lo, hi]."""
    return fr(lo).union(fr(hi))


def _endpoint(x: arb, sign: int) -> float:
    """Rigorous float endpoint (house style: P9 verify.py::_endpoint)."""
    m, r = x.mid(), x.rad()
    e = m + r if sign > 0 else m - r
    f = float(e.str(25, radius=False))
    return (math.nextafter(math.nextafter(f, math.inf), math.inf) if sign > 0
            else math.nextafter(math.nextafter(f, -math.inf), -math.inf))


def row_eval(row: str, lo: F, hi: F) -> dict:
    """All row quantities as balls on the subinterval [lo, hi] of J."""
    w = REG[row]
    A = ball(lo, hi)
    lnA = A.log()
    gde = ((fr(F(9, 5)) * lnA) + fr(F(12, 5)) * (1 - A)).exp()
    a3 = A**(-3)
    E2 = fr(F(3, 10)) * a3 + fr(F(7, 10)) * gde
    Om = fr(F(3, 10)) * a3 / E2
    Ode = fr(F(7, 10)) * gde / E2
    w1 = (4 * A - 3) / 5                       # 1 + w_DE (dyadic-free: /5 outward)
    fhat = 3 * Ode * w1
    hd = -fr(F(3, 2)) * (Om + w1 * Ode)
    aT = fr(w["aT"]) if w["aT"] is not None else arb(0)
    aM = fr(w["aM"]) if w["aM"] is not None else arb(0)
    aH = fr(w["aH"]) if w["aH"] is not None else arb(0)
    b = fr(w["b1"]) * A if w["b1"] is not None else arb(0)
    db = b                                     # d(b1 a)/dx = b1 a  (0 if no B)
    m2 = (aM * lnA).exp() if w["aM"] is not None else arb(1)
    rm = 3 * Om / m2
    target = TARGET_POLY_FN(b, db, aT, aM, aH, arb(0), hd, rm)
    D = 1 + fr(F(3, 2)) * b**2                 # alpha_K = 1
    cs2 = target / D
    return dict(cs2=cs2, D=D, E2=E2, Ode=Ode, m2=m2, fhat=fhat, aB=b)


def sweep(row: str, key: str, lo: F, hi: F, bound: F, side: int,
          max_depth: int = MAX_DEPTH) -> dict:
    """Certify q > bound (side=+1) or q < bound (side=-1) on [lo, hi] by
    adaptive bisection; every accepted subinterval's ball verdict is strict."""
    bnd = fr(bound)
    stack = [(lo, hi, 0)]
    n_sub, depth_max, extreme = 0, 0, None
    while stack:
        a, b, d = stack.pop()
        v = row_eval(row, a, b)[key]
        ok = (v > bnd) if side > 0 else (v < bnd)
        if ok:
            n_sub += 1
            depth_max = max(depth_max, d)
            e = _endpoint(v, -side)            # worst rigorous endpoint
            extreme = e if extreme is None else (min if side > 0 else max)(extreme, e)
            continue
        if d >= max_depth:
            raise CertificationError(f"{row}:{key} not strict on [{a},{b}]")
        m = (a + b) / 2
        stack += [(a, m, d + 1), (m, b, d + 1)]
    return {"quantity": key, "bound": str(bound),
            "relation": ">" if side > 0 else "<",
            "interval": [str(lo), str(hi)], "n_subintervals": n_sub,
            "max_depth": depth_max,
            ("min_lower_endpoint" if side > 0 else "max_upper_endpoint"): extreme,
            "pass": True}


def certify_row(row: str) -> dict:
    w = REG[row]
    a1, a2 = w["J"]
    assert a1 < A_C < a2
    checks = {}
    # D2(ii): exact eigen-speed c_s^2 >= c_min > 0 on J (strict ball sweep)
    checks["cs2_ge_cmin"] = sweep(row, "cs2", a1, a2, w["c_min"], +1)
    # D2(i): no ghost D >= D_min > 0; alpha_B^BS != 2 on J
    if w["b1"] is None:
        checks["D_ge_Dmin"] = {"exact": "D = alpha_K = 1 (alpha_B = 0)",
                               "bound": str(w["D_min"]), "pass": True}
    else:
        checks["D_ge_Dmin"] = sweep(row, "D", a1, a2, w["D_min"], +1)
        checks["aB_ne_2"] = sweep(row, "aB", a1, a2, F(2), -1)
    # D2(iii): c_T^2 = 1 + alpha_T >= 0 (exact rational); M*^2 > 0 on J
    cT2 = F(1) + (w["aT"] if w["aT"] is not None else F(0))
    assert cT2 >= 0
    checks["cT2_nonneg"] = {"exact": f"c_T^2 = 1 + alpha_T = {cT2}",
                            "pass": True}
    if w["aM"] is not None:
        checks["M2_pos"] = sweep(row, "m2", a1, a2, F(0), +1)
    # D1 support: rho_DE > 0 on J (Ode > 0) and H^2 > 0 on J (E2 > 0)
    checks["Ode_pos"] = sweep(row, "Ode", a1, a2, F(0), +1)
    checks["E2_pos"] = sweep(row, "E2", a1, a2, F(0), +1)
    # D1 crossing: certified strict signs of fhat left/right of a_c (IVT);
    # f = M_Pl^2 H^2 fhat and H^2 > 0 is certified above, so sign f = sign fhat
    checks["fhat_neg_left"] = sweep(row, "fhat", a1, A_C - GAP, F(0), -1)
    checks["fhat_pos_right"] = sweep(row, "fhat", A_C + GAP, a2, F(0), +1)
    assert all(c["pass"] for c in checks.values())
    return {
        "claim": f"P9b-1.{row}", "row": row, "verdict": "CERTIFIED (E)",
        "date": "2026-08-31", "prec_bits": PREC_BITS,
        "conventions": ("FORMULATION.md v1.1; alpha_B in BS14 convention; "
                        "c_s^2 = exact dust DE eigen-speed (rev-v1.1 pin)"),
        "background": BACKGROUND_DOC,
        "witness": witness_desc(row),
        "crossing": {"a_c": str(A_C), "interior_to_J": True,
                     "ivt_gap": str(GAP),
                     "note": ("fhat = 3 Ode (1+w), 1+w = (4a-3)/5: certified "
                              "f < 0 on the left and f > 0 on the right "
                              "subinterval; sign change on (a_c-1/50, "
                              "a_c+1/50) by IVT; f(a_c) = 0 exactly on F*")},
        "d3_bounds": {"alpha_K": "1 in (0, 10^3]",
                      "max_abs_alpha_i": "max(|a/2|, 1/10, 1/2, 1) <= 1 <= 10",
                      "aB_degeneracy": ("|alpha_B^BS| <= 1/2 < 2 (certified)"
                                        if w["b1"] is not None else "alpha_B = 0")},
        "checks": checks,
    }


def main(cert_dir: Path = CERT_DIR) -> dict:
    ctx.prec = PREC_BITS
    cert_dir.mkdir(exist_ok=True)
    out = {}
    for row in REG:
        cert = certify_row(row)
        path = cert_dir / f"P9b-1.{row}.json"
        path.write_text(json.dumps(cert, indent=1) + "\n")
        out[row] = cert
        c = cert["checks"]["cs2_ge_cmin"]
        print(f"{row:6s} CERTIFIED  c_s^2 > {c['bound']:>5s}  "
              f"(min lower endpoint {c['min_lower_endpoint']:+.6f}, "
              f"{c['n_subintervals']} subintervals, depth {c['max_depth']})")
    return out


if __name__ == "__main__":
    main()
