"""S4: P9b-2 — Melville-Noller positivity (POS) evaluation of the S3 witnesses.

POS pin (frozen, FORMULATION 2(b2); re-verified against the ar5iv render of
arXiv:1904.05874 on 2026-08-31; FULL pin + applicability map in
notes/s4-positivity.md 1-2):
- MN (16) c_ss >= 0, c_sst >= 0 [(5) with Lambda_2 >> Lambda_3; (14)-(15)],
  general Horndeski; MN (17) shift-symmetric: 2 G2X G4X >= -G2XX G4 and
  2 G4XX + 2 G4X^2/G4 <= G3X^2 (overbars: flat-background values).
- MN (10), example model (their (3), shift-symmetric quartic Horndeski):
      pos. prior:  alpha_B <= 2 alpha_T / (1 + alpha_T),
  from the c_sst part of (16) via their (8)-(9), using X > 0, M^2 > 0,
  G4 > 0 <=> 1 + alpha_T > 0 (side condition checked per row).  MN's
  alpha_B IS alpha_B^BS (their (8) matches BS14; they cite BS14) — our
  reporting convention.  Transfer: flat-space bounds assumed to "continue to
  hold for the G_i evaluated on the cosmological <phi>" (MN), pointwise on J.
- Applicability: only (10) has an alpha-form, and only for Horndeski —
  POS is NOT applicable to the 8 alpha_H rows (GLPV operators are outside
  MN's amplitudes); on the 7 applicable rows it is applied exactly as MN
  apply it themselves (pointwise prior on the trajectory, their (12) usage).

Run:  PYTHONPATH=problems/P9/b/src:problems/P9/b uv run python -m p9b.positivity
Writes certificates/P9b-2.<ROW>.json for the 7 applicable rows.
"""
from __future__ import annotations

import json
from fractions import Fraction as F

from flint import arb, ctx

from .certify import (CERT_DIR, MAX_DEPTH, PREC_BITS, CertificationError,
                      _endpoint, ball, fr)
from .witnesses import A_C, GAP, REG, TARGET_POLY_FN, witness_desc

DATE = "2026-08-31"
POS_APPLICABLE = tuple(r for r in REG if "H" not in r)          # 7 rows
POS_NOT_APPLICABLE = {r: ("alpha_H != 0: beyond-Horndeski (GLPV) operator, outside the Horndeski amplitudes MN (14)-(15); POS makes no statement on this row") for r in REG if "H" in r}

POS_PIN = {
    "assumption_set": "POS (FORMULATION 2(b2); frozen)",
    "source": "Melville-Noller arXiv:1904.05874 (PRD 101, 021502(R)), eqs. (16)-(17); example-model alpha-form eq. (10); ar5iv fetch re-verified 2026-08-31",
    "inequality": "alpha_B^BS <= 2 alpha_T / (1 + alpha_T)  at every a in J",
    "side_condition": "1 + alpha_T > 0 (M^2, G4 > 0 in the MN derivation)",
    "convention": "alpha_B in BS14 convention (= MN's); S1 pin alpha_B^BS = -2 alpha_B^GLV",
    "conditionality": "flat-space bounds transferred to FLRW per MN's stated assumption; alpha-form applied pointwise to the trajectory (MN's own (12)-usage); not claimed to follow from UV physics on this background",
}


def pos_rhs(aT: F) -> F:
    """Exact RHS of MN (10); side condition 1 + alpha_T > 0 enforced."""
    assert 1 + aT > 0, "side condition 1 + alpha_T > 0 violated"
    return 2 * aT / (1 + aT)


# --- extended witness specs (b(a) = b1*a + b0; aT, aM constants or None) -----------
def spec_frozen(row: str) -> dict:
    w = REG[row]
    return dict(J=w["J"], b1=w["b1"], b0=F(0) if w["b1"] is not None else None,
                aT=w["aT"], aM=w["aM"], c_min=w["c_min"], D_min=w["D_min"])


# S4 alternative witnesses (ONE bounded attempt per violating row family,
# frozen after float scans, Arb-certified below): B-family -> braiding made
# negative (POS needs alpha_B <= 0 resp. <= -2/9); braidingless-T family ->
# alpha_T made positive (POS with alpha_B = 0 forces alpha_T >= 0).  Row KT
# has NO entry: POS + the exact row identity exclude the whole row
# (KT_IMPOSSIBILITY).  All D3 bounds hold: |b| <= 21/40, |aT| = 1/10, aM = 1/2.
_J = (F(3, 5), F(1))
ALT = {
    "KB":   dict(J=_J, b1=F(1, 2), b0=F(-21, 40), aT=None, aM=None, c_min=F(1, 10), D_min=F(1)),
    "KBT":  dict(J=_J, b1=F(0), b0=F(-1, 4), aT=F(-1, 10), aM=None, c_min=F(1, 20), D_min=F(109, 100)),
    "KBM":  dict(J=_J, b1=F(1, 2), b0=F(-21, 40), aT=None, aM=F(1, 2), c_min=F(3, 5), D_min=F(1)),
    "KBTM": dict(J=_J, b1=F(0), b0=F(-1, 4), aT=F(-1, 10), aM=F(1, 2), c_min=F(11, 20), D_min=F(109, 100)),
    "KTM":  dict(J=_J, b1=None, b0=None, aT=F(1, 10), aM=F(1, 2), c_min=F(1, 10), D_min=F(1)),
}

KT_IMPOSSIBILITY = {
    "statement": "Row {K,T}: {stable crossing} intersect POS is EMPTY — exact, whole row, family B (not just this witness).",
    "argument": [
        "POS with alpha_B == 0: 0 <= 2 aT/(1+aT); side condition 1+aT > 0 forces alpha_T >= 0 on J.",
        "Exact S2 row identity (derivation/rows.py; re-verified in tests/test_s4_positivity.py): alpha_K c_s^2 = fhat - 2 alpha_T, fhat = f/(MPl^2 H^2).",
        "D1 crossing with rho_DE > 0, H^2 > 0 gives fhat < 0 on the phantom side, so alpha_K c_s^2 = fhat - 2 alpha_T < 0 there.",
        "alpha_K > 0 (D3) forces c_s^2 < 0: D2(ii) fails. No bounded alpha_T >= 0 rescues this.",
    ],
    "level": "exact identity + stated sign analysis (level of P9b-0)",
}


def alt_desc(spec: dict) -> dict:
    """JSON-ready exact description of an alternative witness."""
    b = ("0" if spec["b1"] is None else
         f"{spec['b1']}*a + {spec['b0']}" if spec["b1"] != 0 else str(spec["b0"]))
    return {"alpha_K": "1", "alpha_B_BS": b,
            "alpha_T": str(spec["aT"]) if spec["aT"] is not None else "0",
            "alpha_M": str(spec["aM"]) if spec["aM"] is not None else "0",
            "alpha_H": "0",
            "M2_over_MPl2": f"a^({spec['aM']})" if spec["aM"] is not None else "1",
            "J": [str(spec["J"][0]), str(spec["J"][1])]}


def ext_eval(spec: dict, lo: F, hi: F) -> dict:
    """Row quantities as balls on [lo, hi]; mirrors certify.row_eval but for
    affine braiding b(a) = b1*a + b0 (b' = dB/dx = b1*a) and either-sign aT;
    alpha_H = 0 always (POS is only evaluated on Horndeski rows)."""
    A = ball(lo, hi)
    lnA = A.log()
    gde = ((fr(F(9, 5)) * lnA) + fr(F(12, 5)) * (1 - A)).exp()
    a3 = A**(-3)
    E2 = fr(F(3, 10)) * a3 + fr(F(7, 10)) * gde
    Om = fr(F(3, 10)) * a3 / E2
    Ode = fr(F(7, 10)) * gde / E2
    w1 = (4 * A - 3) / 5
    fhat = 3 * Ode * w1
    hd = -fr(F(3, 2)) * (Om + w1 * Ode)
    aT = fr(spec["aT"]) if spec["aT"] is not None else arb(0)
    aM = fr(spec["aM"]) if spec["aM"] is not None else arb(0)
    if spec["b1"] is None:
        b, db = arb(0), arb(0)
    else:
        b, db = fr(spec["b1"]) * A + fr(spec["b0"]), fr(spec["b1"]) * A
    m2 = (aM * lnA).exp() if spec["aM"] is not None else arb(1)
    rm = 3 * Om / m2
    target = TARGET_POLY_FN(b, db, aT, aM, arb(0), arb(0), hd, rm)
    D = 1 + fr(F(3, 2)) * b**2
    rhs = fr(pos_rhs(spec["aT"] if spec["aT"] is not None else F(0)))
    return dict(cs2=target / D, D=D, E2=E2, Ode=Ode, m2=m2, fhat=fhat, aB=b,
                qpos=rhs - b)


def sweep(spec: dict, key: str, lo: F, hi: F, bound: F, side: int,
          max_depth: int = MAX_DEPTH) -> dict:
    """certify.sweep on the extended evaluator: certify q > bound (side=+1)
    or q < bound (side=-1) on [lo, hi]; strict ball verdicts throughout."""
    bnd = fr(bound)
    stack, n_sub, depth_max, extreme = [(lo, hi, 0)], 0, 0, None
    while stack:
        a, b, d = stack.pop()
        v = ext_eval(spec, a, b)[key]
        if (v > bnd) if side > 0 else (v < bnd):
            n_sub += 1
            depth_max = max(depth_max, d)
            e = _endpoint(v, -side)
            extreme = e if extreme is None else (min if side > 0 else max)(extreme, e)
            continue
        if d >= max_depth:
            raise CertificationError(f"{key} not strict on [{a},{b}]")
        m = (a + b) / 2
        stack += [(a, m, d + 1), (m, b, d + 1)]
    return {"quantity": key, "bound": str(bound),
            "relation": ">" if side > 0 else "<",
            "interval": [str(lo), str(hi)], "n_subintervals": n_sub,
            "max_depth": depth_max,
            ("min_lower_endpoint" if side > 0 else "max_upper_endpoint"): extreme,
            "pass": True}


def pos_verdict(spec: dict) -> tuple[str, dict]:
    """Certified POS verdict for one witness spec on its J.  Returns
    (verdict, check): 'satisfies_exact' (alpha_B = alpha_T = 0: (10) reads
    0 <= 0, exact), 'satisfies_strict' (Arb: qpos > 0 on J), or
    'violates_everywhere' (Arb: qpos < 0 on J, i.e. the certified violating
    subinterval is ALL of J)."""
    a1, a2 = spec["J"]
    if spec["b1"] is None and spec["aT"] is None:
        return "satisfies_exact", {
            "exact": "alpha_B = alpha_T = 0 identically: (10) reads 0 <= 0 (weak inequality saturated at every a in J)",
            "pass": True}
    try:
        return "satisfies_strict", sweep(spec, "qpos", a1, a2, F(0), +1)
    except CertificationError:
        return "violates_everywhere", sweep(spec, "qpos", a1, a2, F(0), -1)


def certify_alt(row: str) -> dict:
    """Full Arb re-certification (stability D1-D3 + strict POS) of the
    alternative witness for a violating row; mirrors certify.certify_row."""
    s = ALT[row]
    a1, a2 = s["J"]
    assert a1 < A_C < a2
    ck = {"cs2_ge_cmin": sweep(s, "cs2", a1, a2, s["c_min"], +1)}
    if s["b1"] is None:
        ck["D_ge_Dmin"] = {"exact": "D = alpha_K = 1 (alpha_B = 0)",
                           "bound": str(s["D_min"]), "pass": True}
    else:
        ck["D_ge_Dmin"] = sweep(s, "D", a1, a2, s["D_min"], +1)
        ck["aB_ne_2"] = sweep(s, "aB", a1, a2, F(2), -1)
    cT2 = F(1) + (s["aT"] if s["aT"] is not None else F(0))
    assert cT2 >= 0
    ck["cT2_nonneg"] = {"exact": f"c_T^2 = 1 + alpha_T = {cT2}", "pass": True}
    if s["aM"] is not None:
        ck["M2_pos"] = sweep(s, "m2", a1, a2, F(0), +1)
    ck["Ode_pos"] = sweep(s, "Ode", a1, a2, F(0), +1)
    ck["E2_pos"] = sweep(s, "E2", a1, a2, F(0), +1)
    ck["fhat_neg_left"] = sweep(s, "fhat", a1, A_C - GAP, F(0), -1)
    ck["fhat_pos_right"] = sweep(s, "fhat", A_C + GAP, a2, F(0), +1)
    ck["pos_strict"] = sweep(s, "qpos", a1, a2, F(0), +1)
    assert all(c["pass"] for c in ck.values())
    return {"witness": alt_desc(s), "c_min": str(s["c_min"]), "D_min": str(s["D_min"]), "checks": ck,
            "note": "same background F*, same J, all D3 bounds hold; witness is joint (every row alpha nonzero on J)"}


def certify_pos_row(row: str) -> dict:
    """The P9b-2.<ROW> certificate (applicable rows only)."""
    assert row in POS_APPLICABLE, POS_NOT_APPLICABLE.get(row, row)
    s = spec_frozen(row)
    verdict, check = pos_verdict(s)
    aT = s["aT"] if s["aT"] is not None else F(0)
    cert = {
        "claim": f"P9b-2.{row}", "row": row, "conditional_on": "POS",
        "date": DATE, "prec_bits": PREC_BITS, "pos_pin": POS_PIN,
        "frozen_witness": witness_desc(row),
        "pos_rhs": f"2 aT/(1+aT) = {pos_rhs(aT)}",
        "side_condition": f"1 + alpha_T = {1 + aT} > 0: holds",
        "verdict": verdict, "pos_check": check,
    }
    if verdict == "violates_everywhere":
        cert["violating_subinterval"] = [str(s["J"][0]), str(s["J"][1])]
    if row == "KM":
        cert["status"] = ("(i) certified stable-crossing witness satisfying POS exists (the frozen witness; exact saturation; strict satisfaction unattainable in-row since alpha_B = alpha_T = 0 identically)")
    elif row == "KT":
        cert["alternative"] = {"found": False, "impossible": KT_IMPOSSIBILITY}
        cert["status"] = ("(ii+) frozen witness violates POS on all of J, and the whole row is excluded: no {K,T} witness satisfies POS and stable crossing (exact identity)")
    else:
        cert["alternative"] = dict(found=True, **certify_alt(row))
        cert["status"] = ("(i) frozen witness violates POS on all of J; POS-compatible alternative witness certified (stability + strict POS)")
    return cert


def main(cert_dir=CERT_DIR) -> dict:
    ctx.prec = PREC_BITS
    cert_dir.mkdir(exist_ok=True)
    out = {}
    for row in POS_APPLICABLE:
        cert = certify_pos_row(row)
        (cert_dir / f"P9b-2.{row}.json").write_text(
            json.dumps(cert, indent=1) + "\n")
        out[row] = cert
        q = cert["pos_check"]
        m = q.get("min_lower_endpoint", q.get("max_upper_endpoint", "exact 0"))
        alt = ("row excluded" if row == "KT" else "n/a" if row == "KM" else
               "alt certified")
        print(f"{row:5s} {cert['verdict']:20s} margin {m}  [{alt}]")
    for row, why in POS_NOT_APPLICABLE.items():
        print(f"{row:5s} POS not applicable ({why.split(':')[0]})")
    return out


if __name__ == "__main__":
    main()
