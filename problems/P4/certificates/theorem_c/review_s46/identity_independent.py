#!/usr/bin/env python3
"""S4-6 review: independent exact check of the centre nB identity.

Claim under review (tc_layer_centre.check_nb_identity, notes/s4-6-centre.md
section 3.2): for kappa^1..kappa^7 the y^8 and y^9 Taylor coefficients of
U6 = -nB vanish on every centre background A = 1 + y^2 Ahat, N = n/y,
W = y^2 w, V = y v, so nB's kappa^1..kappa^7 coefficients have Laurent
offset 10, not 8.

This script does not call tc_layer_centre.  It
  1. builds U6 a second way, as U1 / Dp^6 from the unreduced reduction pair
     (U1 = Bc*UN3c - Ac*dBc*D2, UN3c = dAc*D2 + SD*(Bc*Ae - Ac*Be): the
     numerator the S4-5 review verified in sympy, I3), by an exact fmpq_mpoly
     division with zero remainder, and checks it equals the formula
     build_pairs evaluates (Bcr*UN3r - Acr*EBc*D2r);
  2. substitutes truncated SYMBOLIC series for the background,
        n = n0 + n1 y + n2 y^2 + ...,  w = w0 + w1 y + ...,  v = v0 + v1 y + ...,
        Ahat = a0 + a1 y + ...,  A = 1 + y^2 Ahat,
     odd coefficients included, directly into every monomial of U6 (no layer
     bucketing), and collects the exact coefficient of y^k kappa^j for
     k <= 10 as a polynomial in all the background symbols;
  3. checks: for j = 1..7 the y^8 and y^9 coefficients are the zero
     polynomial (for every background, even or not) and the y^10 coefficient
     is not; for j = 0, 8, 9 nothing below y^10; no odd-parity monomial
     (e_N + e_V odd); the same scan for denU = dB (y^8 nonzero for 1..7);
  4. checks that every atom build_pairs reads (value, d-, dd- slots) involves
     only A, N, W, V, k (tc_energy.compile_mp keeps only those exponents).

Pure exact arithmetic.  Run from the repository root:
    PYTHONPATH=problems/P4/src python problems/P4/results/theorem_c/review_s46/identity_independent.py [out.json]
"""
import json
import os
import sys
import time

from flint import fmpq, fmpq_mpoly_ctx

HERE = os.path.dirname(os.path.abspath(__file__))
KTOP = 10          # collect y^k for k <= KTOP


def main(out_path):
    t0 = time.time()
    import p4.validated.tc_energy as te
    at = te.Atoms(full=True)
    mp = at.mp
    red = at.red
    rp = red.reduce_pair()
    Ac, Bc = rp["chi"][0], rp["chi"][1]
    Ae, Be = rp["eta"][0], rp["eta"][1]
    D2, SD, Dp, D1 = rp["D2"], red.SD, red.Dpoly, red.D1
    names = Dp.context().names()
    rec = dict(schema="review-s46-identity-independent/1", ctx_names=list(names))

    # ---- 1. U6 two ways
    dAc, dBc = D1(Ac), D1(Bc)
    UN3c = dAc * D2 + SD * (Bc * Ae - Ac * Be)
    U1 = Bc * UN3c - Ac * dBc * D2
    q, r = divmod(U1, Dp ** 6)
    rec["U1_divisible_by_Dp6"] = bool(r.is_zero())
    U6 = q
    EAc = 2 * mp["dDp"] * mp["Acr"] + Dp * mp["dAcr"]
    EBc = 2 * mp["dDp"] * mp["Bcr"] + Dp * mp["dBcr"]
    UN3r = EAc * mp["D2r"] + mp["SDr"] * (mp["Bcr"] * mp["Aer"] - mp["Acr"] * mp["Ber"])
    U6_formula = mp["Bcr"] * UN3r - mp["Acr"] * EBc * mp["D2r"]
    rec["U6_equals_build_pairs_formula"] = bool((U6 - U6_formula).is_zero())
    denU = Dp ** 3 * mp["SDr"] * mp["D2r"] ** 2 * mp["Bcr"]
    rec["U6_terms"] = len(U6)
    rec["denU_terms"] = len(denU)

    # ---- 4. variables of every atom build_pairs reads
    iA, iN, iW, iV, ik = (names.index(v) for v in ("A", "N", "W", "V", "k"))
    keep = {iA, iN, iW, iV, ik}
    used_atoms = ["SD", "dSD"]
    for nm in ("TN1r", "Dp", "SDr", "D2r", "Bcr", "Acr", "Aer", "Ber", "dDp", "dAcr", "dBcr"):
        used_atoms += [nm, "d" + nm, "dd" + nm]
    extra = {}
    for nm in used_atoms:
        if nm not in mp:
            extra[nm] = "MISSING"
            continue
        bad = sorted({names[i] for e, c in mp[nm].terms() for i, x in enumerate(e)
                      if int(x) and i not in keep})
        if bad:
            extra[nm] = bad
    rec["atoms_checked"] = sorted(set(used_atoms))
    rec["atoms_with_other_variables"] = extra

    # ---- 2./3. series substitution
    def scan(P, label):
        terms = []
        pmin = None
        odd = 0
        for e, c in P.terms():
            e = [int(x) for x in e]
            if any(e[i] for i in range(len(e)) if i not in keep):
                raise RuntimeError("%s: a monomial in another variable" % label)
            eA, eN, eW, eV, ek = e[iA], e[iN], e[iW], e[iV], e[ik]
            p = -eN + 2 * eW + eV
            pmin = p if pmin is None else min(pmin, p)
            if (eN + eV) % 2:
                odd += 1
            terms.append((eA, eN, eW, eV, ek, p, c))
        return terms, pmin, odd

    def layered(P, label, a_linear=False):
        """{kappa power: {y power: exact coefficient}} and the scan data, for
        P on the symbolic background (a_linear: A = 1 + alpha y + y^2 Ahat,
        a control that breaks the A = 1 + O(y^2) structure)."""
        terms, pmin, odd = scan(P, label)
        D = KTOP - pmin                     # series order needed for the lowest layer
        syms = []
        for s in ("n", "w", "v", "a"):
            syms += ["%s%d" % (s, i) for i in range(D + 1)]
        if a_linear:
            syms.append("alpha")
        R = fmpq_mpoly_ctx.get(tuple(syms))
        g = dict(zip(syms, R.gens()))
        one = R.constant(1)
        zero = R.constant(0)
        ser = {s: [g["%s%d" % (s, i)] for i in range(D + 1)] for s in ("n", "w", "v", "a")}
        # A = 1 + y^2 Ahat as a series (control: + alpha y)
        Aser = [one, g["alpha"] if a_linear else zero] + [ser["a"][i] for i in range(D - 1)]
        base = dict(n=ser["n"], w=ser["w"], v=ser["v"], A=Aser)

        def mul(x, y_, order):
            out = [zero] * (order + 1)
            for i in range(min(order, len(x) - 1) + 1):
                if x[i].is_zero():
                    continue
                for j in range(min(order - i, len(y_) - 1) + 1):
                    out[i + j] = out[i + j] + x[i] * y_[j]
            return out

        powc = {}

        def power(b, k, order):
            key = (b, k, order)
            if key in powc:
                return powc[key]
            if k == 0:
                r_ = [one] + [zero] * order
            else:
                h = power(b, k // 2, order)
                r_ = mul(h, h, order)
                if k % 2:
                    r_ = mul(r_, base[b], order)
            powc[key] = r_
            return r_

        acc = {}
        for eA, eN, eW, eV, ek, p, c in terms:
            if p > KTOP:
                continue
            o = KTOP - p
            s_ = power("n", eN, o)
            s_ = mul(s_, power("w", eW, o), o)
            s_ = mul(s_, power("v", eV, o), o)
            s_ = mul(s_, power("A", eA, o), o)
            row = acc.setdefault(ek, {})
            for i in range(o + 1):
                if not s_[i].is_zero():
                    row[p + i] = row.get(p + i, zero) + s_[i] * c
        return acc, pmin, odd, D, syms, zero

    res = {}
    for label, P in (("U6", U6), ("denU", denU)):
        acc, pmin, odd, D, syms, zero = layered(P, label)
        out = dict(min_layer=pmin, odd_parity_terms=odd, series_order=D,
                   symbols=syms, kappa={})
        for j in sorted(acc):
            row = acc[j]
            ks = sorted(row)
            nz = [k for k in ks if not row[k].is_zero()]
            # symbols that appear in the y^10 coefficient (if any)
            vars10 = None
            if KTOP in row and not row[KTOP].is_zero():
                vars10 = sorted({syms[i] for e, cc in row[KTOP].terms()
                                 for i, x in enumerate(e) if int(x)})
            out["kappa"][j] = dict(
                y_powers_with_terms=ks,
                y_powers_nonzero=nz,
                lowest_nonzero=(nz[0] if nz else None),
                y8_zero=(row.get(8, zero).is_zero()),
                y9_zero=(row.get(9, zero).is_zero()),
                y10_terms=(len(row[KTOP]) if KTOP in row else 0),
                y10_symbols=vars10)
        res[label] = out
    rec["scan"] = res

    u6 = res["U6"]["kappa"]
    du = res["denU"]["kappa"]
    checks = dict(
        U6_built_two_ways=rec["U1_divisible_by_Dp6"] and rec["U6_equals_build_pairs_formula"],
        no_other_variables=not extra,
        no_odd_parity_U6=res["U6"]["odd_parity_terms"] == 0,
        kappa1to7_y8_y9_vanish=all(u6[j]["y8_zero"] and u6[j]["y9_zero"] for j in range(1, 8)),
        kappa1to7_y10_nonzero=all(u6[j]["lowest_nonzero"] == 10 for j in range(1, 8)),
        kappa1to7_lowest_term_layer_8=all(min(u6[j]["y_powers_with_terms"]) == 8
                                          for j in range(1, 8)),
        other_kappa_start_at_10=all((u6[j]["lowest_nonzero"] or 99) >= 10
                                    for j in u6 if j not in range(1, 8)),
        denU_kappa1to7_y8_nonzero=all(du[j]["lowest_nonzero"] == 8 for j in range(1, 8)),
    )
    # ---- controls: the check must fail where the claim would be false
    ctl = {}
    # (a) A = 1 + alpha y + y^2 Ahat: the y^9 coefficient picks up
    #     alpha * dL8/dA, so the A = 1 + O(y^2) structure is load-bearing
    acc, *_ = layered(U6, "U6 (control a)", a_linear=True)
    ctl["a_linear_A_y9_nonzero_kappa"] = [j for j in range(1, 8)
                                          if not acc[j].get(9, None) is None
                                          and not acc[j][9].is_zero()]
    # (b) a mutated U6: add one layer-8 monomial N^2 W^5 kappa^3 (coefficient 1)
    ctxU = U6.context()
    gens = dict(zip(names, ctxU.gens()))
    mut = U6 + gens["N"] ** 2 * gens["W"] ** 5 * gens["k"] ** 3
    acc, *_ = layered(mut, "U6 (control b)")
    ctl["mutated_y8_nonzero_kappa"] = [j for j in range(1, 8)
                                       if 8 in acc[j] and not acc[j][8].is_zero()]
    # (c) a mutated U6 with one odd-parity layer-9 monomial N W^5 kappa^2
    mut = U6 + gens["N"] * gens["W"] ** 5 * gens["k"] ** 2
    acc, _, odd_c, *_ = layered(mut, "U6 (control c)")
    ctl["mutated_odd_parity_terms"] = odd_c
    ctl["mutated_y9_nonzero_kappa"] = [j for j in range(1, 8)
                                       if 9 in acc[j] and not acc[j][9].is_zero()]
    rec["controls"] = ctl
    checks["controls_fail_as_expected"] = (
        bool(ctl["a_linear_A_y9_nonzero_kappa"]) and ctl["mutated_y8_nonzero_kappa"] == [3]
        and ctl["mutated_odd_parity_terms"] == 1 and ctl["mutated_y9_nonzero_kappa"] == [2])
    rec["checks"] = checks
    rec["ok"] = all(checks.values())
    rec["note"] = ("y^9 vanishes for every background with A = 1 + O(y^2), not only even "
                   "ones: the y^9 coefficient of the layer-8 part is the derivative of the "
                   "layer-8 polynomial in the direction (n1, w1, v1), and that polynomial is "
                   "identically 0; there is no layer-9 monomial (no odd parity).")
    rec["seconds"] = time.time() - t0
    with open(out_path, "w") as fh:
        json.dump(rec, fh, indent=1, sort_keys=True, default=str)
    print(json.dumps(dict(ok=rec["ok"], checks=checks, seconds=rec["seconds"]), indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "identity_independent.json"))
