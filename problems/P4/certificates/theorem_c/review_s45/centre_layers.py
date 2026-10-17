"""S4-5 review: exact centre-layer structure of A0 = -TN1r/denT and
B0 = -U6/denU (tc_commutator.build_pairs) at the regular centre, in the
grading of E2 (tc_commutator.check_layer_identity): t-weight
p = -eN + 2 eW + eV, A -> 1 (A - 1 = O(t^2)).

Claims checked (exact fmpq_mpoly):
  C-T  the lowest layer of TN1r + 3 denT lies strictly above that of denT,
       so T = TN1r/denT -> -3 and A0 -> 3 for every kappa at which the
       lowest layer of denT is nonzero at the centre values;
  C-B  the lowest layer of U6 lies strictly above that of denU, so B0 -> 0;
       k = val(U6) - val(denU) is reported (B0 = O(t^k)).
  C-tt the kappa-top coefficient: tt = TN1r_6/denT_5 -> 0 (positive t-weight).
Each test is done first as a polynomial identity in (n, w, v, kappa); if a
layer is not identically zero it is reduced modulo 2nv + 1 (the centre
fixed-point relation) and reported.
"""
import json
import sys
import time

from flint import fmpq, fmpq_mpoly_ctx

from p4.validated.tc_energy import Atoms

t0 = time.time()
at = Atoms(full=True)
mp = at.mp
print("Atoms built in %.1f s" % (time.time() - t0))
EAc = 2 * mp["dDp"] * mp["Acr"] + mp["Dp"] * mp["dAcr"]
EBc = 2 * mp["dDp"] * mp["Bcr"] + mp["Dp"] * mp["dBcr"]
UN3r = EAc * mp["D2r"] + mp["SDr"] * (mp["Bcr"] * mp["Aer"] - mp["Acr"] * mp["Ber"])
U6 = mp["Bcr"] * UN3r - mp["Acr"] * EBc * mp["D2r"]
denT = mp["Dp"] * mp["SDr"] * mp["D2r"] * mp["Bcr"]
denU = mp["Dp"] ** 3 * mp["SDr"] * mp["D2r"] ** 2 * mp["Bcr"]
TN1r = mp["TN1r"]

ctx3 = fmpq_mpoly_ctx.get(("n", "w", "v", "k"))
n_, w_, v_, k_ = ctx3.gens()


def layers(f):
    lay = {}
    for e, cf in f.terms():
        eN, eW, eV, ek = int(e[2]), int(e[3]), int(e[4]), int(e[9])
        p = -eN + 2 * eW + eV
        lay.setdefault(p, ctx3.constant(0))
        lay[p] += ctx3.from_dict({(eN, eW, eV, ek): cf})
    return {p: l for p, l in lay.items() if not l.is_zero()}


def kcoef(f, j):
    ctx = at.red.ctx
    out = ctx.constant(0)
    for e, cf in f.terms():
        if int(e[9]) == j:
            ee = list(map(int, e))
            ee[9] = 0
            out += ctx.from_dict({tuple(ee): cf})
    return out


def mod_fixed_point(l):
    """Reduce a layer polynomial modulo 2nv + 1 by v -> -1/(2n); returns the
    numerator after clearing powers of n (zero iff l in the ideal)."""
    dv = l.degrees()[2]
    num = ctx3.constant(0)
    for e, cf in l.terms():
        a, b, c, d = map(int, e)
        # n^a w^b v^c k^d with v = -1/(2n): n^(a-c) (-1/2)^c, times n^dv to clear
        num += ctx3.from_dict({(a - c + dv, b, 0, d): cf * fmpq(-1, 2) ** c})
    return num


out = {}
LT, LD = layers(TN1r), layers(denT)
pT, pD = min(LT), min(LD)
Z = TN1r + 3 * denT
LZ = layers(Z)
pZ = min(LZ) if LZ else None
out["T"] = dict(val_TN1r=pT, val_denT=pD, val_TN1r_plus_3denT=pZ,
                limit_minus3_exact=(pZ is None or pZ > pD),
                lowest_layer_mod_2nv1_zero=(pZ is not None and pZ <= pD and mod_fixed_point(LZ[pZ]).is_zero()))
print("T:", out["T"])
LU, LDU = layers(U6), layers(denU)
pU, pDU = min(LU), min(LDU)
out["B0"] = dict(val_U6=pU, val_denU=pDU, k=pU - pDU, B0_to_0=(pU > pDU),
                 U6_lowest_layer_mod_2nv1_zero=bool(mod_fixed_point(LU[pU]).is_zero()))
print("B0:", out["B0"])
# kappa-structure of the lowest layers (which kappa powers occur)
out["denT_lowest_layer_kappa_degrees"] = sorted({int(e[3]) for e, _ in LD[pD].terms()})
out["denU_lowest_layer_kappa_degrees"] = sorted({int(e[3]) for e, _ in LDU[pDU].terms()})
out["U6_lowest_layer_kappa_degrees"] = sorted({int(e[3]) for e, _ in LU[pU].terms()})
# parity: odd layers of Z and U6
out["Z_T_odd_layers"] = sorted(p for p in LZ if p % 2)
out["U6_odd_layers"] = sorted(p for p in LU if p % 2)
out["denT_lowest_layer"] = str(LD[pD])[:400]
out["denU_lowest_layer"] = str(LDU[pDU])[:400]
out["U6_lowest_layer"] = str(LU[pU])[:600]
# tt = TN1r_6 / denT_5
t6, d5 = kcoef(TN1r, 6), kcoef(denT, 5)
L6, L5 = layers(t6), layers(d5)
out["tt"] = dict(val_TN1r_6=min(L6), val_denT_5=min(L5), tt_to_0=min(L6) > min(L5))
print("tt:", out["tt"])
for k in ("denT_lowest_layer_kappa_degrees", "denU_lowest_layer_kappa_degrees",
          "U6_lowest_layer_kappa_degrees", "Z_T_odd_layers", "U6_odd_layers"):
    print(k, out[k])
print("denT lowest layer:", out["denT_lowest_layer"])
print("denU lowest layer:", out["denU_lowest_layer"])
print("U6 lowest layer:", out["U6_lowest_layer"])
out["time_s"] = time.time() - t0
json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "centre_layers.json", "w"), indent=1)
