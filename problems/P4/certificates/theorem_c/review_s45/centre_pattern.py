"""S4-5 review, G1 centre half (exact): zero pattern of P_0 and G_0(kappa) of the
scaled 4D system at t = 0 (linscaled.full_system; T = 0 and theta z = 0), as
polynomials in (n, w, v).  Shows det(n P_0 - G_0(kappa)) is kappa-free."""
import json, sys
from p4.validated import linscaled
from p4.validated.linsys import LinSystem
S = linscaled.full_system(LinSystem())
gens = linscaled.CTX.gens()
names = [str(g) for g in gens]
n, w, v, T, dn, dw, dv = gens[1:]
zero = linscaled.CTX.constant(0)
at0 = lambda f: f.compose(gens[0], n, w, v, zero, zero, zero, zero)
d = S.d
out = {}
for nm, M in (("P", S.P), ("G", S.G)):
    for j, Mj in enumerate(M):
        pat = [[("0" if at0(Mj[r][c]).is_zero() else "*") for c in range(d)] for r in range(d)]
        out["%s_kappa^%d at t=0" % (nm, j)] = ["".join(row) for row in pat]
        if nm == "G" and j == 1:
            out["G_kappa^1 nonzero entries"] = [str(at0(Mj[r][c])) for r in range(d) for c in range(d) if not at0(Mj[r][c]).is_zero()]
        if nm == "P" and j == 0:
            out["P_0 rows 1,3"] = [[str(at0(Mj[r][c])) for c in range(d)] for r in (1, 3)]
print("gens:", names)
for k, val in out.items():
    print(k, val)
json.dump(out, open(sys.argv[1], "w"), indent=1)
