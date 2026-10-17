"""S4-5 review (pre-flight for the centre): kappa-roots of the lowest centre
layer (t^10) of U6 at the certified centre values (n, w, v)(t=0), i.e. the
kappa at which B0 = O(t^4) instead of t^2 at the centre.  A root in the strip
{0 <= Re k <= 6, |Im k| >= 14} would put a zero of B0 near (t = 0, k) and break
the commutation there.  Float diagnostic (numpy roots of a degree-9 polynomial
with ball-evaluated coefficients)."""
import json, sys
import numpy as np
from flint import arb, fmpq_mpoly_ctx
from p4.validated.arbseries import precision
from p4.validated import ec_centre
from p4.validated.tc_energy import Atoms
at = Atoms(full=True); mp = at.mp
EAc = 2 * mp["dDp"] * mp["Acr"] + mp["Dp"] * mp["dAcr"]
EBc = 2 * mp["dDp"] * mp["Bcr"] + mp["Dp"] * mp["dBcr"]
UN3r = EAc * mp["D2r"] + mp["SDr"] * (mp["Bcr"] * mp["Aer"] - mp["Acr"] * mp["Ber"])
U6 = mp["Bcr"] * UN3r - mp["Acr"] * EBc * mp["D2r"]
ctx3 = fmpq_mpoly_ctx.get(("n", "w", "v", "k"))
lay = {}
for e, cf in U6.terms():
    eN, eW, eV, ek = int(e[2]), int(e[3]), int(e[4]), int(e[9])
    p = -eN + 2 * eW + eV
    lay.setdefault(p, ctx3.constant(0))
    lay[p] += ctx3.from_dict({(eN, eW, eV, ek): cf})
lay = {p: l for p, l in lay.items() if not l.is_zero()}
pmin = min(lay)
terms = [(int(e[0]), int(e[1]), int(e[2]), int(e[3]), c) for e, c in lay[pmin].terms()]
with precision(256):
    ce = ec_centre.certified_centre_expansion(256, 30)
    ser = ce.series()                        # (n, w, v) series in t
    n0, w0, v0 = (s[0] for s in ser[:3])
    print("centre values n, w, v:", n0.str(10), w0.str(10), v0.str(10), " 2nv+1 =", (2 * n0 * v0 + 1).str(5, radius=True))
    coef = [arb(0)] * 10
    for eN, eW, eV, ek, cf in terms:
        coef[ek] += arb(cf.p) / arb(cf.q) * n0 ** eN * w0 ** eW * v0 ** eV
cf = [float(c.mid()) for c in coef]
roots = np.roots(cf[::-1])
print("lowest layer t^%d; kappa coefficients:" % pmin, ["%.6g" % c for c in cf])
print("roots:", ["%.6g%+.6gi" % (r.real, r.imag) for r in roots])
strip = [r for r in roots if 0 <= r.real <= 6 and abs(r.imag) >= 14]
print("roots in the strip:", strip)
json.dump(dict(pmin=pmin, coef=cf, roots=[[r.real, r.imag] for r in roots], in_strip=[[r.real, r.imag] for r in strip]),
          open(sys.argv[1], "w"), indent=1)
