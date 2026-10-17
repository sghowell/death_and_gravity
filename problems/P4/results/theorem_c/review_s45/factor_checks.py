"""S4-5 review: exact factor structure behind S3/S12 (fmpq_mpoly, via tc_energy.Atoms):
every polynomial the strip endgame needs nonzero divides dT = Dp SDr D2r Bcr up to Dp powers."""
import json, sys
from p4.validated.tc_energy import Atoms
at = Atoms(full=True)                       # runs the import-time _dpow checks (raise on failure)
red, mp = at.red, at.mp
Dp, SDr, D2r, Bcr = mp["Dp"], mp["SDr"], mp["D2r"], mp["Bcr"]
S, W, N, k = red.S, red.W, red.N, red.k
eDp, r = divmod(red.e, Dp)
out = {}
out["e = Dp * (e/Dp) exact"] = r.is_zero()
out["e/Dp == 4 S W (Q1 + k N)"] = (eDp - 4 * S * W * (red.Q[1] + k * N)).is_zero()
c16 = 16 * S**4 * W**2
out["D2r == 16 S^4 W^2 * dn * (e/Dp)"] = (D2r - c16 * red.dn * eDp).is_zero()
out["SDr == 4 S^2 W  (kappa-free)"] = (SDr - 4 * S**2 * W).is_zero()
out["SD == Dp * SDr"] = (red.SD - Dp * SDr).is_zero()
q, r = divmod(red.Q[1], N)
out["N divides Q1 (so Q1 + kN = N (F_N + k))"] = r.is_zero()
rp = red.reduce_pair()
out["D2 == Dp^3 D2r"] = (rp["D2"] - Dp**3 * D2r).is_zero()
out["Bc == Dp^2 Bcr"] = (rp["chi"][1] - Dp**2 * Bcr).is_zero()
out["S Delta~ == det P  (P invertible iff SD != 0)"] = (red.SD - red.S * (red.P[2][2]*red.P[3][3] - red.P[2][3]*red.P[3][2])).is_zero()
out["kappa-degrees (Dp, SDr, D2r, Bcr)"] = [int(f.degrees()[9]) for f in (Dp, SDr, D2r, Bcr)]
out["import-time _dpow checks and U1 = Dp^6 U6"] = "passed (Atoms(full=True) constructed)"
for kk, v in out.items():
    print(f"{kk}: {v}")
json.dump(out, open(sys.argv[1], "w"), indent=1)
