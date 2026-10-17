"""S4-5 review, G1: the order-0 sonic solve M0(kappa) of linsonic4._run4 and the
Frobenius exponent sigma(kappa), on the certified V0 box (balls).

det M0(kappa) is affine in kappa (only the ell.G row depends on kappa).  We
enclose its two coefficients, its root kappa*, and compare with sigma(kappa).
Also: D, E of the level equations, the rank-1 structure of E (the basis of
det(nD + E) = det D n^3 (n - sigma)), and Re sigma on Re kappa >= 0.

Run from a scratch copy:  PYTHONPATH=<copy>/problems/P4/src python g1_order0.py
"""
import json
import sys

from flint import acb, acb_mat, arb

from p4.validated import linmatch, modecount
from p4.validated.arbseries import precision
from p4.validated.linsonic4 import _Structure4, polys4, series_parts4, sigma
from p4.validated.linsys import LinSystem
from p4.validated.shootsys import eval_box

out = {}
with precision(256):
    bg = linmatch.box_background(modecount.V0_EC, modecount.W_V0, K=8)
    L = LinSystem()
    polys = polys4(L)
    Pser, Gser = series_parts4(polys, bg.series())
    ell, row = bg.point["ell"], bg.point["row"]
    rows = [[(arb(1), 0, -1)], [(arb(1), 1, -1)], [(arb(1), row, -1)],
            [(ell[0], 2, 0), (ell[1], 3, 0)]]
    u0 = bg.balls[0]
    base = acb(0)
    st = _Structure4(Pser, Gser, base)
    args = [arb(0)] + list(u0)
    A0, S0 = u0[0], eval_box(L.S, args)
    CW, CV = eval_box(L.dC[2], args), eval_box(L.dC[3], args)
    lG = lambda dg, l: ell[0] * st.G[dg][0][2, l] + ell[1] * st.G[dg][0][3, l]   # noqa: E731
    # M0(kappa) = M00 + kappa M01 (base = 0), exactly the matrix of _run4
    M00 = acb_mat([[lG(0, 2), lG(0, 3)], [A0 * CW, A0 * CV]])
    M01 = acb_mat([[lG(1, 2), lG(1, 3)], [acb(0), acb(0)]])
    d0 = M00.det()
    d1 = M01[0, 0] * M00[1, 1] - M01[0, 1] * M00[1, 0]     # det is affine: d0 + d1 kappa
    kstar = -d0 / d1
    D, E = st.DE(rows)
    s0 = sigma(D, E[0])                                      # sigma(0)
    s1 = sigma(D, E[1]) - acb(0)                             # d sigma / d kappa (E[1] = dE/dkappa)
    # sigma(kappa) = -tr(D^-1 (E0 + kappa E1)) = s0 + kappa * (-tr(D^-1 E1))
    s1 = -sum(((D.inv() * E[1])[i, i] for i in range(4)), acb(0))
    sroot = -s0 / s1
    # proportionality test: det M0 = lambda sigma  <=>  d0 s1 - d1 s0 = 0
    prop = d0 * s1 - d1 * s0
    # rank-1 structure of E: rows 0..2 of E0, E1 exactly zero
    rank1 = all(E[k][i, l] == 0 for k in (0, 1) for i in range(3) for l in range(4))
    # the right null vector of the fluid block P_fl(u0) and C~ . r
    P0 = st.P[0]
    # fluid block rows/cols 2,3
    Pfl = acb_mat([[P0[2, 2], P0[2, 3]], [P0[3, 2], P0[3, 3]]])
    rnull = (P0[2, 3] * -1, P0[2, 2])                          # (-P23, P22) annihilated by row 2
    resid_row3 = P0[3, 2] * rnull[0] + P0[3, 3] * rnull[1]
    Cr = CW * rnull[0] + CV * rnull[1]
    for k, v in dict(A0=A0, S0=S0, CW=CW, CV=CV, d0=d0, d1=d1, kstar=kstar, sigma0=s0,
                     dsigma=s1, sigma_root=sroot, prop_d0s1_minus_d1s0=prop,
                     detPfl=Pfl.det(), Pfl_null_resid=resid_row3, Ctilde_dot_rnull=Cr,
                     detD=D.det()).items():
        out[k] = v.str(12, radius=True) if hasattr(v, "str") else str(v)
        print(f"{k:24s} {out[k]}")
    print("E rank-1 (rows 0-2 exactly 0):", rank1)
    out["E_rank1"] = rank1
    print("kstar real part < 0:", (kstar.real < 0), " |Im kstar| ball:", kstar.imag.str(5, radius=True))
    out["kstar_re_negative"] = bool(kstar.real < 0)
json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "g1_order0.json", "w"), indent=1)
