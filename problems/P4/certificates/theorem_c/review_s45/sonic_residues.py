"""S4-5 review: fixed-kappa residues at the sonic point x = 0 (balls, certified
V0 box), for the boundary-term lemma S8 of the strip and s4-energy section 2:

  R_T(kappa) := lim x T = TN1r(u0; k) / (E_D SDr(u0) D2r(u0; k) Bcr(u0; k)),
  R_B(kappa) := lim x^2 B0 = -(dU6/dx)(0; k) / (E_D^3 SDr D2r^2 Bcr)(u0; k),
  Res tt     := lim x tt = TN1r_6(u0) / (E_D SDr D2r_2 Bcr_3)(u0),

with E_D = D'(0) and dU6/dx(0) = grad U6(u0) . u'(0) (U6(u0) = 0: E5).
Claims: R_T = 1 + sigma(kappa), R_B = 1 + sigma(kappa), Res tt = -(1 - sigma0);
hence Res A1 = R_T... (sign conventions: A0 = -T, B0 = -U6/denU).
Checked as kappa-polynomial identities with ball coefficients (each
coefficient of the difference must contain 0), and at sample strip kappa.
"""
import json
import sys

from flint import acb, arb

from p4.validated.arbseries import precision
from p4.validated.tc_energy import Atoms, kval, pmul, pscale, sonic_expansion_cached, u_prime_balls, u_to_pw

out = {}
with precision(256):
    at = Atoms(full=True)
    ex = sonic_expansion_cached(256)
    x0 = arb(0)
    u = ex.eval(x0)
    up = u_prime_balls(ex, x0, 256)
    pw = u_to_pw(u, at.dmax)
    ED = sum((kval(at.cmp["gDp_" + nm], pw, 0)[0] * d for nm, d in
              (("A", up[0]), ("N", up[1]), ("W", up[2]), ("V", up[3]))), arb(0))
    SDr = kval(at.cmp["SDr"], pw, 0)[0]
    D2r, Bcr = kval(at.cmp["D2r"], pw), kval(at.cmp["Bcr"], pw)
    TN1r = kval(at.cmp["TN1r"], pw)
    Dp0 = kval(at.cmp["Dp"], pw, 0)[0]
    # sigma(kappa) = -sigma0 - (1 - sigma0) kappa from linsonic4 (balls, G1 script)
    s0 = arb("0.0901059705079") + arb(0, 6e-14)
    # R_T identity: TN1r(u0;k) - (1+sigma(k)) E_D SDr D2r Bcr  == 0 coefficientwise
    den = pscale(pmul(D2r, Bcr), ED * SDr)
    onepsig = [1 - s0, -(1 - s0)]            # 1 + sigma(k) = (1 - s0)(1 - k)
    rhs = pmul(onepsig, den)
    diffT = [a - b for a, b in zip(TN1r + [arb(0)] * (len(rhs) - len(TN1r)), rhs)]
    out["Dp(u0)"] = Dp0.str(5, radius=True)
    out["E_D"] = ED.str(12, radius=True)
    out["R_T_identity_coeff_diffs"] = [d.str(3, radius=True) for d in diffT]
    out["R_T_identity_holds_in_balls"] = all(d.contains(0) for d in diffT)
    out["R_T_identity_rel_coeff_radius"] = max(float(d.rad()) / max(1e-300, float(abs(r).upper()))
                                               for d, r in zip(diffT, rhs) if not r.contains(0))
    # tt residue
    out["Res_tt"] = (TN1r[6] / (ED * SDr * D2r[2] * Bcr[3])).str(12, radius=True)
    out["minus_1_minus_sigma0"] = (-(1 - s0)).str(12, radius=True)
    # U6 and its x-derivative at the sonic point
    mp = at.mp
    EAc = 2 * mp["dDp"] * mp["Acr"] + mp["Dp"] * mp["dAcr"]
    EBc = 2 * mp["dDp"] * mp["Bcr"] + mp["Dp"] * mp["dBcr"]
    UN3r = EAc * mp["D2r"] + mp["SDr"] * (mp["Bcr"] * mp["Aer"] - mp["Acr"] * mp["Ber"])
    U6 = mp["Bcr"] * UN3r - mp["Acr"] * EBc * mp["D2r"]
    from p4.validated.tc_energy import compile_mp
    cU6 = compile_mp(U6)
    gU6 = [compile_mp(U6.derivative(i)) for i in (1, 2, 3, 4)]
    dmax = [max(a, b) for a, b in zip(at.dmax, cU6[1])]
    pw2 = u_to_pw(u, [d + 1 for d in dmax])
    U6v = kval(cU6, pw2)
    out["U6(u0) coefficient balls contain 0"] = all(v.contains(0) for v in U6v)
    dU6 = [arb(0)] * len(U6v)
    for g, d in zip(gU6, up):
        gv = kval(g, pw2, len(U6v) - 1)
        for j in range(len(dU6)):
            dU6[j] += gv[j] * d
    denB = pscale(pmul(pmul(D2r, D2r), Bcr), ED ** 3 * SDr)
    rhsB = pmul(onepsig, denB)
    lhsB = [-x for x in dU6]
    diffB = [a - b for a, b in zip(lhsB + [arb(0)] * (len(rhsB) - len(lhsB)), rhsB + [arb(0)] * (len(lhsB) - len(rhsB)))]
    out["R_B_identity_coeff_diffs"] = [d.str(3, radius=True) for d in diffB]
    out["R_B_identity_holds_in_balls"] = all(d.contains(0) for d in diffB)
    # sample strip kappa
    samples = {}
    for k in (acb(0, 14), acb(3, 20), acb(6, 14), acb(0.5, 1000), acb(6, 1e6)):
        def ev(p):
            acc = acb(0)
            for c_ in reversed(p):
                acc = acc * k + c_
            return acc
        RT = ev(TN1r) / ev(den)
        RB = ev(lhsB) / ev(denB)
        sig = -s0 - (1 - s0) * k
        A1res = 2 - RT                       # Res A1 = Res A0 - Res (log B0)' = -R_T + 2
        rho_res = A1res.real - arb(-1) / 32 * (-(1 - s0))
        samples[str(k)] = dict(R_T=RT.str(8, radius=True), one_plus_sigma=(1 + sig).str(8, radius=True),
                               R_B=RB.str(8, radius=True), gamma=rho_res.str(8, radius=True))
    out["samples"] = samples
    gam_min = 1 + s0 - (1 - s0) / 32
    out["gamma_min_eps0"] = gam_min.str(10, radius=True)
for k, v in out.items():
    print(k, ":", v)
json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "sonic_residues.json", "w"), indent=1)
