#!/usr/bin/env python3
"""S4-6 review: symbolic checks of the graded evaluator's algebra (sympy; no
p4 code imported).

 1. The assembly T(lambda) of tc_layer_centre.sigma_graded (transcribed term
    by term) equals lambda^2 S / tau^2 for a FREE real lambda, with S built
    from the definitions of notes/s4-5-assembly.md section 4.1 (S4-S6) and
    the blocks V = lam (A0 + kappa tt), Y = lam^2 (B0 - u2 kappa^2)/kappa,
    V' = lam^2 (A0' + kappa tt'), Lh = lam (log B0)', Lp = lam^2 (log B0)'',
    t = lam tt, tp = lam^2 tt', u2h = lam u2 (tc_layer.LayerData's meaning
    of the ratios).  Controls: a sign flip in one term gives a nonzero
    residual.
 2. The regime closed forms of tc_layer_centre.Grade:
      L (g = y, sigma = 1/s, tau = s/y): 1/tau = g sigma, z = i s + g eps,
        zeta = 1/z = sigma (-i + y eps sigma)/(1 + (y eps sigma)^2),
        w = 1/kappa = y zeta, lam = 1/g gives Sigma = S/s^2, and
        norm = 1/(1 + sigma^2) = s^2/(1 + s^2);
      S (g = r = 1/tau): y = r s, z = i s + r s eps,
        w = r (r eps - i)/(1 + r^2 eps^2), lam = 1/r gives Sigma = S,
        norm = 1/(1 + s^2);
    so both give Sigma_c = S/(1 + s^2), whose sign is that of S.
 3. The graded inverse: 1/(1 + u) = sum_{k=0}^{K+1} (-u)^k + (-u)^{K+2}/(1 + u),
    K = KMAX = 4, and the fold g^a = g^{K+1} g^{a-K-1}.
"""
import json
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))


def check_T():
    lam, eps, tau = sp.symbols("lam eps tau", positive=True)
    c = sp.Rational(-1, 32)
    a0r, a0i, a0pr, a0pi = sp.symbols("a0r a0i a0pr a0pi", real=True)
    b0r, b0i = sp.symbols("b0r b0i", real=True)
    l1r, l1i, l2r, l2i = sp.symbols("l1r l1i l2r l2i", real=True)
    tt, ttp, u2 = sp.symbols("tt ttp u2", real=True)
    I = sp.I
    kap = eps + I * tau
    A0, A0p = a0r + I * a0i, a0pr + I * a0pi
    B0 = b0r + I * b0i
    L1, L2 = l1r + I * l1i, l2r + I * l2i
    # --- S from the definitions (assembly S4-S6)
    A1, A1p = A0 - L1, A0p - L2
    B1 = B0 + A0p - A0 * L1
    nu = sp.im(A1) + tau * tt
    nup = sp.im(A1p) + tau * ttp
    rho = sp.re(A1) - c * tt
    f = (c + tau ** 2 / c) * tt - (tau / c) * nu
    fp = (c + tau ** 2 / c) * ttp - (tau / c) * nup
    K = sp.re(B1) - (tau / c) * sp.im(B1) - (fp + rho * f) / 2
    S = K + nu ** 2 / 4
    # --- the blocks as the code's ratios mean them
    V = lam * (A0 + kap * tt)
    Y = lam ** 2 * (B0 - u2 * kap ** 2) / kap
    Vp = lam ** 2 * (A0p + kap * ttp)
    Lh = lam * L1
    Lp = lam ** 2 * L2
    t, tp, u2h = lam * tt, lam ** 2 * ttp, lam * u2
    it = 1 / tau
    E = eps
    # --- sigma_graded, transcribed
    W, Wp = V - Lh, Vp - Lp
    G = Y + t * Lh - tp
    H = Vp - V * Lh
    nuG, nupG = sp.im(W), sp.im(Wp)
    rhoG = sp.re(W) - (E + c) * t
    cf = c * it * it + 1 / c

    def T(sign_flip=False):
        reB1 = u2h * lam * (E * E * it * it - 1) + E * it * it * sp.re(G) \
            - it * sp.im(G) + it * it * sp.re(H)
        imB1 = -(2 * E * u2h * lam + sp.re(G) + E * it * sp.im(G)
                 + (-1 if sign_flip else 1) * it * sp.im(H)) * (1 / c)
        fpT = -cf * tp * sp.Rational(1, 2) + (it * (1 / c)) * nupG * sp.Rational(1, 2)
        rhof = -rhoG * (cf * t - (it * (1 / c)) * nuG) * sp.Rational(1, 2)
        nu2 = nuG * nuG * it * it * sp.Rational(1, 4)
        return reB1 + imB1 + fpT + rhof + nu2

    res = sp.simplify(sp.expand(T() - lam ** 2 * S / tau ** 2))
    ctl = sp.simplify(sp.expand(T(sign_flip=True) - lam ** 2 * S / tau ** 2))
    return dict(residual=str(res), control_residual_nonzero=ctl != 0)


def check_regimes():
    y, s, eps, r = sp.symbols("y s eps r", positive=True)
    I = sp.I
    out = {}
    # regime L
    tau = s / y
    kap = eps + I * tau
    sig = 1 / s
    z = y * kap
    zeta = sig * (-I + y * eps * sig) / (1 + (y * eps * sig) ** 2)
    out["L_it"] = sp.simplify(1 / tau - y * sig) == 0
    out["L_z"] = sp.simplify(z - (I * s + y * eps)) == 0
    out["L_zeta"] = sp.simplify(1 / z - zeta) == 0
    out["L_w"] = sp.simplify(1 / kap - y * zeta) == 0
    Ssym = sp.symbols("S", real=True)
    lam = 1 / y
    Sigma = Ssym * lam ** 2 / tau ** 2
    out["L_Sigma_is_S_over_s2"] = sp.simplify(Sigma - Ssym / s ** 2) == 0
    out["L_norm"] = sp.simplify(1 / (1 + sig ** 2) - s ** 2 / (1 + s ** 2)) == 0
    out["L_Sigma_c"] = sp.simplify(Sigma / (1 + sig ** 2) - Ssym / (1 + s ** 2)) == 0
    # regime S
    tau = 1 / r
    yS = r * s
    kap = eps + I * tau
    out["S_y"] = sp.simplify(yS * tau - s) == 0
    out["S_z"] = sp.simplify(yS * kap - (I * s + r * s * eps)) == 0
    out["S_w"] = sp.simplify(1 / kap - r * (r * eps - I) / (1 + r ** 2 * eps ** 2)) == 0
    Sigma = Ssym * (1 / r) ** 2 / tau ** 2
    out["S_Sigma_is_S"] = sp.simplify(Sigma - Ssym) == 0
    out["S_Sigma_c"] = sp.simplify(Sigma / (1 + s ** 2) - Ssym / (1 + s ** 2)) == 0
    out["sign_factor_positive"] = True       # 1/(1+s^2) > 0 for real s
    return out


def check_inverse(Kmax=4):
    u = sp.symbols("u")
    lhs = 1 / (1 + u)
    rhs = sum((-u) ** k for k in range(Kmax + 2)) + (-u) ** (Kmax + 2) / (1 + u)
    g = sp.symbols("g", positive=True)
    a = 9
    return dict(neumann=sp.simplify(lhs - rhs) == 0,
                fold=sp.simplify(g ** a - g ** (Kmax + 1) * g ** (a - Kmax - 1)) == 0)


def main(out):
    rec = dict(schema="review-s46-sympy-forms/1", T=check_T(), regimes=check_regimes(),
               inverse=check_inverse())
    rec["ok"] = (rec["T"]["residual"] == "0" and rec["T"]["control_residual_nonzero"]
                 and all(rec["regimes"].values()) and all(rec["inverse"].values()))
    json.dump(rec, open(out, "w"), indent=1)
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "sympy_forms.json"))
