"""From-scratch symbolic derivation of the KHA CSS system from Einstein-Euler.

No sympy.diffgeom, no tensor package: plain index loops over 4x4 sympy matrices.

Chain (every step is a function returning sympy objects; ``derive()`` collects them):
  1. metric  ds^2 = -alpha(t,r)^2 dt^2 + a(t,r)^2 dr^2 + r^2 dOmega^2
  2. Gamma^l_{mn} = 1/2 g^{ls}(d_m g_{sn} + d_n g_{sm} - d_s g_{mn})
  3. R_{mn} = d_l Gamma^l_{mn} - d_n Gamma^l_{ml} + Gamma^l_{ls} Gamma^s_{mn}
              - Gamma^l_{ns} Gamma^s_{ml};  G_{mn} = R_{mn} - g_{mn} R/2   (MTW signs)
  4. T^{mn} = (rho + p) u^m u^n + p g^{mn}, p = rho/3, u^t = 1/(alpha sqrt(1-V^2)),
     u^r = V/(a sqrt(1-V^2))  (V = radial velocity seen by the static observer)
  5. E_{mn} := G_{mn} - 8 pi T_{mn};  U^n := nabla_m T^{mn}
  6. first-jet form (fields and first derivatives as symbols), then the CSS ansatz
     t = -e^{-s}, r = e^{x-s} with the KHA unknowns
         N = alpha/(a e^x),  A = a^2,  W = 4 pi r^2 a^2 rho,  V,
     i.e. alpha = N a e^x, rho = W e^{2s}/(4 pi e^{2x} a^2);  d_t = e^s(d_s + d_x),
     d_r = e^{s-x} d_x.
  7. exact identities with the KHA rows of ``p4.css.coeffs`` (``identities``).
Results are exact: every residual is decided by ``sympy.cancel`` on rational functions
in the symbols (the only transcendental objects, e^s and e^x, are replaced by the
positive symbols es, ex; no square roots occur because a, not A, is the primitive field).
"""
from __future__ import annotations

import sympy as sp

# --- coordinates, fields ------------------------------------------------------------
t, r, th, ph = sp.symbols("t r theta phi", real=True)
COORDS = (t, r, th, ph)
alpha = sp.Function("alpha")(t, r)
a = sp.Function("a")(t, r)
rho = sp.Function("rho")(t, r)
V = sp.Function("V")(t, r)
FIELDS = (alpha, a, rho, V)

# first-jet symbols (values and first partials of the four fields)
J = sp.symbols("al al_t al_r aa aa_t aa_r rh rh_t rh_r Vv Vv_t Vv_r", real=True)
al, al_t, al_r, aa, aa_t, aa_r, rh, rh_t, rh_r, Vv, Vv_t, Vv_r = J

# similarity variables and KHA unknowns (functions of (s, x)), plus their jet symbols
s, x = sp.symbols("s x", real=True)
es, ex = sp.symbols("es ex", positive=True)               # e^s, e^x
Nf, af, Wf, Vf = (sp.Function(n)(s, x) for n in ("N", "a", "W", "V"))
K = sp.symbols("N N_s N_x a_ a_s a_x W W_s W_x V V_s V_x", real=True)
N_, N_s, N_x, a_, a_s, a_x, W_, W_s, W_x, V_, V_s, V_x = K


# --- 1-3: geometry -------------------------------------------------------------------
def metric():
    g = sp.diag(-alpha**2, a**2, r**2, r**2 * sp.sin(th) ** 2)
    return g, g.inv()


def christoffel(g, ginv):
    n = 4
    Gam = [[[0] * n for _ in range(n)] for _ in range(n)]
    for l in range(n):
        for m in range(n):
            for nu in range(n):
                tot = 0
                for sg in range(n):
                    tot += ginv[l, sg] * (sp.diff(g[sg, nu], COORDS[m]) + sp.diff(g[sg, m], COORDS[nu])
                                          - sp.diff(g[m, nu], COORDS[sg]))
                Gam[l][m][nu] = sp.cancel(tot / 2)
    return Gam


def ricci(Gam):
    R = sp.zeros(4)
    for m in range(4):
        for nu in range(4):
            tot = 0
            for l in range(4):
                tot += sp.diff(Gam[l][m][nu], COORDS[l]) - sp.diff(Gam[l][m][l], COORDS[nu])
                for sg in range(4):
                    tot += Gam[l][l][sg] * Gam[sg][m][nu] - Gam[l][nu][sg] * Gam[sg][m][l]
            R[m, nu] = sp.cancel(tot)
    return R


def einstein(g, ginv, R):
    Rs = sp.cancel(sum(ginv[i, j] * R[i, j] for i in range(4) for j in range(4)))
    return (R - g * Rs / 2).applyfunc(sp.cancel), Rs


# --- 4-5: fluid ----------------------------------------------------------------------
def stress_tensor(g, ginv):
    """T^{mn}, T_{mn} and the check u.u (must be -1).  Only u^m u^n products occur."""
    p = rho / 3
    W2 = 1 / (1 - V**2)
    uu = sp.zeros(4)
    uu[0, 0] = W2 / alpha**2
    uu[0, 1] = uu[1, 0] = V * W2 / (alpha * a)
    uu[1, 1] = V**2 * W2 / a**2
    T_up = ((rho + p) * uu + p * ginv).applyfunc(sp.cancel)
    T_dn = (g * T_up * g).applyfunc(sp.cancel)
    norm = sp.cancel(sum(g[i, j] * uu[i, j] for i in range(4) for j in range(4)))
    return T_up, T_dn, norm


def divergence(T_up, Gam):
    """U^n = nabla_m T^{mn} = d_m T^{mn} + Gamma^m_{ml} T^{ln} + Gamma^n_{ml} T^{ml}."""
    out = []
    for nu in range(4):
        tot = 0
        for m in range(4):
            tot += sp.diff(T_up[m, nu], COORDS[m])
            for l in range(4):
                tot += Gam[m][m][l] * T_up[l, nu] + Gam[nu][m][l] * T_up[m, l]
        out.append(sp.cancel(tot))
    return out


# --- 6: jets and the CSS ansatz -----------------------------------------------------------
def to_jet(expr):
    """Replace fields and their first derivatives by the symbols J; assert first order."""
    for D in expr.atoms(sp.Derivative):
        assert sum(c for _, c in D.variable_count) == 1, f"second derivative present: {D}"
    rep = {}
    for f, (v, vt, vr) in zip(FIELDS, [J[i:i + 3] for i in (0, 3, 6, 9)]):
        rep[sp.Derivative(f, t)] = vt
        rep[sp.Derivative(f, r)] = vr
    expr = expr.subs(rep)
    return expr.subs({f: v for f, v in zip(FIELDS, (al, aa, rh, Vv))})


def css_fields():
    """alpha, a, rho, V as functions of (s, x) through the KHA unknowns."""
    alpha_S = Nf * af * sp.exp(x)
    rho_S = Wf * sp.exp(2 * s) / (4 * sp.pi * sp.exp(2 * x) * af**2)   # r^2 = e^{2x-2s}
    return alpha_S, af, rho_S, Vf


def chain(f):
    """(d_t f, d_r f) for f = F(s, x):  d_t = e^s (d_s + d_x),  d_r = e^{s-x} d_x."""
    return sp.exp(s) * (sp.diff(f, s) + sp.diff(f, x)), sp.exp(s - x) * sp.diff(f, x)


def to_similarity(jet_expr):
    """Jet expression in (t, r) -> rational function of the (s, x)-jet symbols K, es, ex."""
    rep = {t: -sp.exp(-s), r: sp.exp(x - s)}
    for (v, vt, vr), F in zip([J[i:i + 3] for i in (0, 3, 6, 9)], css_fields()):
        Ft, Fr = chain(F)
        rep.update({v: F, vt: Ft, vr: Fr})
    e = jet_expr.subs(rep, simultaneous=True).doit()
    drep = {}
    for F, (v, vs, vx) in zip((Nf, af, Wf, Vf), [K[i:i + 3] for i in (0, 3, 6, 9)]):
        drep[sp.Derivative(F, s)] = vs
        drep[sp.Derivative(F, x)] = vx
    e = e.subs(drep).subs({Nf: N_, af: a_, Wf: W_, Vf: V_})
    return sp.cancel(e.subs({s: sp.log(es), x: sp.log(ex)}))


# --- driver -----------------------------------------------------------------------------
def derive():
    """Every intermediate object of the derivation, as a dict of sympy expressions."""
    g, ginv = metric()
    Gam = christoffel(g, ginv)
    R = ricci(Gam)
    G, Rs = einstein(g, ginv, R)
    T_up, T_dn, norm = stress_tensor(g, ginv)
    U = divergence(T_up, Gam)
    E = (G - 8 * sp.pi * T_dn).applyfunc(sp.cancel)
    eqs_tr = {"E_tt": to_jet(E[0, 0]), "E_tr": to_jet(E[0, 1]), "E_rr": to_jet(E[1, 1]),
              "U_t": to_jet(U[0]), "U_r": to_jet(U[1])}
    eqs_sx = {k: to_similarity(v) for k, v in eqs_tr.items()}
    return {"g": g, "ginv": ginv, "Gamma": Gam, "Ricci": R, "Rscalar": Rs, "Einstein": G,
            "T_up": T_up, "T_dn": T_dn, "u_norm": norm, "Euler": U, "E": E, "jets_tr": eqs_tr,
            "jets_sx": eqs_sx, "U_angular": (U[2], U[3]),
            "E_offdiag": (E[0, 2], E[0, 3], E[1, 2], E[1, 3], E[2, 3])}


# --- 7: exact identities with the KHA rows of p4.css.coeffs ------------------------------------
A_ = sp.Symbol("A", real=True)


def kha_rows_symbolic():
    """Rows of ``p4.css.coeffs`` on symbols (A -> a_^2).  The two float literals of css.py
    (8.0/3.0 in G, 4.0/3.0 in e) are the only non-exact objects; they are restored to
    8/3, 4/3 and this is asserted (no other Float may occur)."""
    from p4 import css
    q = {}
    for k, v in css.coeffs(A_, N_, W_, V_).items():
        v = sp.sympify(v)
        for fl in v.atoms(sp.Float):
            assert abs(sp.nsimplify(fl, rational=True)) in (sp.Rational(8, 3), sp.Rational(4, 3)), fl
        q[k] = sp.cancel(sp.nsimplify(v, rational=True).subs(A_, a_**2))
    rows = {
        "row1": 2 * a_x / a_ - q["FA"],                        # A_x/A - F_A   (Hamiltonian)
        "row2": N_x / N_ - q["FN"],                            # N_x/N - F_N   (slicing)
        "rowM": 2 * (a_s + a_x) / a_ - q["G"],                 # (A_s+A_x)/A - G (momentum)
        "row3": q["sa"] * W_s + q["sb"] * V_s + q["a"] * W_x + q["b"] * V_x + q["e"],
        "row4": q["sc"] * W_s + q["sd"] * V_s + q["c"] * W_x + q["d"] * V_x + q["f"],
    }
    return q, rows


def identities(d=None):
    """Exact residuals (all must be 0) relating the derived equations to the KHA rows.

    I1: E_tt = N^2 e^{2s} row1;  I2: E_tr = e^{2s-x} rowM;  I3: E_rr = e^{2s-2x} (2 row2 + row1);
    I4: with rowM = row1 = row2 = 0 used to eliminate a_s, a_x, N_x, the Euler pair equals
        M (row3, row4) with M the 2x2 matrix fixed by the (W_s, V_s) coefficients; the
        residual is the *identical* vanishing of the remaining (W_x, V_x, source) parts.
    """
    d = d or derive()
    q, rows = kha_rows_symbolic()
    E = d["jets_sx"]
    res = {"I1": sp.cancel(E["E_tt"] - N_**2 * es**2 * rows["row1"]),
           "I2": sp.cancel(E["E_tr"] - es**2 / ex * rows["rowM"]),
           "I3": sp.cancel(E["E_rr"] - es**2 / ex**2 * (2 * rows["row2"] + rows["row1"]))}
    elim = {a_s: a_ * q["G"] / 2 - a_x}                       # momentum constraint
    U = [sp.cancel(E[k].subs(elim).subs(a_x, a_ * q["FA"] / 2).subs(N_x, N_ * q["FN"]))
         for k in ("U_t", "U_r")]
    R = [rows["row3"], rows["row4"]]
    LU = sp.Matrix([[sp.diff(u, v) for v in (W_s, V_s)] for u in U])
    LR = sp.Matrix([[sp.diff(rw, v) for v in (W_s, V_s)] for rw in R])
    M = (LU * LR.inv()).applyfunc(sp.cancel)
    UR = M * sp.Matrix(R)
    res["I4_t"] = sp.cancel(U[0] - UR[0])
    res["I4_r"] = sp.cancel(U[1] - UR[1])
    return {"residuals": res, "M": M, "detM": sp.factor(M.det()), "rows": rows, "coeffs": q}


def sonic_system_rows():
    """Rows P.u' - Q of the certified polynomial system ``validated.systems.sonic_system``
    as sympy expressions in (a_, N_, W_, V_, a_x, N_x, W_x, V_x); requires python-flint."""
    from p4.validated.systems import sonic_system
    sys_ = sonic_system()
    gens = (sp.Symbol("t"), A_, N_, W_, V_)

    def conv(p):
        tot = sp.Integer(0)
        for exps, c in p.terms():
            term = sp.Rational(int(c.p), int(c.q))
            for g, e in zip(gens, exps):
                term *= g ** int(e)
            tot += term
        return tot.subs(A_, a_**2)

    up = (2 * a_ * a_x, N_x, W_x, V_x)                          # u' = (A', N', W', V')
    return [sp.expand(sum(conv(sys_.P[i][j]) * up[j] for j in range(4)) - conv(sys_.Q[i]))
            for i in range(4)]


def identities_with_certified_system(rows=None):
    """P.u' - Q of sonic_system == (A S row1, N row2, 3 W S row3, W S row4) at W_s = V_s = 0."""
    rows = rows or kha_rows_symbolic()[1]
    css_ = {k: v.subs({W_s: 0, V_s: 0}) for k, v in rows.items()}
    S = 1 - V_**2
    target = [a_**2 * S * css_["row1"], N_ * css_["row2"], 3 * W_ * S * css_["row3"],
              W_ * S * css_["row4"]]
    return [sp.cancel(p - q) for p, q in zip(sonic_system_rows(), target)]
