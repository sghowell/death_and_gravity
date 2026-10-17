"""Symbolic checks (sympy) of the identities used in notes/s2-m3.md and in m3.py (reviewer script).

Each item prints True when the identity reduces to 0.  u = u1 + i u2, v = v1 + i v2 with real
functions; phi, G, c, q real functions of s; Q, V complex.
Run (repo root):  .venv/bin/python problems/P2/certificates/review_s2/indep_identities_sympy.py
"""
import sympy as sp

s = sp.symbols("s", real=True)
eta, a, b, th, sgn = sp.symbols("eta a b theta sgn", real=True)
u1, u2, v1, v2, phi, G, c, q, Q1, Q2 = [sp.Function(n, real=True)(s) for n in
                                        ("u1", "u2", "v1", "v2", "phi", "G", "c", "q", "Q1", "Q2")]
I = sp.I
u = u1 + I * u2
ub = u1 - I * u2
res = {}


def zero(e):
    return sp.simplify(sp.expand(e)) == 0


# 1. integration-by-parts identity behind section 3
lhs = sp.exp(-I * phi) * (-sp.diff(u, s, 2)) * ub
rhs = (-sp.diff(sp.exp(-I * phi) * sp.diff(u, s) * ub, s) + sp.exp(-I * phi) * sp.diff(u, s) * sp.diff(ub, s)
       - I * sp.diff(phi, s) * sp.exp(-I * phi) * sp.diff(u, s) * ub)
res["1 e^{-i phi}(-u'')ubar = -(e^{-i phi}u'ubar)' + e^{-i phi}|u'|^2 - i phi' e^{-i phi} u' ubar"] = zero(lhs - rhs)

# 2. Young with eps1 = eta|phi'|:  eta P A^2 + P B^2/(4 eta) - P A B = P (sqrt(eta) A - B/(2 sqrt(eta)))^2
P, Aa, Bb = sp.symbols("P A B", nonnegative=True)
etap = sp.symbols("etap", positive=True)
res["2 Young: eta P A^2 + P B^2/(4 eta) - P A B is a square times P"] = zero(
    etap * P * Aa ** 2 + P * Bb ** 2 / (4 * etap) - P * Aa * Bb - P * (sp.sqrt(etap) * Aa - Bb / (2 * sp.sqrt(etap))) ** 2)
# Re(-i phi' e^{-i phi} u' ubar) >= -|phi'||u'||u|  (|z| bound; nothing to simplify) -- recorded
res["2b Re(-i phi' e^{-i phi} u'ubar) >= -|phi'||u'||u| (modulus bound)"] = True

# 3. Lemma 4 ground-state identity, u = g v, g = e^{-G}:
#    c|u'|^2 + q|u|^2 = c g^2 |v'|^2 + E g^2 |v|^2 + (c g g' |v|^2)',  E = q + c'G' - c G'^2 + c G''
g = sp.exp(-G)
v = v1 + I * v2
vb = v1 - I * v2
uu = g * v
uub = g * vb
abs_up2 = sp.diff(uu, s) * sp.diff(uub, s)
E = q + sp.diff(c, s) * sp.diff(G, s) - c * sp.diff(G, s) ** 2 + c * sp.diff(G, s, 2)
lhs = c * abs_up2 + q * uu * uub
rhs = c * g ** 2 * sp.diff(v, s) * sp.diff(vb, s) + E * g ** 2 * v * vb + sp.diff(c * g * sp.diff(g, s) * v * vb, s)
res["3 Lemma 4: c|u'|^2 + q|u|^2 = c g^2|v'|^2 + E|u|^2 + (c g g'|v|^2)'"] = zero(lhs - rhs)

# 4. E' as coded in m3.py (G''' = 0 on a cell)
Gp, Gpp = sp.symbols("Gp Gpp", real=True)          # G' = Gp + Gpp (s - s0) on the cell
s0 = sp.symbols("s0", real=True)
Gp_s = Gp + Gpp * (s - s0)
Ecell = q + sp.diff(c, s) * Gp_s - c * Gp_s ** 2 + c * Gpp
Ep_code = (sp.diff(q, s) + sp.diff(c, s, 2) * Gp_s + 2 * sp.diff(c, s) * Gpp - sp.diff(c, s) * Gp_s ** 2
           - 2 * c * Gp_s * Gpp)
res["4 E' formula of Checker.E_parts"] = zero(sp.diff(Ecell, s) - Ep_code)

# 5. c, c', c'' with |phi'| = sgn phi'  (sgn constant on a cell)
cc = sp.cos(phi) - eta * sgn * sp.diff(phi, s)
cp_code = -sp.sin(phi) * sp.diff(phi, s) - eta * sgn * sp.diff(phi, s, 2)
cpp_code = -sp.cos(phi) * sp.diff(phi, s) ** 2 - sp.sin(phi) * sp.diff(phi, s, 2) - eta * sgn * sp.diff(phi, s, 3)
res["5a c' as coded"] = zero(sp.diff(cc, s) - cp_code)
res["5b c'' as coded"] = zero(sp.diff(cc, s, 2) - cpp_code)

# 6. q' as coded: q = Re(e^{-i phi} Q) - sgn phi'/(4 eta);  q' = Re(e^{-i phi}(Q' - i phi' Q)) - sgn phi''/(4 eta)
# (explicit real parts: sp.re does not know that Derivative(Q1(s), s) is real)
def re_rot(phase, X1, X2):
    """Re(e^{-i phase}(X1 + i X2)) for real X1, X2."""
    return sp.cos(phase) * X1 + sp.sin(phase) * X2


qq = re_rot(phi, Q1, Q2) - sgn * sp.diff(phi, s) / (4 * eta)
# Q' - i phi' Q = (Q1' + phi' Q2) + i (Q2' - phi' Q1)
qp_code = (re_rot(phi, sp.diff(Q1, s) + sp.diff(phi, s) * Q2, sp.diff(Q2, s) - sp.diff(phi, s) * Q1)
           - sgn * sp.diff(phi, s, 2) / (4 * eta))
res["6 q' as coded"] = zero(sp.diff(qq, s) - qp_code)
res["6b Re(e^{-i phi}(Q1 + i Q2)) = cos(phi) Q1 + sin(phi) Q2"] = zero(
    sp.expand_complex(sp.re(sp.expand_complex(sp.exp(-I * phi) * (Q1 + I * Q2)))) - re_rot(phi, Q1, Q2))

# 7. dQ/ds = e^{3 i theta} V'(r) F(r) when dr/ds = e^{i theta} F(r), Q = e^{2 i theta}(V(r) - w^2)
r = sp.Function("r")(s)
Vf, Ff = sp.Function("V"), sp.Function("F")
w2 = sp.symbols("w2")
Qexpr = sp.exp(2 * I * th) * (Vf(r) - w2)
dQ = sp.diff(Qexpr, s).subs(sp.Derivative(r, s), sp.exp(I * th) * Ff(r))
res["7 Q' = e^{3 i theta} V'(r) F(r)"] = zero(dQ - sp.exp(3 * I * th) * sp.diff(Vf(r), r) * Ff(r))

# 8. V = D K / r^6 (K = L r^2 + r D' - 2D) equals F (L/r^2 + F'/r), F = D/r^2; and V' as coded
rr, L, Lam, Q2s = sp.symbols("r L Lambda Qsq", positive=True)
D = -Lam / 3 * rr ** 4 + rr ** 2 - 2 * rr + Q2s
F = D / rr ** 2
V_def = F * (L / rr ** 2 + sp.diff(F, rr) / rr)
Dp = sp.diff(D, rr)
K = L * rr ** 2 + rr * Dp - 2 * D
V_code = D * K / rr ** 6
Kp_code = 2 * L * rr + Dp + rr * sp.diff(D, rr, 2) - 2 * Dp
N = D * K
Np_code = Dp * K + D * Kp_code
Vp_code = (rr * Np_code - 6 * N) / rr ** 7
res["8a V_l = D K / r^6"] = zero(V_def - V_code)
res["8b K' as coded"] = zero(sp.diff(K, rr) - Kp_code)
res["8c V' as coded"] = zero(sp.diff(V_def, rr) - Vp_code)
res["8d F as coded (u = 1/r): -(Lam/3) r^2 + 1 - 2/r + Q^2/r^2"] = zero(F - (-Lam / 3 * rr ** 2 + 1 - 2 / rr + Q2s / rr ** 2))

# 9. rotated equation: u(s) = Phi(x_P + e^{i theta} s), Phi'' = (V - w^2) Phi  =>  u'' = e^{2 i theta}(V - w^2) u
Phi = sp.Function("Phi")
xP = sp.symbols("x_P")
uS = Phi(xP + sp.exp(I * th) * s)
res["9 d^2/ds^2 Phi(x_P + e^{i theta}s) = e^{2 i theta} Phi''"] = zero(
    sp.diff(uS, s, 2) - sp.exp(2 * I * th) * sp.Subs(sp.diff(Phi(sp.Symbol("y")), sp.Symbol("y"), 2),
                                                     sp.Symbol("y"), xP + sp.exp(I * th) * s).doit())

# 10. decay rate on the contour: |e^{-i w e^{i pi/4} s}| = e^{s (a + b)/sqrt 2}
w = a + I * b
expo = sp.re(sp.expand_complex(-I * w * sp.exp(I * sp.pi / 4)))
res["10 Re(-i w e^{i pi/4}) = (a + b)/sqrt(2)"] = zero(expo - (a + b) / sp.sqrt(2))

# 11. window hull corners: a^2 - b^2 increasing in a (a > 0) and in b (b <= 0); 2ab decreasing in a, increasing in b
res["11 hull monotonicity (signs of partials)"] = all([
    sp.diff(a ** 2 - b ** 2, a) == 2 * a, sp.diff(a ** 2 - b ** 2, b) == -2 * b,
    sp.diff(2 * a * b, a) == 2 * b, sp.diff(2 * a * b, b) == 2 * a])

# 12. tail lemma: d/ds [s^3/(s1^4 + s^4)^2] = s^2 (3 s1^4 - 5 s^4)/(s1^4 + s^4)^3  (|phi'| decreasing for s^4 > 3 s1^4/5)
s1 = sp.symbols("s1", positive=True)
res["12 |phi'| monotone for |s| > (3/5)^{1/4} s1"] = zero(
    sp.diff(s ** 3 / (s1 ** 4 + s ** 4) ** 2, s) - s ** 2 * (3 * s1 ** 4 - 5 * s ** 4) / (s1 ** 4 + s ** 4) ** 3)

# 13. HO model behind mu_HO: V ~ V0 - k x^2/2 on x - x_P = e^{i pi/4} s gives Q = i(V - w^2) = k s^2/2 + i(V0 - w^2)
k, V0 = sp.symbols("k V0", positive=True)
Qm = I * (V0 - k * (sp.exp(I * sp.pi / 4) * s) ** 2 / 2 - w ** 2)
res["13 Re Q = k s^2/2 + Im(w^2) = k s^2/2 + 2ab (barrier-top model)"] = zero(
    sp.re(sp.expand_complex(Qm)) - (k * s ** 2 / 2 + 2 * a * b))

ok = all(res.values())
for k_, v_ in res.items():
    print(("PASS " if v_ else "FAIL ") + k_)
print("ALL PASS" if ok else "SOME FAIL")
