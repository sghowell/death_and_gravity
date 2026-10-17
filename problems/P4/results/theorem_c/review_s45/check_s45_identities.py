"""S4-5 review: independent sympy checks of the prose identities of the
Theorem C chain (no repository code imported).

  I1  elimination of eta from the (chi, eta) pair: chi'' - T chi' - U chi = 0,
      T = a + d + b'/b, U = a' + b c - a (T - a)            (s4-reduction 1.3)
  I2  Liouville normal form: v = chi/f, 2 f'/f = T  =>  v'' - Qv = (chi'' - T chi' - U chi)/f,
      Q = U + T^2/4 - T'/2                                   (s4-energy 1)
  I3  the polynomial forms used by tc_energy.Atoms / tc_commutator.build_pairs:
      T = TN1/(SD Bc D2), U = U1/(SD Bc D2^2), Q = QNm/(SD^3 Bc^2 D2^2),
      with D1 = SD d/dx and TN1, TD1, TN1d, UN3c, U1, QNm as in tc_energy.py:100-108
  I4  commutation (m = 1): chi1 = chi' solves chi1'' + A1 chi1' + B1 chi1 = 0,
      A1 = A0 - B0'/B0, B1 = B0 + A0' - A0 B0'/B0              (s4-commutator 2)
  I5  the chi1 energy identity (s4-commutator 3, assembly S6) as a pointwise
      differential identity:  BT' = g|u'|^2 + g nu Im(conj(u) u') - g K |u|^2
      for every solution u of u'' + A1 u' + B1 u = 0, with g' = rho g,
      rho = Re A1 - c tt, nu = Im A1 + tau tt, f = (c + tau^2/c) tt - (tau/c) nu,
      K = Re B1 - (tau/c) Im B1 - (f' + rho f)/2,
      BT = Re[(1 + i tau/c) g conj(u) u'] + g f |u|^2/2
  I6  the layer normalization (tc_layer.sigma_box): Sigma = S lam^2/tau^2,
      S = K + nu^2/4, transcribed term by term from the code
  I7  negative controls: perturbing one sign in I5 and in I6 gives a nonzero residual

Run:  python check_s45_identities.py      (sympy only; about 10-60 s)
"""
import sympy as sp

x = sp.symbols('x', real=True)
ok_all = True


def report(name, expr):
    global ok_all
    r = sp.simplify(sp.expand(expr))
    good = (r == 0)
    ok_all &= good
    print(f"{name}: residual {'0' if good else r}")
    return good


# ---------------------------------------------------------------- I1, I2
a, b, c_, d = (sp.Function(n)(x) for n in ('a', 'b', 'c', 'd'))
chi, eta = sp.Function('chi')(x), sp.Function('eta')(x)
sysrules = {sp.Derivative(chi, x): a * chi + b * eta,
            sp.Derivative(eta, x): c_ * chi + d * eta}
T = a + d + sp.diff(b, x) / b
U = sp.diff(a, x) + b * c_ - a * (T - a)
chi1 = a * chi + b * eta                      # chi'
chi2 = sp.diff(chi1, x).subs(sysrules)        # chi''
report("I1 elimination  chi'' - T chi' - U chi", chi2 - T * chi1 - U * chi)

Tf, Uf = sp.Function('T')(x), sp.Function('U')(x)
X = sp.Function('X')(x)                       # generic chi
f = sp.Function('f')(x)
v = X / f
frule = {sp.Derivative(f, (x, 2)): sp.diff(Tf * f / 2, x).subs(sp.Derivative(f, x), Tf * f / 2),
         sp.Derivative(f, x): Tf * f / 2}
lhs = (sp.diff(v, x, 2) - (Uf + Tf**2 / 4 - sp.diff(Tf, x) / 2) * v)
lhs = lhs.subs(sp.Derivative(f, (x, 2)), frule[sp.Derivative(f, (x, 2))]).subs(sp.Derivative(f, x), Tf * f / 2)
report("I2 normal form  v'' - Q v - (chi'' - T chi' - U chi)/f",
       lhs - (sp.diff(X, x, 2) - Tf * sp.diff(X, x) - Uf * X) / f)

# ---------------------------------------------------------------- I3
SD, Ac, Bc, Ae, Be, D2 = (sp.Function(n)(x) for n in ('SD', 'Ac', 'Bc', 'Ae', 'Be', 'D2'))
D1 = lambda F: SD * sp.diff(F, x)             # tc_reduce.Reduction.D1 = S Delta~ d/dx
aa, bb, cc, dd = Ac / D2, Bc / D2, Ae / D2, Be / D2
T3 = aa + dd + sp.diff(bb, x) / bb
U3 = sp.diff(aa, x) + bb * cc - aa * (T3 - aa)
Q3 = U3 + T3**2 / 4 - sp.diff(T3, x) / 2
dAc, dBc, dBe, dD2, dSD = D1(Ac), D1(Bc), D1(Be), D1(D2), D1(SD)
ddBc, ddD2 = D1(dBc), D1(dD2)
TN1 = SD * Bc * (Ac + Be) + dBc * D2 - Bc * dD2
TD1 = dD2 * Bc + D2 * dBc
UN3c = dAc * D2 + SD * (Bc * Ae - Ac * Be)
U1 = Bc * UN3c - Ac * dBc * D2
TN1d = SD * SD * (dBc * (Ac + Be) + Bc * (D1(Ac) + dBe)) + SD * (ddBc * D2 - Bc * ddD2) \
    - dSD * (dBc * D2 - Bc * dD2)
QNm = U1 * SD * SD * Bc + TN1 * TN1 * SD / 4 - (TN1d * D2 * Bc - SD * TN1 * TD1) / 2
report("I3a T  = TN1/(SD Bc D2)", sp.together(T3 - TN1 / (SD * Bc * D2)))
report("I3b U  = U1/(SD Bc D2^2)", sp.together(U3 - U1 / (SD * Bc * D2**2)))
report("I3c Q  = QNm/(SD^3 Bc^2 D2^2)", sp.together(Q3 - QNm / (SD**3 * Bc**2 * D2**2)))

# ---------------------------------------------------------------- I4
A0, B0 = sp.Function('A0')(x), sp.Function('B0')(x)
Y = sp.Function('Y')(x)
A1 = A0 - sp.diff(B0, x) / B0
B1 = B0 + sp.diff(A0, x) - A0 * sp.diff(B0, x) / B0
Lop = lambda F: sp.diff(F, x, 2) + A0 * sp.diff(F, x) + B0 * F
y1 = sp.diff(Y, x)
comm = sp.diff(y1, x, 2) + A1 * sp.diff(y1, x) + B1 * y1
# operator identity: chi1'' + A1 chi1' + B1 chi1 = (Lchi)' - (B0'/B0) Lchi for every Y
report("I4 commutation  (operator identity)",
       sp.together(comm - (sp.diff(Lop(Y), x) - sp.diff(B0, x) / B0 * Lop(Y))))

# ---------------------------------------------------------------- I5
# Values and first derivatives as real symbols; d/dx by the chain rule with
#   g' = rho g, u'' = -A1 u' - B1 u, tt' = tt1, Im A1' = pi1 (f' needs them).
tau, cc0 = sp.symbols('tau c', real=True, nonzero=True)
ur, ui, ur1, ui1, prr, pii, qr, qi, tt, tt1, pi1, g = sp.symbols(
    'ur ui ur1 ui1 pr pi qr qi tt tt1 pi1 g', real=True)
rho = prr - cc0 * tt
nu = pii + tau * tt
fF = (cc0 + tau**2 / cc0) * tt - (tau / cc0) * nu
fF1 = (cc0 + tau**2 / cc0) * tt1 - (tau / cc0) * (pi1 + tau * tt1)      # f'
# u'' = -(pr + i pi)(ur1 + i ui1) - (qr + i qi)(ur + i ui)
upp_r = -(prr * ur1 - pii * ui1) - (qr * ur - qi * ui)
upp_i = -(prr * ui1 + pii * ur1) - (qr * ui + qi * ur)
reb = ur * ur1 + ui * ui1               # Re(conj(u) u') = (|u|^2)'/2
imb = ur * ui1 - ui * ur1               # Im(conj(u) u')


def I5_resid(Kexpr):
    BT = g * (reb - (tau / cc0) * imb) + g * fF * (ur**2 + ui**2) / 2
    deriv = {g: rho * g, ur: ur1, ui: ui1, ur1: upp_r, ui1: upp_i, tt: tt1, pii: pi1}
    dBT = sum(sp.diff(BT, s) * ds for s, ds in deriv.items())
    rhs = g * (ur1**2 + ui1**2) + g * nu * imb - g * Kexpr * (ur**2 + ui**2)
    return sp.expand(dBT - rhs)


K = qr - (tau / cc0) * qi - (fF1 + rho * fF) / 2
report("I5 chi1 energy identity  BT' - [g|u'|^2 + g nu b - g K |u|^2]", I5_resid(K))

# ---------------------------------------------------------------- I6
eps, lam = sp.symbols('epsilon lambda', positive=True)
Vr, Vi, Vpr, Vpi, Lr, Li, Lpr, Lpi, Yr, Yi = sp.symbols('Vr Vi Vpr Vpi Lr Li Lpr Lpi Yr Yi', real=True)
ttv, ttp, u2 = sp.symbols('tt ttp u2', real=True)
kap = eps + sp.I * tau
V_, Vp_, L_, Lp_, Y_ = Vr + sp.I * Vi, Vpr + sp.I * Vpi, Lr + sp.I * Li, Lpr + sp.I * Lpi, Yr + sp.I * Yi
# true blocks:  A0 = V - kappa tt,  A0' = V' - kappa tt',  B0 = u2 kappa^2 + kappa Y,
#               (log B0)' = L, (log B0)'' = L'
A0v, A0pv, B0v = V_ - kap * ttv, Vp_ - kap * ttp, u2 * kap**2 + kap * Y_
A1v, A1pv = A0v - L_, A0pv - Lp_
B1v = B0v + A0pv - A0v * L_
re = lambda z: sp.re(sp.expand(z))
im = lambda z: sp.im(sp.expand(z))
nuv, nupv = im(A1v) + tau * ttv, im(A1pv) + tau * ttp
rhov = re(A1v) - cc0 * ttv
fv = (cc0 + tau**2 / cc0) * ttv - (tau / cc0) * nuv
fpv = (cc0 + tau**2 / cc0) * ttp - (tau / cc0) * nupv
Kv = re(B1v) - (tau / cc0) * im(B1v) - (fpv + rhov * fv) / 2
Sv = Kv + nuv**2 / 4
# tc_layer.sigma_box, transcribed (hatted blocks, iota = 1/tau, E = eps)
it = 1 / tau
Vh, Vph, Lh, Lph, Yh = lam * V_, lam**2 * Vp_, lam * L_, lam**2 * Lp_, lam**2 * Y_
t, tp, u2h = lam * ttv, lam**2 * ttp, lam * u2
W, Wp = Vh - Lh, Vph - Lph
G = Yh + t * Lh - tp
H = Vph - Vh * Lh
nuh, nuph = im(W), im(Wp)
rhoh = re(W) - (eps + cc0) * t
cf = cc0 * it * it + 1 / cc0


def sigma_code(sgn=1):
    reB1 = u2h * lam * (eps * eps * it * it - 1) + eps * it * it * re(G) - it * im(G) + it * it * re(H)
    imB1 = -(2 * eps * u2h * lam + re(G) + eps * it * im(G) + it * im(H)) / cc0
    fp = -cf * tp / 2 + sgn * (it / cc0) * nuph / 2
    rhof = -rhoh * (cf * t - (it / cc0) * nuh) / 2
    nu2 = nuh * nuh * it * it / 4
    return reB1 + imB1 + fp + rhof + nu2


report("I6 layer normalization  Sigma - S lam^2/tau^2", sp.together(sigma_code() - Sv * lam**2 / tau**2))

# ---------------------------------------------------------------- I7 controls
Kbad = qr + (tau / cc0) * qi - (fF1 + rho * fF) / 2
r7a = I5_resid(Kbad)
r7b = sp.simplify(sp.together(sigma_code(-1) - Sv * lam**2 / tau**2))
print("I7 controls (must be nonzero):", r7a != 0, r7b != 0)
ok_all &= (r7a != 0) and (r7b != 0)

print("ALL IDENTITIES HOLD" if ok_all else "SOME CHECK FAILED")
