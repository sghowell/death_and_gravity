"""S3 review: negative controls for sympy_identities.py (a wrong partial-fraction coefficient and the
wrong sign of 2 i w r^2 in the factored Wronskian must both return False)."""
import time, sympy as sp
t0=time.time()
r = sp.symbols("r", real=True); w, L = sp.symbols("w L")
r1, r2, r3, r4, lam = sp.symbols("r1 r2 r3 r4 Lam", positive=True)
rts=[r1,r2,r3,r4]
Dq = -(lam/3)*sp.Mul(*[(r-x) for x in rts]); dDq = sp.diff(Dq, r)
A = [x**2/dDq.subs(r, x) for x in rts]
print("S1 true:", sp.simplify(sum(A[i]/(r-rts[i]) for i in range(4)) - r**2/Dq) == 0, time.time()-t0)
Abad = [x**3/dDq.subs(r, x) for x in rts]
print("S1 wrong (r_i^3):", sp.simplify(sum(Abad[i]/(r-rts[i]) for i in range(4)) - r**2/Dq) == 0, time.time()-t0)
D = sp.Function("D")(r); Ep = sp.Function("Ep")(r); Rp_ = sp.Function("Rp")(r); Rc_ = sp.Function("Rc")(r)
cp, cc = sp.symbols("c_p c_c", real=True); rp_, rc_ = sp.symbols("r_p r_c", positive=True)
dEp = sp.I*w*r**2/D*Ep
Zp = Ep*sp.exp(-sp.I*w*cp)*Rp_; Zc = sp.exp(sp.I*w*cc)*Rc_/Ep
dZp = sp.diff(Zp, r).subs(sp.Derivative(Ep, r), dEp); dZc = sp.diff(Zc, r).subs(sp.Derivative(Ep, r), dEp)
Wf = D*(Rp_*sp.diff(Rc_, r) - sp.diff(Rp_, r)*Rc_)*sp.exp(-sp.I*w*cp + sp.I*w*cc)/(rp_*rc_)
for sgn in (2, -2):
    Wz = (D*(Zp*dZc - dZp*Zc) + sgn*sp.I*w*r**2*Zp*Zc)/(rp_*rc_)
    print("S3 with", sgn, "i w r^2:", sp.simplify(sp.expand(Wz - Wf)) == 0)
print(time.time()-t0)
