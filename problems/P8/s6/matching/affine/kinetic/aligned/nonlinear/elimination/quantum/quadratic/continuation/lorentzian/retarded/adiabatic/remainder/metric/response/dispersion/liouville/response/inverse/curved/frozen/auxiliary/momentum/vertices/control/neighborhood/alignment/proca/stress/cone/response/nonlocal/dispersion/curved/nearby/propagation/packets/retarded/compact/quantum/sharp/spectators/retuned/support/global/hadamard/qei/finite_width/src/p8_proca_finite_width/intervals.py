"""Rigorous ball jets on a WHOLE real interval, not sampled times."""
from fractions import Fraction
from functools import cache
from math import factorial

import sympy as sp
from flint import acb, arb, ctx, fmpq
from p8_proca_global_hadamard import canonical
from p8_proca_global_support import model


def scalar(value):
    if isinstance(value,(arb,acb)):
        return acb(value)
    real,imag=sp.sympify(value).as_real_imag()
    real,imag=sp.Rational(real),sp.Rational(imag)
    return acb(fmpq(int(real.p),int(real.q)),fmpq(int(imag.p),int(imag.q)))

class BallJets:
    def __init__(self,degree,radius="1/100"):
        if isinstance(degree,bool) or not isinstance(degree,int) or degree<1:
            raise ValueError("Require a positive exact integer time-jet degree")
        if isinstance(radius,bool) or not isinstance(radius,(str,int,Fraction)):
            raise TypeError("Require an exact rational radius, not binary floating point")
        radius=Fraction(radius)
        if radius<=0:
            raise ValueError("Require a positive exact radius")
        radius=str(radius)
        self.degree=degree
        self.zero=tuple(acb(0) for _ in range(degree+1))
        self.one=(acb(1),)+self.zero[1:]
        self.variable=(acb(arb(0,radius)),acb(1))+self.zero[2:]

    def constant(self,value):
        return (scalar(value),)+self.zero[1:]

    def add(self,*values):
        return tuple(sum((value[j] for value in values),acb(0)) for j in range(self.degree+1))

    def scale(self,value,multiplier):
        return tuple(scalar(multiplier)*v for v in value)

    def mul(self,left,right):
        return tuple(sum((left[j]*right[n-j] for j in range(n+1)),acb(0))
            for n in range(self.degree+1))

    def inv(self,value):
        if value[0].contains(0):
            raise ValueError("A whole-interval pivot enclosure includes zero")
        out=[1/value[0]]
        for n in range(1,self.degree+1):
            out.append(-sum((value[j]*out[n-j] for j in range(1,n+1)),acb(0))/value[0])
        return tuple(out)

    def div(self,left,right):
        return self.mul(left,self.inv(right))

    def derivative(self,value):
        return tuple((j+1)*value[j+1] for j in range(self.degree))+(acb(0),)

    def truncate(self,value,degree):
        return tuple(value[j] if j<=degree else acb(0) for j in range(self.degree+1))

    def polynomial(self,expr):
        out=self.zero
        for c in sp.Poly(expr,model.u).all_coeffs():
            out=self.add(self.mul(out,self.variable),self.constant(c))
        return out

    def rational(self,value):
        num,den=sp.fraction(sp.cancel(value))
        return self.div(self.polynomial(num),self.polynomial(den))

    def sqrt(self,value):
        if not value[0].imag.is_zero() or not value[0].real>0:
            raise ValueError("Require a strictly positive real root enclosure")
        root=acb(value[0].real.sqrt())
        out=[root]
        for n in range(1,self.degree+1):
            out.append((value[n]-sum((out[j]*out[n-j] for j in range(1,n)),acb(0)))/(2*root))
        return tuple(out)

    def matrix(self,rows):
        return tuple(tuple(rows[i][j] for j in range(2)) for i in range(2))

    def mzero(self):
        return self.matrix([[self.zero,self.zero],[self.zero,self.zero]])

    def madd(self,*values):
        return self.matrix([[self.add(*(v[i][j] for v in values)) for j in range(2)] for i in range(2)])

    def mscale(self,value,multiplier):
        return self.matrix([[self.scale(value[i][j],multiplier) for j in range(2)] for i in range(2)])

    def mmul(self,left,right):
        return self.matrix([[self.add(*(self.mul(left[i][r],right[r][j]) for r in range(2)))
            for j in range(2)] for i in range(2)])

    def transpose(self,value):
        return self.matrix([[value[j][i] for j in range(2)] for i in range(2)])

    def mderivative(self,value):
        return self.matrix([[self.derivative(value[i][j]) for j in range(2)] for i in range(2)])

    def mtruncate(self,value,degree):
        return self.matrix([[self.truncate(value[i][j],degree) for j in range(2)] for i in range(2)])

def ceiling_positive(value):
    z=value.upper().ceil().unique_fmpz()
    if z is None:
        raise ValueError("Non-finite or non-unique integer upper enclosure")
    return int(z)

def abs_upper(value):
    return value.real.abs_upper()+value.imag.abs_upper()

def matrix_norm(value,derivative=0):
    a=[[factorial(derivative)*abs_upper(value[i][j][derivative]) for j in range(2)] for i in range(2)]
    return max(ceiling_positive(sum(row,arb(0))) for row in a+list(zip(*a)))

@cache
@ctx.workprec(160)
def build(order=6,derivatives=3,radius="1/100"):
    T=BallJets(order+derivatives,radius)
    bg=model.background()
    ell=T.rational(bg["ell"])
    a=T.rational(bg["a"])
    mass=T.rational(2*bg["J_new"]/bg["lam"]**2)
    c=T.sqrt(T.rational(bg["clock_physical_squared_speed"]))
    omega=(T.div(c,a),T.inv(a))
    Vi=T.matrix([[T.one,T.zero],[ell,T.one]])
    V=T.matrix([[T.one,T.zero],[T.scale(ell,-1),T.one]])
    B=T.matrix([[T.add(T.mul(mass,omega[0]),T.div(T.mul(ell,ell),a)),T.div(ell,a)],
        [T.div(ell,a),T.inv(a)]])
    rs=[T.mscale(B,-sp.I)]
    d=canonical.data()
    hs=[]
    for j in (1,0,-1):
        H=(-d["constant_symplectic_form"]*d["complete_Laurent_coefficients"][j]).subs(
            model.substitution(),simultaneous=True)
        hs.append(tuple(tuple(T.rational(H[i,j]) for j in range(4)) for i in range(4)))
    block=lambda H,row,col:T.matrix([[H[i+row][j+col] for j in range(2)] for i in range(2)])
    Hpp=[block(H,2,2) for H in hs]
    Hpq=[block(H,2,0) for H in hs]
    Hqp=[block(H,0,2) for H in hs]
    Hqq=[block(H,0,0) for H in hs]

    def coefficient(n,current):
        forcing=T.mderivative(current[n-1]) if 1<=n<=len(current) else T.mzero()
        if n<len(Hqq):
            forcing=T.madd(forcing,Hqq[n])
        for r in range(len(hs)):
            for p in range(len(current)):
                if p+r==n:
                    forcing=T.madd(forcing,T.mmul(current[p],Hpq[r]),T.mmul(Hqp[r],current[p]))
                for s in range(len(current)):
                    if p+r+s==n:
                        forcing=T.madd(forcing,T.mmul(T.mmul(current[p],Hpp[r]),current[s]))
        return forcing

    for n in range(1,order+1):
        forcing=coefficient(n,rs)
        transformed=T.mmul(T.mmul(T.transpose(V),forcing),V)
        solved=T.matrix([[T.scale(T.div(transformed[i][j],T.add(omega[i],omega[j])),-sp.I)
            for j in range(2)] for i in range(2)])
        rs.append(T.mtruncate(T.mmul(T.mmul(T.transpose(Vi),solved),Vi),T.degree-n))

    # Configuration generator in the moving V frame is -ik*omega + B.
    # Store B as a finite polynomial in rho=1/k.
    bc={}
    for n in range(order+2):
        X=T.mzero()
        if n==0:
            X=T.madd(X,Hpq[1])
        for p in range(len(rs)):
            for r in (0,2):
                if p+r-1==n:
                    X=T.madd(X,T.mmul(Hpp[r],rs[p]))
        bc[n]=T.mmul(T.mmul(Vi,X),V)
        if n==0:
            bc[n]=T.madd(bc[n],T.mscale(T.mmul(Vi,T.mderivative(V)),-1))
    defects={n:coefficient(n,rs) for n in range(order+1,2*order+3)}
    return {"T":T,"R":rs,"V":V,"Vi":Vi,"omega":omega,"Hpp":Hpp,"Hpq":Hpq,
        "B_coefficients":bc,"defect_coefficients":defects,
        "R_derivative_upper_integers":[[matrix_norm(R,j) for j in range(derivatives+1)] for R in rs],
        "B_coefficient_derivative_upper_integers":{n:[matrix_norm(v,j) for j in range(derivatives+1)]
            for n,v in bc.items()},
        "defect_coefficient_upper_integers":{n:matrix_norm(v) for n,v in defects.items()},
        "omega_derivative_upper_integers":[max(ceiling_positive(factorial(j)*abs_upper(w[j]))
            for w in omega) for j in range(derivatives+1)],
        "frequency_lower_gt_99_over_100":all(w[0].imag.is_zero() and w[0].real>arb("99/100") for w in omega)}

if __name__=="__main__":
    d=build()
    for key,value in d.items():
        if key.endswith("integers") or key.startswith("frequency_lower"):
            print(key,value,flush=True)
