"""Exact finite time jets of the ACTUAL global all-order Riccati induction."""
from functools import cache

import sympy as sp
from p8_proca_global_hadamard import canonical
from p8_proca_global_hadamard import jets as parent_jets
from p8_proca_global_support import model

FIELD=sp.QQ.algebraic_field(sp.sqrt(17985),sp.I)


class TimeJets:
    def __init__(self,degree):
        if isinstance(degree,bool) or not isinstance(degree,int) or degree<0:
            raise ValueError("Require a nonnegative exact integer time-jet degree")
        self.degree=degree
        self.zero=tuple(FIELD.zero for _ in range(degree+1))
        self.one=(FIELD.one,)+self.zero[1:]

    def constant(self,value):
        return (FIELD.convert(value),)+self.zero[1:]

    def add(self,*values):
        return tuple(sum((value[j] for value in values),FIELD.zero) for j in range(self.degree+1))

    def scale(self,value,scalar):
        scalar=FIELD.convert(scalar)
        return tuple(scalar*v for v in value)

    def mul(self,left,right):
        return tuple(sum((left[j]*right[n-j] for j in range(n+1)),FIELD.zero)
            for n in range(self.degree+1))

    def inv(self,value):
        if value[0]==FIELD.zero:
            raise ValueError("Cannot invert a vanishing time-jet pivot")
        out=[FIELD.one/value[0]]
        for n in range(1,self.degree+1):
            out.append(-sum((value[j]*out[n-j] for j in range(1,n+1)),FIELD.zero)/value[0])
        return tuple(out)

    def div(self,left,right):
        return self.mul(left,self.inv(right))

    def derivative(self,value):
        return tuple(FIELD.convert(j+1)*value[j+1] for j in range(self.degree))+(FIELD.zero,)

    def truncate(self,value,degree):
        return tuple(value[j] if j<=degree else FIELD.zero for j in range(self.degree+1))

    def rational(self,value):
        num,den=sp.fraction(sp.cancel(value))
        pn,pd=sp.Poly(num,model.u),sp.Poly(den,model.u)
        n=tuple(FIELD.convert(pn.nth(j)) for j in range(self.degree+1))
        d=tuple(FIELD.convert(pd.nth(j)) for j in range(self.degree+1))
        return self.div(n,d)

    def sqrt(self,value,root):
        root=FIELD.convert(root)
        if root*root!=value[0] or root==FIELD.zero:
            raise ValueError("The supplied root does not match the nonzero constant jet")
        out=[root]
        for n in range(1,self.degree+1):
            out.append((value[n]-sum((out[j]*out[n-j] for j in range(1,n)),FIELD.zero))/(2*root))
        return tuple(out)

    def matrix(self,rows):
        return tuple(tuple(rows[i][j] for j in range(2)) for i in range(2))

    def mzero(self):
        return self.matrix([[self.zero,self.zero],[self.zero,self.zero]])

    def madd(self,*values):
        return self.matrix([[self.add(*(v[i][j] for v in values)) for j in range(2)] for i in range(2)])

    def mmul(self,left,right):
        return self.matrix([[self.add(*(self.mul(left[i][r],right[r][j]) for r in range(2)))
            for j in range(2)] for i in range(2)])

    def transpose(self,value):
        return self.matrix([[value[j][i] for j in range(2)] for i in range(2)])

    def mderivative(self,value):
        return self.matrix([[self.derivative(value[i][j]) for j in range(2)] for i in range(2)])

    def mtruncate(self,value,degree):
        return self.matrix([[self.truncate(value[i][j],degree) for j in range(2)] for i in range(2)])


@cache
def build(degree=6):
    T=TimeJets(degree)
    bg=model.background()
    ell=T.rational(bg["ell"])
    a=T.rational(bg["a"])
    mass=T.rational(2*bg["J_new"]/bg["lam"]**2)
    c=T.sqrt(T.rational(bg["clock_physical_squared_speed"]),sp.sqrt(17985)/135)
    omega=(T.div(c,a),T.inv(a))
    Vi=T.matrix([[T.one,T.zero],[ell,T.one]])
    V=T.matrix([[T.one,T.zero],[T.scale(ell,-1),T.one]])
    B=T.matrix([[T.add(T.mul(mass,omega[0]),T.div(T.mul(ell,ell),a)),T.div(ell,a)],
        [T.div(ell,a),T.inv(a)]])
    R0=T.matrix([[T.scale(B[i][j],-sp.I) for j in range(2)] for i in range(2)])
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

    def coefficient(n,rs):
        forcing=T.mderivative(rs[n-1]) if n>=1 else T.mzero()
        if n<len(Hqq):
            forcing=T.madd(forcing,Hqq[n])
        for r in range(len(hs)):
            for p in range(len(rs)):
                if p+r==n:
                    forcing=T.madd(forcing,T.mmul(rs[p],Hpq[r]),T.mmul(Hqp[r],rs[p]))
                for s in range(len(rs)):
                    if p+r+s==n:
                        forcing=T.madd(forcing,T.mmul(T.mmul(rs[p],Hpp[r]),rs[s]))
        return forcing

    rs=[R0]
    for n in range(1,degree+1):
        forcing=coefficient(n,rs)
        transformed=T.mmul(T.mmul(T.transpose(V),forcing),V)
        solved=T.matrix([[T.scale(T.div(transformed[i][j],T.add(omega[i],omega[j])),-sp.I)
            for j in range(2)] for i in range(2)])
        result=T.mmul(T.mmul(T.transpose(Vi),solved),Vi)
        rs.append(T.mtruncate(result,degree-n))

    for n in range(degree+1):
        residual=T.mtruncate(coefficient(n,rs),degree-n)
        assert all(value==FIELD.zero for row in residual for jet in row for value in jet),n
        assert rs[n][0][1]==rs[n][1][0],n
        assert all(rs[n][i][j][p]==FIELD.zero for i in range(2) for j in range(2)
            for p in range(degree-n+1) if (p+n)%2),n
    if degree>=2:
        prior=parent_jets.data()["center_Riccati_second_coefficient"]
        assert all(rs[2][i][j][0]==FIELD.convert(prior[i,j]) for i in range(2) for j in range(2))
    return {"degree":degree,"jets":rs,"at_center":[sp.Matrix(2,2,lambda i,j,n=n:FIELD.to_sympy(rs[n][i][j][0]))
        for n in range(degree+1)]}
