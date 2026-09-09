"""Exact rational compact-interval and state-center bounds."""
from fractions import Fraction as F
from functools import cache
from math import isqrt

import sympy as sp
from p8_proca_global_support import model

from .canonical import data as regular_data
from .time_jets import build


def fraction(x):
    x=sp.Rational(x)
    return F(int(x.p),int(x.q))

def poly_range(p,radius):
    lo=hi=F(0)
    for (j,),c in sp.Poly(p,model.u).terms():
        c=fraction(c)
        if j==0:
            lo+=c
            hi+=c
        elif j%2:
            d=abs(c)*radius**j
            lo-=d
            hi+=d
        elif c>0:
            hi+=c*radius**j
        else:
            lo+=c*radius**j
    return lo,hi

def rational_range(expr,radius):
    if radius<=0:
        raise ValueError("Positive rational radius required")
    n,d=sp.fraction(sp.cancel(expr))
    a,b=poly_range(n,radius)
    c,e=poly_range(d,radius)
    if c<=0<=e:
        raise ValueError("Denominator enclosure includes zero")
    vals=(a/c,a/e,b/c,b/e)
    return min(vals),max(vals)

def absv(interval):
    return max(abs(v) for v in interval)

def matrix_norm(M,r):
    return max(sum(absv(rational_range(M[i,j],r)) for j in range(M.cols)) for i in range(M.rows))

def gershgorin(M,r):
    return min(rational_range(M[i,i],r)[0]-
        sum(absv(rational_range(M[i,j],r)) for j in range(M.cols) if j!=i) for i in range(M.rows))

def ceil(x):
    return -((-x.numerator)//x.denominator)

@cache
def center_bounds():
    d=build(8)
    scale=10**50
    root=F(isqrt(17985*scale**2),scale)
    assert root**2<17985<(root+F(1,scale))**2
    intervals={}
    for n in (0,2,4,6,8):
        vals=[]
        for x in d["at_center"][n]/sp.I:
            p=sp.Poly(sp.expand(x),sp.sqrt(17985))
            assert p.degree()<=1
            a,b=fraction(p.nth(0)),fraction(p.nth(1))
            ends=(a+b*root,a+b*(root+F(1,scale)))
            vals.append((min(ends),max(ends)))
        intervals[n]=tuple(vals)
    bounds={n:ceil(max(sum(absv(intervals[n][2*i+j]) for j in range(2)) for i in range(2)))
        for n in (2,4,6,8)}
    assert bounds[2]<32 and bounds[4]<506 and bounds[6]<9000 and bounds[8]<600000
    B0=-d["at_center"][0]/sp.I
    assert B0[1,1]-abs(B0[1,0])==sp.Rational(9,10)
    # The first Gershgorin row is bounded with the same exact root enclosure.
    assert -intervals[0][0][1]-absv(intervals[0][1])>F(9,10)
    norm0=max(sum(absv(intervals[0][2*i+j]) for j in range(2)) for i in range(2))
    assert norm0<F(64,5)
    total=F(32,16**2)+F(506,16**4)+F(9000,16**6)+F(600001,16**8)
    assert total<F(3,20)
    total64=F(32,64**2)+F(506,64**4)+F(9000,64**6)+F(600001,64**8)
    assert total64<F(1,100)
    return {"root_isolation":(root,root+F(1,scale)),"coefficient_norm_upper_integers":bounds,
        "tail_budget_at_16":total,"tail_budget_at_64":total64,
        "positive_graph_threshold":16,"positive_graph_lower":F(3,4),"positive_graph_upper":13}

@cache
def transfer_bounds(radius=F(1,100)):
    d=regular_data()
    omega=sp.zeros(4)
    omega[:2,2:]=sp.eye(2)
    omega[2:,:2]=-sp.eye(2)
    H=-omega*d["actual_coefficients"][1]
    l0=gershgorin(H,radius)
    h0=matrix_norm(H,radius)
    h1=matrix_norm(H.diff(model.u),radius)
    # Each real Hamiltonian Laurent term has the same 1 and infinity norm
    # up to conjugation by a signed permutation, but use both for safety.
    b=sum(max(matrix_norm(d["actual_coefficients"][j],radius),
        matrix_norm(d["actual_coefficients"][j].T,radius)) for j in range(-3,1))
    assert l0>0
    lower=F(1,ceil(1/l0))
    upper=ceil(h0)
    derivative=ceil(h1)
    remainder=ceil(b)
    rate=ceil(F(derivative+2*upper*remainder,2*lower))
    return {"radius":radius,"energy_lower":lower,"energy_upper":upper,
        "energy_derivative_norm_upper":derivative,"remainder_norm_upper":remainder,
        "norm_log_rate":rate,"exact_raw":(l0,h0,h1,b)}

if __name__=="__main__":
    print("CENTER",center_bounds(),flush=True)
    d=transfer_bounds()
    print("TRANSFER",{k:v for k,v in d.items() if k!="exact_raw"},flush=True)
