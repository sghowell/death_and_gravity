"""Complete Euler systems and uniform complex-large-frequency chart remainders."""
from functools import cache

import sympy as sp

from . import model

o=model.old


def clean(value):
    return value.applyfunc(sp.factor) if isinstance(value,sp.MatrixBase) else sp.factor(value)


def limit_matrix(matrix,power=0):
    return clean(matrix.applyfunc(lambda v:sp.limit(v/o.q**power,o.q,sp.oo)))


@cache
def data():
    lapse=model.margin.scalar()["lapse"]
    L=(o.action()["base"]+model.margin.delta_J*o.n**2).subs(o.n,lapse)
    x,xd=sp.Matrix([o.v,o.matter]),sp.Matrix([o.vd,o.sd])
    alpha=clean(sp.hessian(L,tuple(xd)))
    beta=clean(sp.Matrix(2,2,lambda i,j:sp.diff(L,xd[i],x[j])))
    gamma=clean(sp.hessian(L,tuple(x)))
    g=model.margin.gamma()
    raw={"unitary":(alpha,beta,gamma),
        "gamma":tuple(clean(g[name]/o.a**3) for name in ("alpha","beta","gamma"))}
    out={}
    for name,(A,B,C) in raw.items():
        K=limit_matrix(A)/2
        Bq,Cq=limit_matrix(B,1),limit_matrix(C,1)
        G=clean((model.dt(Bq)+o.H*Bq-Cq)/2)
        friction=clean(A.inv()*(model.dt(A)+3*o.H*A+B-B.T))
        force=clean(A.inv()*(model.dt(B)+3*o.H*B-C))
        out[name]={"alpha_over_a_cubed":A,"beta_over_a_cubed":B,"gamma_over_a_cubed":C,
            "K_principal":K,"G_principal":G,"beta_q_principal":Bq,
            "complete_velocity_friction":friction,"complete_velocity_force":force,
            "complete_bounded_force_remainder":clean(force-o.q*K.inv()*G),
            "complete_leading_spatial_force":clean(K.inv()*G),
            "finite_q_chart_denominators":sorted({str(sp.factor(sp.denom(sp.cancel(v))))
                for matrix in (A,B,C,friction,force) for v in matrix})}
    return out


@cache
def checks():
    out={}
    for name,d in data().items():
        K,G=d["K_principal"],d["G_principal"]
        den=o.theta if name=="unitary" else o.lam
        sign=1 if name=="unitary" else -1
        K0=sp.Matrix([[(o.J+o.w**2/2)/den**2,sign*o.w/(2*den)],
            [sign*o.w/(2*den),sp.Rational(1,2)]])
        out[name+"_actual_complete_Euler_gradient_is_old_positive_form"]=clean((G-K0).subs(model.substitution(),simultaneous=True))
        out[name+"_retuned_principal_kinetic_adds_only_actual_margin"]=clean(K-K0-sp.diag(model.margin.delta_J/den**2,0))
        out[name+"_leading_velocity_mixing_is_symmetric"]=clean(d["beta_q_principal"]-d["beta_q_principal"].T)
        out[name+"_complete_friction_has_no_positive_q_growth"]=limit_matrix(d["complete_velocity_friction"],1)
        out[name+"_complete_force_remainder_has_no_positive_q_growth"]=limit_matrix(d["complete_bounded_force_remainder"],1)
    out["unitary_complete_force_remainder_is_exactly_q_independent"]=clean(
        data()["unitary"]["complete_bounded_force_remainder"].diff(o.q))
    out["unitary_complete_friction_is_exactly_q_independent"]=clean(
        data()["unitary"]["complete_velocity_friction"].diff(o.q))
    return out
