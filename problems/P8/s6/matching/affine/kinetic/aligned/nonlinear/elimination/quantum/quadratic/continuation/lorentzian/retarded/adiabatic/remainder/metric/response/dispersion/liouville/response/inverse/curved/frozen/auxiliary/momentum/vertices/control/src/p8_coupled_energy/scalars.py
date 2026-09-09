"""Two-chart scalar Lagrangian energy, retaining antisymmetric mixing."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_m1_physical import quadratic as old

z=sp.Symbol("inverse_squared_physical_momentum",nonnegative=True)
H,ell,theta,lam,Je=old.H,old.l,old.theta,old.lam,old.J
BASE=(H,ell,theta,lam,Je)
FIRST=sp.symbols("H_first ell_first Theta_first Lambda_first Je_first",real=True)
SECOND=sp.symbols("H_second ell_second Theta_second Lambda_second Je_second",real=True)


def canonical(value):
    return sp.factor(sp.cancel(value))


def derivative(value):
    result=2*H*z*sp.diff(value,z)
    result+=sum(after*sp.diff(value,before) for before,after in zip(BASE+FIRST,FIRST+SECOND))
    return canonical(result)


def clean(matrix):
    return sp.ImmutableMatrix(matrix.applyfunc(canonical))


def chart(value):
    if type(value) is not str or value not in ("unitary","gamma"):
        raise ValueError("Require the declared scalar chart")
    return value


def data(value):
    return _data(chart(value))


@cache
def _data(value):
    original=old.symbolic(value)
    mapping={old.q:1/z,old.w:-ell*lam}
    alpha=clean(original["alpha"].subs(mapping,simultaneous=True))
    beta=clean(original["beta"].subs(mapping,simultaneous=True))
    C=sp.hessian(original["density"],old.Q)
    positive_gamma=clean((C-original["B"].T*original["alpha"]*original["B"]).subs(mapping,simultaneous=True))
    symmetric=clean((beta+beta.T)/2)
    antisymmetric=clean((beta-beta.T)/2)
    V=clean(positive_gamma+symmetric.applyfunc(derivative)+3*H*symmetric)
    leading=clean(V.applyfunc(lambda entry:canonical(z*entry).subs(z,0)))
    remainder=clean(V-leading/z)
    # The leading beta is symmetric, so no q-enhanced antisymmetric source remains.
    beta_leading=clean(beta.applyfunc(lambda entry:canonical(z*entry).subs(z,0)))
    beta_remainder=clean(beta-beta_leading/z)
    for matrix in (alpha,remainder,antisymmetric,beta_remainder):
        for entry in matrix:
            if canonical(entry).as_numer_denom()[1].subs(z,0)==0:
                raise ValueError("An uncancelled high-frequency denominator remains")
    return {"alpha":alpha,"beta":beta,"positive_gamma":positive_gamma,
            "symmetric_boundary":symmetric,"antisymmetric_mixing":antisymmetric,
            "potential":V,"principal_gradient":leading,"potential_remainder":remainder,
            "beta_leading":beta_leading,"beta_remainder":beta_remainder,
            "energy_velocity_coefficient":clean(alpha.applyfunc(derivative)+3*H*alpha),
            "energy_principal_coefficient":clean(leading.applyfunc(derivative)+H*leading),
            "energy_remainder_coefficient":clean(remainder.applyfunc(derivative)+3*H*remainder),
            "energy_mixing_coefficient":clean(antisymmetric.applyfunc(derivative)+3*H*antisymmetric)}


@cache
def actual_background_jets():
    bg=model.coefficients()["background"]
    actual=(bg["H"],bg["ell"],bg["theta"],bg["lam"],bg["J"]+4*model.MARGIN/bg["h"]**2)
    return {symbol:sp.factor(sp.diff(value,model.u,degree))
            for degree,symbols in enumerate((BASE,FIRST,SECOND))
            for symbol,value in zip(symbols,actual)}


@cache
def checks():
    actual=actual_background_jets()
    bg=model.coefficients()["background"]
    out={}
    for name in ("unitary","gamma"):
        d=data(name)
        denominator=theta if name=="unitary" else lam
        cross=-ell*lam/theta if name=="unitary" else ell
        expected=sp.Matrix([[2*bg["J"]/actual[denominator]**2+actual.get(cross,cross.subs(actual))**2,
                             actual.get(cross,cross.subs(actual))],
                            [actual.get(cross,cross.subs(actual)),1]])
        out[name+"_actual_principal_gradient"]=clean(
            d["principal_gradient"].subs(actual,simultaneous=True)-expected)
        out[name+"_leading_beta_symmetric"]=clean(d["beta_leading"]-d["beta_leading"].T)
        out[name+"_mixing_antisymmetric"]=clean(d["antisymmetric_mixing"]+d["antisymmetric_mixing"].T)
    udata=data("unitary")
    gdata=data("gamma")
    alpha0=clean(gdata["alpha"].subs(z,0))
    J0=Je+ell**2*lam**2/2
    D=lam**2-J0*z
    v=sp.Matrix((sp.sqrt(2)*J0/lam,ell*lam/sp.sqrt(2)))
    out["gamma_finite_q_positive_rank_one_kinetic_increment"]=clean(
        gdata["alpha"]-alpha0-z*(v*v.T)/D)
    vec=sp.Matrix(sp.symbols("coordinate_1 coordinate_2",real=True))
    for name,d in (("unitary",udata),("gamma",gdata)):
        den=theta if name=="unitary" else lam
        cross=-ell*lam/theta if name=="unitary" else ell
        principal=clean(d["alpha"].subs(z,0))
        out[name+"_principal_kinetic_completed_square"]=canonical(
            (vec.T*principal*vec)[0]-(vec[1]+cross*vec[0])**2-2*Je*vec[0]**2/den**2)
    return out
