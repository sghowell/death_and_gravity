"""Differentiated ordinary-energy subtraction and independently checked tails."""
from functools import cache

import sympy as sp
from p8_vector_hadamard import series
from p8_vector_regularity import spectral
from p8_vector_state import wkb
from p8_vector_subtraction import tail

from . import readouts


def adiabatic_rows(kind,order):
    readouts.validate(kind,order)
    return _adiabatic_rows(kind,order)


@cache
def _adiabatic_rows(kind,order):
    generic=tail.algebra()
    data=wkb.frequency(kind)
    mapping={tail.wa:1,tail.wb:1,tail.c:readouts.weights(kind)["c1"],
             tail.p:data["P2"],tail.r:data["P4"],tail.b2:data["B2"],tail.b4:data["B4"]}
    H,lam=wkb.background()["H"],wkb.background()["lambda"]
    result={}
    for n,name in enumerate(("zero","second","fourth")):
        value=sp.factor(generic[name].subs(mapping,simultaneous=True))
        for _ in range(order):
            value=sp.factor(spectral.D0(value)+((1-2*n)*lam-3*H)*value)
        result[n]=value
    return result


def reference_tail(kind,order,wkb_order=4):
    readouts.validate(kind,order)
    if type(wkb_order) is not int or not 2<=wkb_order<=4:
        raise ValueError("Require native reference WKB coefficient order 2 through 4")
    return _reference_tail(kind,order,wkb_order)


@cache
def _reference_tail(kind,order,wkb_order):
    data=readouts.row_bounds(kind,order)
    lam=wkb.background()["lambda"]
    S,R,Sb,Rb={0:sp.Integer(1)},{},{0:sp.Integer(1)},{}
    for n in range(1,wkb_order+1):
        P=series.coefficient(kind,n)
        B=sp.factor(spectral.D0(P)-2*n*lam*P)
        S[n],R[n]=P,B
        Sb[n],Rb[n]=wkb.box_bound(P)["absolute_upper"],wkb.box_bound(B)["absolute_upper"]
    adiabatic=adiabatic_rows(kind,order)
    polynomial=spectral.numerator(S,R,data["coefficients_in_t"],readouts.weights(kind)["c1"],adiabatic)
    low={n:sp.factor(polynomial.get(n,0)) for n in range(5)}
    upper=spectral.numerator(Sb,Rb,data["coefficient_majorants"],sp.Integer(2),
        {n:wkb.box_bound(value)["absolute_upper"] for n,value in adiabatic.items()},absolute=True)
    mass=wkb.MASS_TIME_MIN
    bound=2*sum(value/mass**(2*(n-5)) for n,value in upper.items() if n>=5)
    return {"low_coefficient_residuals":low,
            "integrable_tail_certified":all(value==0 for value in low.values()),
            "reference_bracket_tail_over_t_cubed_upper":bound,
            "reference_bracket_tail_integer_upper":sp.ceiling(bound),
            "no_lower_Laurent_power":min(polynomial)>=0,
            "all_majorant_coefficients_nonnegative":all(value>=0 for value in upper.values())}


@cache
def checks():
    H=wkb.background()["H"]
    out={}
    for kind in ("transverse","longitudinal"):
        for n,name in enumerate(("zero","second","fourth")):
            pressure=tail.reference_terms()[kind+"_pressure"][name]
            out[kind+f"_actual_adiabatic_Ward_conservation_{2*n}"]=sp.factor(
                adiabatic_rows(kind,1)[n]+3*H*(adiabatic_rows(kind,0)[n]+pressure))
        for order in range(6):
            for n,value in reference_tail(kind,order)["low_coefficient_residuals"].items():
                out[kind+f"_ordinary_energy_time{order}_nonintegrable_coefficient_{n}"]=value
    failure=reference_tail("longitudinal",4,2)["low_coefficient_residuals"][4]
    out["fourth_order_low_reference_failure_control"]=sp.factor(failure.subs({readouts.u:0,readouts.z:1})-23808)
    return out
