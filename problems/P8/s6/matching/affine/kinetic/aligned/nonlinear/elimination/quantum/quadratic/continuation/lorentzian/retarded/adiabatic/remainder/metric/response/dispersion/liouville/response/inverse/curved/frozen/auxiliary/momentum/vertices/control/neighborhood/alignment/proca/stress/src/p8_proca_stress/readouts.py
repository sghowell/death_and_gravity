"""Actual ordinary Proca energy derivatives in the unchanged clock modes."""
from functools import cache

import sympy as sp
from p8_constant_proca import quantum
from p8_vector_regularity import spectral
from p8_vector_state import wkb
from p8_vector_subtraction import tail
from p8_vector_variation import modes

u,z,omega2=spectral.u,spectral.z,spectral.omega2


def validate(kind,order):
    if type(kind) is not str or kind not in ("transverse","longitudinal"):
        raise ValueError("Require a named physical Proca polarization kind")
    if type(order) is not int or not 0<=order<=5:
        raise ValueError("Require native physical derivative order 0 through 5")


def weights(kind):
    validate(kind,0)
    c=tail.physical_weights()[kind+"_energy"]["c1"]
    return {"mode":kind,"A":sp.Integer(1),"B":sp.Integer(1),"c1":c}


def rows(kind,order):
    validate(kind,order)
    return _rows(kind,order)


@cache
def _rows(kind,order):
    H,lam=wkb.background()["H"],wkb.background()["lambda"]
    d=weights(kind)["c1"]-lam/2
    matrix=sp.Matrix([[2*d,2,0],[-omega2,0,1],[0,-2*omega2,-2*d]])
    row=sp.Matrix([[omega2,0,1]])
    for _ in range(order):
        row=spectral.clean(row.applyfunc(spectral.D)-3*H*row+row*matrix)
    return spectral.clean(row)


def row_bounds(kind,order):
    validate(kind,order)
    return _row_bounds(kind,order)


@cache
def _row_bounds(kind,order):
    coefficients,majorants,reconstructions,degrees=[],[],[],[]
    for column,component in enumerate(rows(kind,order)):
        exact,upper={},{}
        for (power,),value in sp.Poly(component,omega2).terms():
            box=wkb.box_bound(value)
            exact[-power]=value
            upper[-power]=box["absolute_upper"]
            reconstructions.append(box["reconstruction"])
            if value!=0:
                degrees.append(2*power+(-1,0,1)[column]<=order+1)
        coefficients.append(exact)
        majorants.append(upper)
    envelope=sum(sum(poly.values())*(1,2,3)[column] for column,poly in enumerate(majorants))
    return {"coefficients_in_t":coefficients,"coefficient_majorants":majorants,
            "normalized_reference_product_envelope":envelope,
            "all_frequency_powers_at_most_derivative_order_plus_one":all(degrees),
            "box_reconstructions":reconstructions}


@cache
def checks():
    H=wkb.background()["H"]
    out={}
    for kind in ("transverse","longitudinal"):
        out[kind+"_actual_canonical_ordinary_energy_row"]=rows(kind,0)-sp.Matrix([[omega2,0,1]])
        # The unchanged pressure row yields exact ordinary Proca energy
        # conservation in each polarization before momentum integration.
        pressure=spectral.rows(kind+"_pressure",0)
        out[kind+"_actual_mode_energy_Ward_identity"]=spectral.clean(rows(kind,1)+3*H*(rows(kind,0)+pressure))
        weight=tail.physical_weights()[kind+"_pressure"]
        target=sp.diag(weight["B"]*(modes.q+modes.mass2),weight["A"])/modes.scale**3
        target=target.subs(z,modes.q/(modes.q+modes.mass2))
        out[kind+"_actual_new_pressure_is_the_old_covariant_readout"]=spectral.clean(
            quantum.matrices()[kind]["new_canonical_pressure"]-target)
    return out
