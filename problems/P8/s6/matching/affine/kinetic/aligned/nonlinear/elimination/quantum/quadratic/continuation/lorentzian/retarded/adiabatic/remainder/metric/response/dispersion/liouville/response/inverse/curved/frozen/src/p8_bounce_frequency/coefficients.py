"""Freeze the physical Euler operator before a constant source congruence."""
import json
from functools import cache

import sympy as sp
from p8_affine import verify as affine
from p8_prepared_volterra import tree
from p8_vector_metric_dispersion import principal, spectral
from p8_vector_metric_local import verify as local
from p8_vector_metric_response import tadpole
from p8_vector_metric_response import verify as response

LOCAL_SHA="4ab3213760790f7309b69596c347d21c2a20552f7d010a12b500f6552b1cbaa0"
RESPONSE_SHA="b67525a295507e7c6e9a05a01c5f51e36abca7f004bd8d2b570bcb81bb9aabe7"
L=sp.Integer(10)**24
MASS=sp.Integer(1000)
MARGIN=sp.Rational(1,10**6)
T=sp.ImmutableMatrix([[1,0],[sp.Rational(1,2),1]])
mass=sp.Symbol("constant_reference_mass",positive=True)


def degree(value):
    if type(value) is not int or value not in range(5):
        raise ValueError("Require native derivative order 0 through 4")
    return value


def pinned(path,expected):
    if affine.sha(path)!=expected:
        raise ValueError("A frozen diagnostic input report changed")
    return json.loads(path.read_text())


@cache
def background_jets():
    old=tree.system.old
    bg=old.background()
    actual={tree.H:bg["H"],tree.theta:bg["theta"],tree.ell:bg["ell"],
            tree.w:bg["w"],tree.Je:bg["J"]+4*MARGIN/(1+old.u**2)**6,
            tree.delta:bg["delta"]}
    result={}
    for field,expression in actual.items():
        for j in range(3):
            result[sp.diff(field,tree.u,j)]=sp.diff(expression,old.u,j).subs(old.u,0)
    return result


@cache
def tree_rows():
    physical=tree.data()["physical_metric_Euler_currents"]
    return tuple(sp.ImmutableMatrix((T.T*sp.Matrix([
        [sp.diff(row,sp.diff(field,tree.u,j)) for field in (tree.n,tree.zeta)]
        for row in physical]).xreplace(background_jets())*T).applyfunc(sp.factor))
        for j in range(3))


def local_row(value):
    return _local_row(degree(value))


@cache
def _local_row(value):
    rows=pinned(local.REPORT,LOCAL_SHA)["finite_Euler_current_coefficient_rows"]
    u=sp.Symbol("u",real=True)
    C=sp.Matrix([[sum(sp.sympify(rows[str(order)][left+right][value],
        locals={"u":u,"constant_reference_mass":mass}) for order in range(3))
        for right in ("N","Z")] for left in ("N","Z")])
    return sp.ImmutableMatrix((T.T*C.subs(u,0)*T).applyfunc(sp.factor))


@cache
def profile_box():
    data=pinned(response.REPORT,RESPONSE_SHA)["complete_background_cancelled_response_scale_example"]
    return {key:sp.Rational(value) for key,value in data["fixed_profile_C0_upper_bounds"].items()}


@cache
def tadpole_matrix():
    data=tadpole.physical_vertices()
    return sp.ImmutableMatrix((T.T*data["density_hessian"]*T).applyfunc(sp.expand))


@cache
def pairs():
    return {sector:sp.ImmutableMatrix((T.T*spectral.actual(spectral.pair(sector)).subs(spectral.h,1)).applyfunc(sp.factor))
            for sector in ("T","L")}


@cache
def gram():
    pair=pairs()
    return sp.ImmutableMatrix((2*pair["T"]*pair["T"].T+pair["L"]*pair["L"].T).applyfunc(sp.factor))


@cache
def finite():
    return sp.ImmutableMatrix((T.T*principal.finite().subs(spectral.h,1)*T).applyfunc(sp.factor))


@cache
def checks():
    m=mass
    expected_zero=sp.Matrix([
        [(697897*m**4+855392*m**2+10020608)/52488,
         (6867*m**4+39776*m**2-84864)/324],
        [(763*m**4+1760*m**2+8960)/36,3*(15*m**4+160*m**2-384)/2]])
    expected_two=sp.Matrix([
        [(35917*m**2-7740)/6561,2*(1243*m**2-8676)/243],
        [2*(1243*m**2-972)/243,4*(5*m**2-36)]])
    expected_four=sp.Matrix([[-sp.Rational(7357,6561),-sp.Rational(562,243)],
                             [-sp.Rational(562,243),-4]])
    expected_tree=(sp.Matrix([[sp.Rational(749377,250000),-sp.Rational(3,200)],
                             [-sp.Rational(3,200),-sp.Rational(9,100)]]),
                   sp.zeros(2),sp.diag(0,6))
    out={"constant_chart_determinant":T.det()-1}
    for j in range(3):
        out["actual_physical_tree_frozen_coefficient_"+str(j)]=tree_rows()[j]-expected_tree[j]
    for j,target in enumerate((expected_zero,sp.zeros(2),expected_two,sp.zeros(2),expected_four)):
        out["actual_finite_local_frozen_coefficient_"+str(j)]=sp.ImmutableMatrix((local_row(j)-target).applyfunc(sp.factor))
    out["frozen_fourth_contact_equals_dispersion_contact"]=local_row(4)-finite()
    data=tadpole.physical_vertices()
    rho,p=data["rho"],data["pressure"]
    out["fixed_tadpole_constant_chart_Hessian"]=tadpole_matrix()-sp.Matrix(
        [[2*rho-sp.Rational(13,4)*p,3*rho-sp.Rational(9,2)*p],
         [3*rho-sp.Rational(9,2)*p,-9*p]])
    z=spectral.z
    for sector,target in {
        "T":sp.Matrix([sp.Rational(109,81)*(1-z),2*(1-z)]),
        "L":sp.Matrix([(109+117*z)/81,2*(1+z)])}.items():
        out[sector+"_actual_pair_in_constant_chart"]=pairs()[sector]-target
    out["positive_pair_Gram_replays_actual_spectral_matrix"]=sp.ImmutableMatrix(
        (gram()-T.T*spectral.actual(spectral.matrix()).subs(spectral.h,1)*T).applyfunc(sp.factor))
    out["curved_bounce_jets_not_set_to_zero"]=sp.diff(tree.theta,tree.u).xreplace(background_jets())-3
    out["time_dependent_chart_second_jet_retained"]=sp.diff(tree.delta,tree.u,2).xreplace(background_jets())+3
    return {key:sp.ImmutableMatrix(value.applyfunc(sp.factor)) if isinstance(value,sp.MatrixBase)
            else sp.factor(value) for key,value in out.items()}
