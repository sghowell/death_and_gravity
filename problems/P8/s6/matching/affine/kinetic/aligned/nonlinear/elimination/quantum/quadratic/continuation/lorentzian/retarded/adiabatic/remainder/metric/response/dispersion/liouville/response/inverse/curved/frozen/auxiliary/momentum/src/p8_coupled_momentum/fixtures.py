"""Mixed scalar, tensor, matter and vector fixtures at exact symbolic time."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_physical import geometry
from p8_physical import jets as j

from . import constraints


def labels(value):
    if type(value) is not int or value not in (2,3,4):
        raise ValueError("Require native labelled degree two, three or four")
    return value


def build(value):
    return _build(labels(value))


@cache
def _build(value):
    momenta={2:((1,2,1),(-1,-2,-1)),
             3:((1,0,0),(0,1,0),(-1,-1,0)),
             4:((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))}[value]
    ctx=constraints.context(momenta,symbolic_time=True,symbolic_scale=True)
    bg=model.coefficients()["background"]
    fields={name:ctx.jet() for name in ("zeta","scalar_p","chi","chi_p")}
    tensor,tensor_p=j.zeros(ctx),j.zeros(ctx)
    vector,vector_p=([ctx.jet() for _ in range(3)] for _ in range(2))
    polarizations=[]
    for n,wave in enumerate(momenta):
        k=sp.Matrix(wave)
        # A fixed diagonal seed can project to zero for an allowed k.
        # Use two nonvanishing rational tangents, with no square roots.
        tangent=next(k.cross(sp.eye(3)[:,a]) for a in range(3)
                     if k.cross(sp.eye(3)[:,a])!=sp.zeros(3,1))
        other=k.cross(tangent)
        E=tangent*tangent.T/(tangent.T*tangent)[0]-other*other.T/(other.T*other)[0]
        norm=sp.trace(E*E)
        polarizations.append(E)
        fields["zeta"]+=ctx.leg(n,sp.Rational(n+1,100))
        fields["scalar_p"]+=ctx.leg(n,sp.Rational(n+1,11))
        fields["chi"]+=ctx.leg(n,sp.Rational(n+2,19))
        fields["chi_p"]+=ctx.leg(n,sp.Rational(n+3,23))
        for a in range(3):
            vector[a]+=ctx.leg(n,sp.Rational(n+a+2,29))
            vector_p[a]+=ctx.leg(n,sp.Rational((n+1)*(a+1)+1,37))
            for b in range(3):
                tensor[a][b]+=ctx.leg(n,E[a,b]*sp.Rational(n+1,31))
                tensor_p[a][b]+=ctx.leg(n,E[a,b]*sp.Rational(n+1,41)/norm)
    geo=geometry.derive(ctx,fields["zeta"],tensor)
    result=constraints.derive(ctx,geo,fields["zeta"],tensor,fields["scalar_p"],tensor_p,
        fields["chi"],fields["chi_p"],vector,vector_p,bg["H"],bg["ell"])
    return {"context":ctx,"geometry":geo,"fields":fields,"tensor":tensor,"tensor_p":tensor_p,
            "vector":vector,"vector_p":vector_p,"polarizations":polarizations,"result":result}


@cache
def checks():
    out={}
    bg=model.coefficients()["background"]
    for n in (2,3,4):
        data=build(n)
        ctx,geo,result=data["context"],data["geometry"],data["result"]
        rest=constraints.residual(result["momentum"],geo,result["physical_lower_current"])
        for degree in range(1,n):
            masks=[mask for mask in range(1,ctx.full) if mask.bit_count()==degree]
            out[f"mixed_{n}_leg_all_three_constraints_degree_{degree}"]=sp.ImmutableMatrix(
                [[rest[a].coefficient(mask) for mask in masks] for a in range(3)])
        target=constraints.linear_target(ctx,data["fields"]["scalar_p"],data["fields"]["chi"],bg["ell"])
        out[f"mixed_{n}_leg_exact_linear_matter_shift"]=sp.ImmutableMatrix(
            [[(result["vector_orders"][0][a]-target[a]).coefficient(1<<leg) for leg in range(n)]
             for a in range(3)])
        for leg,E in enumerate(data["polarizations"]):
            k=sp.Matrix(ctx.momenta[leg])
            out[f"mixed_{n}_leg_TT_trace_{leg}"]=sp.trace(E)
            out[f"mixed_{n}_leg_TT_divergence_{leg}"]=sp.ImmutableMatrix(E*k)
    return out


@cache
def negative_controls():
    data=build(4)
    ctx,geo=data["context"],data["geometry"]
    divergence=sum(data["vector_p"][a].derivative(a) for a in range(3))
    missing=[sum(geo["inverse"][a][b]*data["vector"][b]*divergence for b in range(3))/2
             for a in range(3)]
    masks=[mask for mask in range(1,ctx.full) if mask.bit_count()==2]
    matrix=sp.ImmutableMatrix([[value.coefficient(mask) for mask in masks] for value in missing])
    bg=model.coefficients()["background"]
    actual=data["result"]["vector_orders"][0]
    no_matter=constraints.linear_target(ctx,data["fields"]["scalar_p"],ctx.jet(),bg["ell"])
    matter_defect=sp.ImmutableMatrix([[(actual[a]-no_matter[a]).coefficient(1<<leg)
                                      for leg in range(4)] for a in range(3)])
    return {"omitted_vector_divergence_constraint_defect":matrix,
            "omitted_background_matter_linear_constraint_solution_defect":matter_defect}
