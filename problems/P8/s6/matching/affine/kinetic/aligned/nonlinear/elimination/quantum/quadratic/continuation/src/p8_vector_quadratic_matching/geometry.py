"""Covariant fixed-coefficient continuation of the scalar auxiliary geometry."""
from functools import cache

import sympy as sp
from p8_vector_quadratic import geometry as four
from p8_vector_quadratic import tensors

dimension = sp.Symbol("spatial_dimension", real=True)


@cache
def heat_kernel():
    d = dimension
    g0, gx = four.g0, four.gx
    h0, h1 = tensors.H[:2]
    rate, acceleration = h0+four.v10, h1+four.v20
    n0, n1, n11, v1, v11 = four.n10, four.n01, four.n02, four.v01, four.v02
    temporal = acceleration+rate*(d*rate-n0)
    ricci00 = -d*(acceleration+rate**2-n0*rate)-g0*(n11+n1**2+(d-2)*n1*v1)/gx
    ricci01 = (d-1)*(n1*rate-four.v11)
    ricci11 = -gx*temporal/g0-n11-n1**2+n1*v1-(d-1)*v11
    ricci22 = -gx*temporal/g0-v11-(d-2)*v1**2-n1*v1
    scalar = sp.expand(ricci00/g0+(ricci11+(d-1)*ricci22)/gx)
    ricci_squared = sp.expand(ricci00**2/g0**2+2*ricci01**2/(g0*gx)
                              +(ricci11**2+(d-1)*ricci22**2)/gx**2)
    # Flat (d-1)-dimensional fiber over the (u,x) two-dimensional base.
    base_curvature = -(acceleration+rate**2-n0*rate)/g0-(n11+n1**2-n1*v1)/gx
    hessian00 = acceleration+rate**2-n0*rate+g0*n1*v1/gx
    hessian01 = four.v11-n1*rate
    hessian11 = v11+gx*rate**2/g0
    gradient = rate**2/g0+v1**2/gx
    riemann_squared = sp.expand(4*base_curvature**2+4*(d-1)*(hessian00**2/g0**2
                                           +2*hessian01**2/(g0*gx)+hessian11**2/gx**2)
                                +2*(d-1)*(d-2)*gradient**2)
    scalar_heat = sp.expand((riemann_squared-ricci_squared)/180+scalar**2/72)
    return {"scalar": scalar, "ricci_squared": ricci_squared, "riemann_squared": riemann_squared,
            "scalar_heat": scalar_heat}


@cache
def checks():
    return {"continued_auxiliary_"+name+"_returns_four_dimensions":
            sp.expand(value.subs(dimension, 3)-four.heat_kernel()[name])
            for name, value in heat_kernel().items()}
