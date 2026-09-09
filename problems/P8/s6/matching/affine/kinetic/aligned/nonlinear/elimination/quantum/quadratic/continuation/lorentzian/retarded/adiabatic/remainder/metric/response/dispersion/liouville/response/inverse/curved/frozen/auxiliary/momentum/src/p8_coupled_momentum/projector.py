"""Exact York inverse and its nonzero-wavevector operator bound."""
from functools import cache

import sympy as sp


@cache
def data():
    k=sp.Matrix(sp.symbols("wave_x wave_y wave_z",real=True))
    k2=(k.T*k)[0]
    operator=-k2*sp.eye(3)-k*k.T/3
    inverse=-(sp.eye(3)-k*k.T/(4*k2))/k2
    lift=sp.Matrix([[k[a]*int(b==c)+k[b]*int(a==c)-sp.Rational(2,3)*int(a==b)*k[c]
                    for c in range(3)] for a in range(3) for b in range(3)])
    norm_squared=2*sp.eye(3)/k2-k*k.T/(2*k2**2)
    return {"wavevector":k,"k_squared":k2,"operator":operator,"inverse":inverse,
            "real_York_lift_without_common_i":lift,
            "full_lift_inverse_norm_squared":norm_squared,
            "inverse_norm_upper":"1/|k|^2",
            "full_York_lift_inverse_norm_upper":"sqrt(2)/|k| <= 2/|k|",
            "domain":"Real nonzero wavevector. Zero proper-subset transfer is excluded only from this finite vertex domain, not from the action."}


@cache
def checks():
    d=data()
    E,R,L,k,k2=(d[name] for name in ("operator","inverse","real_York_lift_without_common_i","wavevector","k_squared"))
    return {"exact_three_component_York_inverse":sp.ImmutableMatrix((E*R-sp.eye(3)).applyfunc(sp.cancel)),
            "exact_York_determinant":sp.factor(E.det()+sp.Rational(4,3)*k2**3),
            "exact_full_tensor_lift_inverse_norm_matrix":sp.ImmutableMatrix(
                (R.T*L.T*L*R-d["full_lift_inverse_norm_squared"]).applyfunc(sp.cancel)),
            "sharp_lift_norm_deficit_is_positive_rank_one":sp.ImmutableMatrix(
                (2*sp.eye(3)/k2-d["full_lift_inverse_norm_squared"]-k*k.T/(2*k2**2)).applyfunc(sp.cancel))}
