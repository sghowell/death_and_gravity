"""Independent scalar heat-kernel geometry on a diagonal auxiliary metric."""
from functools import cache, wraps
from itertools import product

import sympy as sp

from . import tensors

g0, gx = sp.symbols("auxiliary_lapse_squared auxiliary_scale_squared", positive=True)
n10, n01, n20, n11, n02 = sp.symbols("log_lapse_t log_lapse_x log_lapse_tt log_lapse_tx log_lapse_xx", real=True)
v10, v01, v20, v11, v02 = sp.symbols("log_scale_t log_scale_x log_scale_tt log_scale_tx log_scale_xx", real=True)


def _indices(arity):
    def decorate(function):
        cached = cache(function)
        @wraps(function)
        def checked(*args):
            if len(args) != arity or any(type(index) is not int or not 0 <= index < 4 for index in args):
                raise ValueError("Require native four-dimensional coordinate indices")
            return cached(*args)
        return checked
    return decorate


def derivative(value, direction):
    if type(direction) is not int or not 0 <= direction < 4:
        raise ValueError("Require a native four-dimensional coordinate direction")
    if direction == 0:
        flows = {g0: 2*n10*g0, gx: 2*(tensors.H[0]+v10)*gx,
                 tensors.H[0]: tensors.H[1], n10: n20, n01: n11,
                 v10: v20, v01: v11}
    elif direction == 1:
        flows = {g0: 2*n01*g0, gx: 2*v01*gx, n10: n11, n01: n02,
                 v10: v11, v01: v02}
    else:
        return sp.Integer(0)
    return sp.expand(sum(sp.diff(value, variable)*flow for variable, flow in flows.items()))


@_indices(1)
def metric(index):
    return g0 if index == 0 else gx


@_indices(3)
def connection(a, b, c):
    return sp.expand(((derivative(metric(a), b) if a == c else 0)
                      +(derivative(metric(a), c) if a == b else 0)
                      -(derivative(metric(b), a) if b == c else 0))/(2*metric(a)))


@_indices(4)
def curvature(a, b, c, d):
    return sp.expand(derivative(connection(a, d, b), c)-derivative(connection(a, c, b), d)
                     +sum(connection(a, c, e)*connection(e, d, b)
                          -connection(a, d, e)*connection(e, c, b) for e in range(4)))


@_indices(2)
def ricci(a, b):
    return sp.expand(sum(curvature(c, a, c, b) for c in range(4)))


@cache
def heat_kernel():
    scalar = sp.expand(sum(ricci(a, a)/metric(a) for a in range(4)))
    ricci_squared = sp.expand(sum(ricci(a, b)**2/(metric(a)*metric(b))
                                  for a, b in product(range(4), repeat=2)))
    riemann_squared = sp.expand(sum(metric(a)*curvature(a, b, c, d)**2
                                    /(metric(b)*metric(c)*metric(d))
                                    for a, b, c, d in product(range(4), repeat=4)))
    # Full fixed four-dimensional scalar coefficient, modulo the covariant
    # total Laplacian only. No Gauss-Bonnet identity is used in this basis.
    scalar_heat = sp.expand((riemann_squared-ricci_squared)/180+scalar**2/72)
    return {"scalar": scalar, "ricci_squared": ricci_squared,
            "riemann_squared": riemann_squared, "scalar_heat": scalar_heat}


@cache
def checks():
    h0, h1 = tensors.H[:2]
    zero = {value: 0 for value in (n10, n01, n20, n11, n02, v10, v01, v20, v11, v02)}
    zero.update({g0: 1, gx: 1})
    data = heat_kernel()
    rate, rate_derivative = h0+v10, h1+v20
    temporal_curvature = rate_derivative+rate*(3*rate-n10)
    expected = {(0, 0): -3*(rate_derivative+rate**2-n10*rate)-g0*(n02+n01**2+n01*v01)/gx,
                (0, 1): 2*(n01*rate-v11),
                (1, 1): -gx*temporal_curvature/g0-n02-n01**2+n01*v01-2*v02,
                (2, 2): -gx*temporal_curvature/g0-v02-v01**2-n01*v01}
    out = {"coordinate_R00": sp.expand(ricci(0, 0).subs(zero)+3*(h1+h0**2)),
            "coordinate_R11": sp.expand(ricci(1, 1).subs(zero)+h1+3*h0**2),
            "coordinate_scalar": sp.expand(data["scalar"].subs(zero)-tensors.scalar_curvature()),
            "coordinate_Riemann_square": sp.expand(data["riemann_squared"].subs(zero)
                                                    -12*((h1+h0**2)**2+h0**4)),
            "coordinate_Ricci_symmetry": sp.expand(ricci(0, 1)-ricci(1, 0))}
    out.update({"general_coordinate_Ricci_"+str(a)+str(b): sp.expand(ricci(a, b)-value)
                for (a, b), value in expected.items()})
    return out
