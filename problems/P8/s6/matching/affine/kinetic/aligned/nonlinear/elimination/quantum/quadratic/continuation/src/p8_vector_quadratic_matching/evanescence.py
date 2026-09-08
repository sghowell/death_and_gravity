"""Why four-dimensional boundary reductions cannot precede continuation."""
from functools import cache

import sympy as sp
from p8_vector_quadratic import tensors

from . import geometry


@cache
def gauss_bonnet_control():
    d = geometry.dimension
    h, hd = tensors.H[:2]
    zero = {value: 0 for value in (geometry.four.n10, geometry.four.n01, geometry.four.n20,
                                  geometry.four.n11, geometry.four.n02, geometry.four.v10,
                                  geometry.four.v01, geometry.four.v20, geometry.four.v11, geometry.four.v02)}
    zero.update({geometry.four.g0: 1, geometry.four.gx: 1})
    curvatures = geometry.heat_kernel()
    actual = sp.expand((curvatures["riemann_squared"]-4*curvatures["ricci_squared"]
                        +curvatures["scalar"]**2).subs(zero))
    prefactor = d*(d-1)*(d-2)
    expected = prefactor*(4*h**2*hd+(d+1)*h**4)
    divergence = 3*h**2*hd+d*h**4
    remainder = -prefactor*(d-3)*h**4/3
    f0, f1, f2 = sp.symbols("conformal_mass0 conformal_mass1 conformal_mass2", real=True)
    amplitude = sp.Symbol("conformal_amplitude", real=True)
    # After the exact D-dimensional FLRW boundary identity, the derivative
    # of -GB/180 at D=3 is (1/90)[lambda^2](H+v')^4 for a scalar mass.
    rate = h+amplitude*f1/2-amplitude**2*f0*f1/2
    local_dimensional_derivative = sp.expand(rate**4).coeff(amplitude, 2)/90
    first = h**2/60
    zero_coefficient = h**2*(hd+h**2)/30
    derivative_first = h*hd/30
    operator = sp.expand(2*zero_coefficient*f0-2*first*f2-2*(derivative_first+3*h*first)*f1)
    rate_variation = sp.diff(local_dimensional_derivative, f1)
    direct_operator = (sp.diff(local_dimensional_derivative, f0)-hd*sp.diff(rate_variation, h)
                       -f1*sp.diff(rate_variation, f0)-f2*sp.diff(rate_variation, f1)-3*h*rate_variation)
    return {"FLRW_Gauss_Bonnet_D": actual,
            "exact_curvature_identity": sp.expand(actual-expected),
            "exact_boundary_identity": sp.expand(actual-sp.Rational(4, 3)*prefactor*divergence-remainder),
            "topological_factor_has_simple_zero_at_D_three": sp.diff(prefactor*(d-3)/540, d).subs(d, 3),
            "scalar_mass_dimension_derivative_density": local_dimensional_derivative,
            "compact_coefficients": {"first_time_derivative_squared": first, "field_squared": zero_coefficient},
            "dimension_derivative_Euler_operator": operator,
            "scalar_compact_variation_identity": sp.expand(direct_operator-operator),
            "nonzero_dimension_derivative_fixture": operator.subs({h: 1, hd: 0, f0: 1, f1: 0, f2: 0})}
