"""Independent local heat coefficient using the physical negative lapse metric."""
from functools import cache

import sympy as sp
from p8_vector_quadratic import bimetric as four
from p8_vector_quadratic import geometry as coordinate
from p8_vector_quadratic import tensors
from p8_vector_quadratic_matching import geometry


@cache
def local_loop_density():
    d = geometry.dimension
    substitutions = {}
    for jet, lapse_jet, scale_jet in zip(four.JETS[1:],
                                        (coordinate.n10, coordinate.n01, coordinate.n20, coordinate.n11, coordinate.n02),
                                        (coordinate.v10, coordinate.v01, coordinate.v20, coordinate.v11, coordinate.v02)):
        temporal, spatial = four.log_jet(four.temporal, jet), four.log_jet(four.spatial, jet)
        substitutions[lapse_jet] = (-temporal+d*spatial)/4
        substitutions[scale_jet] = (temporal+(d-2)*spatial)/4
    output = sp.Integer(0)
    for term in sp.Add.make_args(geometry.heat_kernel()["scalar_heat"]):
        powers = term.as_powers_dict()
        p = sp.sympify(powers.get(coordinate.g0, 0))
        q = sp.sympify(powers.get(coordinate.gx, 0))
        coefficient = term/(coordinate.g0**p*coordinate.gx**q)
        temporal_power = (d-1)/4-p/2+q/2
        spatial_power = d*(d-1)/4+d*p/2+(d-2)*q/2
        density_factor = four.binomial_second(four.temporal[(0, 0)], temporal_power)*four.binomial_second(four.spatial[(0, 0)], spatial_power)
        # g00=-N^2, while the relative absolute volume is unchanged.
        # The Lorentzian loop residue has +a4_scalar.
        output += (-1)**p*sp.expand(coefficient.subs(substitutions)*density_factor).coeff(four.amplitude, 2)
    return sp.expand(output)


@cache
def loop_pole():
    fields = tuple(four.temporal.values())+tuple(four.spatial.values())
    plus, minus = four.fourier_mapping(tensors.plus), four.fourier_mapping(tensors.minus)
    output = sp.Integer(0)
    for powers, coefficient in sp.Poly(local_loop_density(), *fields).terms():
        factors = [field for field, power in zip(fields, powers) for _ in range(power)]
        if len(factors) != 2:
            raise ValueError("Require a homogeneous quadratic mass density")
        left, right = factors
        output += coefficient*(plus[left]*minus[right]+minus[left]*plus[right])/2
    return sp.expand(output)
