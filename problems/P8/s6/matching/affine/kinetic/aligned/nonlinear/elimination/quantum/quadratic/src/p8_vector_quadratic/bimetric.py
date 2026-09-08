"""Quadratic four-derivative pole from the independent scalar determinant."""
from functools import cache

import sympy as sp

from . import geometry, tensors

JETS = ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2))
temporal = {jet: sp.Symbol("temporal_mass_"+"".join(map(str, jet)), real=True) for jet in JETS}
spatial = {jet: sp.Symbol("spatial_mass_"+"".join(map(str, jet)), real=True) for jet in JETS}
amplitude = sp.Symbol("mass_expansion_amplitude", real=True)


def log_jet(field, jet):
    first = amplitude*field[jet]
    if sum(jet) == 1:
        second = -amplitude**2*field[(0, 0)]*field[jet]
    else:
        axes = [(1, 0)]*jet[0]+[(0, 1)]*jet[1]
        second = -amplitude**2*(field[axes[0]]*field[axes[1]]+field[(0, 0)]*field[jet])
    return first+second


def binomial_second(variable, power):
    return 1+power*amplitude*variable+power*(power-1)*amplitude**2*variable**2/2


@cache
def local_quadratic_density():
    # The auxiliary metric is sqrt(det(1+Y/m^2)) (1+Y/m^2)^(-1) g.
    # The pointwise original scale factor can be set to one after coordinate
    # curvature evaluation; all spatial derivative jets below are rescaled
    # by that same pointwise scale factor.
    substitutions = {}
    for jet, lapse_jet, scale_jet in zip(JETS[1:],
                                        (geometry.n10, geometry.n01, geometry.n20, geometry.n11, geometry.n02),
                                        (geometry.v10, geometry.v01, geometry.v20, geometry.v11, geometry.v02)):
        substitutions[lapse_jet] = (-log_jet(temporal, jet)+3*log_jet(spatial, jet))/4
        substitutions[scale_jet] = (log_jet(temporal, jet)+log_jet(spatial, jet))/4
    output = sp.Integer(0)
    for term in sp.Add.make_args(geometry.heat_kernel()["scalar_heat"]):
        powers = term.as_powers_dict()
        p, q = sp.sympify(powers.get(geometry.g0, 0)), sp.sympify(powers.get(geometry.gx, 0))
        coefficient = term/(geometry.g0**p*geometry.gx**q)
        temporal_power = sp.Rational(1, 2)-p/2+q/2
        spatial_power = sp.Rational(3, 2)+3*p/2+q/2
        density_factor = binomial_second(temporal[(0, 0)], temporal_power)*binomial_second(spatial[(0, 0)], spatial_power)
        output += sp.expand(coefficient.subs(substitutions)*density_factor).coeff(amplitude, 2)
    # Gamma_G pole = -a4_scalar/(32*pi^2*epsilon) with epsilon=(4-d)/2.
    return -sp.expand(output)


def fourier_mapping(mode):
    out = {}
    for fields, data in ((temporal, mode.temporal), (spatial, mode.spatial)):
        for (time_order, spatial_order), variable in fields.items():
            out[variable] = (mode.sign*sp.I*tensors.k)**spatial_order*data[time_order]/tensors.mass**2
    return out


@cache
def fourier_bilinear():
    density = local_quadratic_density()
    plus, minus = fourier_mapping(tensors.plus), fourier_mapping(tensors.minus)
    summed = {variable: plus[variable]+minus[variable] for variable in plus}
    return sp.expand((density.subs(summed)-density.subs(plus)-density.subs(minus))/2)
