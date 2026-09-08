"""Independent flat Feynman-kernel and conformal-mass pole checks."""
from functools import cache
from itertools import pairwise

import sympy as sp
from p8_vector_bubble import pole

from . import bimetric, geometry, invariants, kernel, tensors


@cache
def flat():
    alpha, beta, time_frequency = sp.symbols("alpha beta Euclidean_time_frequency", real=True)
    k, m = tensors.k, tensors.mass
    mapping = {value: 0 for value in tensors.H}
    for mode in (tensors.plus, tensors.minus):
        for order in range(5):
            mapping[mode.temporal[order]] = m**2*alpha*(mode.sign*sp.I*time_frequency)**order
            mapping[mode.spatial[order]] = m**2*beta*(mode.sign*sp.I*time_frequency)**order
    frozen = pole.parameter_integral()
    old_mapping = {pole.m2: m**2, pole.P2: time_frequency**2+k**2,
                   pole.TA: alpha+3*beta, pole.TA2: alpha**2+3*beta**2,
                   pole.PA: alpha*time_frequency**2+beta*k**2,
                   pole.PA2: alpha**2*time_frequency**2+beta**2*k**2}
    return {"flat_frozen_Feynman_pole_order_"+str(order): sp.factor(
        kernel.pole(order).subs(mapping)+frozen[name].subs(old_mapping)/2)
        for order, name in ((0, "constant"), (2, "second"), (4, "fourth"))}


def total_derivative(value, field_jets):
    flows = {**dict(pairwise(tensors.H)), tensors.k: -tensors.H[0]*tensors.k,
             **dict(pairwise(field_jets))}
    return sp.expand(sum(sp.diff(value, var)*flows[var] for var in value.free_symbols if var in flows))


def euler_lagrange(value, field_jets):
    out = sp.Integer(0)
    for order in range(5):
        term = sp.diff(value, field_jets[order])
        for _ in range(order):
            term = -total_derivative(term, field_jets)-3*tensors.H[0]*term
        out += term
    return sp.factor(out)


@cache
def scalar_conformal():
    n = sp.symbols("scalar_mass_profile0:9", real=True)
    m, H, k = tensors.mass, tensors.H, tensors.k
    mapping = {mode.temporal[order]: m**2*n[order] for mode in (tensors.plus, tensors.minus) for order in range(5)}
    mapping.update({mode.spatial[order]: m**2*n[order] for mode in (tensors.plus, tensors.minus) for order in range(5)})
    box = n[2]+3*H[0]*n[1]-k**2*n[0]
    gradient_squared = n[1]**2+k**2*n[0]**2
    R = tensors.scalar_curvature()
    expected = {0: -sp.Rational(3, 2)*m**4*n[0]**2,
                2: -sp.Rational(3, 4)*m**2*gradient_squared,
                4: -(9*box**2+6*R*n[0]*box+3*R*gradient_squared)/72}
    raw = {order: sp.expand(kernel.pole(order).subs(mapping)) for order in (0, 2, 4)}
    return {"field_jets": n, "raw_scalar_mass_pole": raw, "independent_conformal_Proca_pole": expected,
            "compact_variation_checks": {"scalar_mass_conformal_action_order_"+str(order): euler_lagrange(raw[order]-expected[order], n)
                                          for order in (0, 2, 4)}}

@cache
def general_second():
    t = sp.symbols("independent_temporal0:9", real=True)
    s = sp.symbols("independent_spatial0:9", real=True)
    mapping = {}
    for mode in (tensors.plus, tensors.minus):
        mapping.update({value: tensors.mass**2*t[order] for order, value in enumerate(mode.temporal)})
        mapping.update({value: tensors.mass**2*s[order] for order, value in enumerate(mode.spatial)})
    difference = sp.expand((kernel.pole(2)-kernel.independent_second()).subs(mapping))
    # Both independent profiles flow under every integration by parts.
    all_jets = (t, s)
    def time(value):
        flows = {**dict(pairwise(tensors.H)), tensors.k: -tensors.H[0]*tensors.k}
        for jets in all_jets:
            flows.update(dict(pairwise(jets)))
        return sp.expand(sum(sp.diff(value, var)*flow for var, flow in flows.items()))
    out = {}
    for name, jets in (("temporal", t), ("spatial", s)):
        value = 0
        for order in range(3):
            term = sp.diff(difference, jets[order])
            for _ in range(order):
                term = -time(term)-3*tensors.H[0]*term
            value += term
        out["independent_curved_second_derivative_"+name] = sp.factor(value)
    return out


@cache
def constant_anisotropic():
    t, s = bimetric.temporal[(0, 0)], bimetric.spatial[(0, 0)]
    mapping = {variable: 0 for group in (bimetric.temporal, bimetric.spatial)
               for jet, variable in group.items() if jet != (0, 0)}
    background = geometry.heat_kernel()["scalar_heat"].subs({
        value: 0 for value in (geometry.n10, geometry.n01, geometry.n20, geometry.n11, geometry.n02,
                               geometry.v10, geometry.v01, geometry.v20, geometry.v11, geometry.v02)})
    background = background.subs({geometry.g0: 1, geometry.gx: 1})
    expected = -background*(3*t**2-18*t*s+15*s**2)/8
    return {"constant_anisotropic_auxiliary_metric_scaling":
            sp.expand(bimetric.local_quadratic_density().subs(mapping)-expected)}


@cache
def withheld_control():
    data = scalar_conformal()
    n = data["field_jets"]
    mapping = {value: tensors.mass**2*n[order] for mode in (tensors.plus, tensors.minus)
               for group in (mode.temporal, mode.spatial) for order, value in enumerate(group)}
    literal = sp.expand(invariants.evaluate(4).subs(mapping)/(960*tensors.mass**4))
    residual = euler_lagrange(literal-data["independent_conformal_Proca_pole"][4], n)
    fixture = {value: 0 for value in tensors.H+n}
    fixture.update({tensors.H[0]: 1, n[1]: 1, tensors.k: 0})
    return {"literal_transcribed_fourth_order_scalar_pole": literal,
            "nonzero_compact_variation_residual": residual,
            "nonzero_exact_fixture": residual.subs(fixture),
            "status": "WITHHELD_LITERAL_REFERENCE_NOT_ACCEPTED_POLE"}
