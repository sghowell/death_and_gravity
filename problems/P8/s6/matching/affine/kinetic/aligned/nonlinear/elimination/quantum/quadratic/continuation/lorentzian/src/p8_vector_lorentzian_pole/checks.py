"""Physical-signature, frozen curvature and Lorentzian Feynman checks."""
from functools import cache
from itertools import combinations

import sympy as sp
from p8_affine_kinetic import scalar
from p8_vector_bubble import pole as bubble
from p8_vector_curvature import geometry as frozen_geometry
from p8_vector_quadratic import tensors
from p8_vector_quadratic_matching import geometry, kernel

from . import bimetric, frame, response, wick


@cache
def frame_checks():
    out = {"Lorentzian_Ricci_contraction_"+str(i)+str(j): sp.expand(
        sum(frame.eta(r)*frame.riemann(r, i, r, j) for r in range(4))-frame.ricci(i, j).subs(geometry.dimension, 3))
        for i in range(4) for j in range(4)}
    for a, b in combinations(range(4), 2):
        out["Lorentzian_mass_commutator_"+str(a)+str(b)] = sp.ImmutableMatrix([[sp.expand(
            frame.covariant(tensors.plus, (a, b), i, j)-frame.covariant(tensors.plus, (b, a), i, j)
            -sum(frame.eta(r)*(frame.riemann(a, b, i, r)*frame.covariant(tensors.plus, (), r, j)
                              +frame.riemann(a, b, j, r)*frame.covariant(tensors.plus, (), i, r)) for r in range(4)))
            for j in range(4)] for i in range(4)])
    return out


@cache
def signature_checks():
    return {"direct_Lorentzian_second_pole_equals_local_continuation_in_all_D":
            sp.expand(frame.second_loop_pole()+wick.continuation(kernel.pole(2))),
            "direct_Lorentzian_fourth_pole_equals_local_continuation_in_all_D":
            sp.expand(bimetric.loop_pole()+wick.continuation(kernel.pole(4)))}


@cache
def physical_geometry():
    old, g = frozen_geometry.rolling(), geometry.four
    values = {field: 0 for field in (g.n10, g.n01, g.n20, g.n11, g.n02,
                                     g.v10, g.v01, g.v20, g.v11, g.v02)}
    values.update({g.g0: -1, g.gx: 1, geometry.dimension: 3})
    values.update({field: sp.diff(scalar.background()["H"], scalar.u, j) for j, field in enumerate(tensors.H)})
    names = {"scalar": "R", "ricci_squared": "Ricci2", "riemann_squared": "Riemann2"}
    return {"physical_clock_"+name+"_replays_S6_49": sp.factor(
        geometry.heat_kernel()[name].subs(values)-old[prior_name])
        for name, prior_name in names.items()}


@cache
def flat():
    alpha, beta, frequency = sp.symbols("alpha beta physical_time_frequency", real=True)
    m, k = tensors.mass, tensors.k
    mapping = {value: 0 for value in tensors.H}
    for mode in (tensors.plus, tensors.minus):
        for order in range(5):
            mapping[mode.temporal[order]] = m**2*alpha*(mode.sign*sp.I*frequency)**order
            mapping[mode.spatial[order]] = m**2*beta*(mode.sign*sp.I*frequency)**order
    old_mapping = {bubble.m2: m**2, bubble.P2: k**2-frequency**2,
                   bubble.TA: alpha+3*beta, bubble.TA2: alpha**2+3*beta**2,
                   bubble.PA: beta*k**2-alpha*frequency**2,
                   bubble.PA2: beta**2*k**2-alpha**2*frequency**2}
    names = {0: "constant", 2: "second", 4: "fourth"}
    actual = {0: frame.zero_loop_control(),
              2: frame.second_loop_pole().subs(geometry.dimension, 3),
              4: bimetric.loop_pole().subs(geometry.dimension, 3)}
    return {"Lorentzian_flat_Feynman_pole_"+str(order): sp.factor(
        actual[order].subs(mapping)-bubble.parameter_integral()[name].subs(old_mapping)/2)
        for order, name in names.items()}


@cache
def history_control():
    actual = scalar.background()["a"]
    changed = sp.factor(actual.subs(scalar.u, sp.I*scalar.u))
    fixture = {scalar.u: sp.Rational(1, 2)}
    return {"original_physical_scale_factor": actual, "naive_time_rotated_scale_factor": changed,
            "nonzero_history_difference": sp.factor(changed-actual),
            "original_half_time_scale_factor": actual.subs(fixture),
            "changed_half_time_scale_factor": changed.subs(fixture),
            "fourth_pole_second_time_square_bounce": response.compact(4, 0)["second_time_derivative_squared"].subs(scalar.u, 0)}
