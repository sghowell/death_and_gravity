"""Independent four-index, variational and cosmic-time Weyl-tensor audits."""

import itertools

import sympy as sp
from p8_m1_weyl import reduction


def test_weyl_tensor_four_index_dynamic_contraction():
    tt, tz, zz = sp.symbols("gamma_tt gamma_tz gamma_zz")
    signs = (1, -1, -1, -1)

    def hessian(i, j, a, b):
        polarization = 1 if (i, j) == (1, 1) else -1 if (i, j) == (2, 2) else 0
        return polarization*{(0, 0): tt, (0, 3): tz, (3, 0): tz, (3, 3): zz}.get((a, b), 0)

    curvature = {}
    for a, b, c, d in itertools.product(range(4), repeat=4):
        curvature[a, b, c, d] = (
            hessian(a, d, c, b)+hessian(b, c, d, a)
            - hessian(a, c, d, b)-hessian(b, d, c, a))/sp.Integer(2)
    ricci = {(b, d): sum(signs[a]*curvature[a, b, a, d] for a in range(4))
             for b, d in itertools.product(range(4), repeat=2)}
    scalar = sum(signs[a]*ricci[a, a] for a in range(4))
    riemann_squared = sum(sp.prod(signs[i] for i in indices)*value**2
                          for indices, value in curvature.items())
    ricci_squared = sum(signs[b]*signs[d]*value**2 for (b, d), value in ricci.items())
    # The literal polarization diag(1,-1,0) has norm^2=2. Divide by 2
    # to compare with the declared unit-normalized physical tensor.
    weyl_squared = sp.expand((riemann_squared-2*ricci_squared+scalar**2/3)/2)
    assert scalar == 0
    assert sp.expand(weyl_squared-((tt+zz)**2/2-2*tz**2)) == 0
    # The difference from (Box gamma)^2/2 is a Hessian-determinant
    # divergence, not a pointwise zero and not a license to drop C^2.
    assert sp.expand(weyl_squared-(tt-zz)**2/2) == 2*(tt*zz-tz**2)
    assert weyl_squared.subs({tt: 0, tz: 0, zz: 1}) == sp.Rational(1, 2)
    assert weyl_squared.subs({tt: 1, tz: -1, zz: 1}) == 0


def test_weyl_spatial_average_and_total_time_derivative():
    t, k = sp.symbols("t k", real=True)
    mode = sp.Function("gamma")(t)
    # For gamma(t,z)=mode(t)*cos(kz), averaging the literal pointwise
    # four-index contraction gives the following density.
    literal_average = (sp.diff(mode, t, 2)-k**2*mode)**2/4-k**2*sp.diff(mode, t)**2
    reduced_average = (sp.diff(mode, t, 2)+k**2*mode)**2/4
    divergence = -k**2*sp.diff(mode*sp.diff(mode, t), t)
    assert sp.expand(literal_average-reduced_average-divergence) == 0


def test_weyl_action_sign_and_full_off_shell_field_map():
    t, k, mass, coefficient = sp.symbols("t k M c_C", real=True)
    a = sp.Function("a")(t)
    mode = sp.Function("gamma")(t)
    hubble = sp.diff(a, t)/a
    momentum_squared = k**2/a**2
    e0 = sp.diff(mode, t, 2)+3*hubble*sp.diff(mode, t)+momentum_squared*mode
    density0 = mass**2*a**3*(sp.diff(mode, t)**2-momentum_squared*mode**2)/8
    variation0 = sp.diff(density0, mode)-sp.diff(sp.diff(density0, sp.diff(mode, t)), t)
    assert sp.simplify(variation0+mass**2*a**3*e0/4) == 0
    density1 = coefficient*a**3*(e0-2*hubble*sp.diff(mode, t))**2/2
    full_map = (2*e0-8*hubble*sp.diff(mode, t))/mass**2
    mapped_first_order = density1+coefficient*variation0*full_map
    expected = 2*coefficient*a**3*hubble**2*sp.diff(mode, t)**2
    assert sp.simplify(mapped_first_order-expected) == 0
    # Omitting 2 beta E0 in an off-shell action identity leaves E0^2.
    incomplete_map = -8*hubble*sp.diff(mode, t)/mass**2
    wrong_residual = density1+coefficient*variation0*incomplete_map-expected
    assert sp.simplify(wrong_residual-coefficient*a**3*e0**2/2) == 0


def test_weyl_operator_factorization_by_direct_cosmic_differentiation():
    t, k = sp.symbols("t k", real=True)
    a = sp.Function("a")(t)
    mode = sp.Function("gamma")(t)
    hubble = sp.diff(a, t)/a
    hubble_dot = sp.diff(hubble, t)
    momentum_squared = k**2/a**2
    b = sp.diff(hubble, t, 2)-2*hubble*hubble_dot
    l0 = lambda value: sp.diff(value, t, 2)+3*hubble*sp.diff(value, t)+momentum_squared*value
    conformal_d = lambda value: a*sp.diff(a*sp.diff(value, t), t)+k**2*value
    direct = conformal_d(conformal_d(mode))/a**4
    factorized = l0(l0(mode))+2*(hubble**2-hubble_dot)*l0(mode)
    factorized += 4*hubble_dot*momentum_squared*mode-2*b*sp.diff(mode, t)
    assert sp.simplify(direct-factorized) == 0
    variables = {reduction.H: hubble, reduction.H1: hubble_dot,
                 reduction.H2: sp.diff(hubble, t, 2), reduction.Q: momentum_squared}
    exported = l0(l0(mode))+reduction.J.subs(variables)*l0(mode)
    exported -= (reduction.F.subs(variables)*sp.diff(mode, t)
                 + reduction.C.subs(variables)*mode)/4
    assert sp.simplify(direct-exported) == 0
