"""Root-authored audit from physical TT action and independent chain rules.

No production residual dictionary is used as the scientific oracle. Exact
identities here are supplemented by the written uniform-estimate review;
checking a matrix at rational fixtures alone does not certify an interval.
"""

from fractions import Fraction
from functools import cache

import pytest
import sympy as sp
from p8_variable_response import (
    bounds,
    connection,
    exact,
    independent,
    operator,
    source,
)


@cache
def physical_derivation():
    u, c, k = sp.symbols("u c k", real=True, positive=True)
    a, b = (1+u**2)**2, 2/(1+u**2)**2
    y = b/a
    hp = 4*(1-u**2)/(1+u**2)**2
    beta1 = y**3*hp/(c*(c-y))
    spring = a**3*y*beta1/4
    kg, kf = a**3/8, b**3/(8*c)
    gg, gf = a/8, c*b/8
    physical = sp.Matrix([
        [0, 1, 0, 0],
        [-(gg*k**2+spring)/kg, -sp.diff(kg, u)/kg, spring/kg, 0],
        [0, 0, 0, 1],
        [spring/kf, 0, -(gf*k**2+spring)/kf, -sp.diff(kf, u)/kf],
    ])
    fs = sp.sqrt(2*(kg+kf))
    fr = sp.sqrt(2*kg*kf/(kg+kf))
    ag, bg = 1/fs, -(kf/(kg+kf))/fr
    af, bf = ag, (kg/(kg+kf))/fr
    w = sp.Matrix([
        [ag, 0, bg, 0],
        [sp.diff(ag, u), ag, sp.diff(bg, u), bg],
        [af, 0, bf, 0],
        [sp.diff(af, u), af, sp.diff(bf, u), bf],
    ])
    return u, c, k, physical, w


def test_canonical_system_is_the_literal_physical_action_system():
    u, c, k, physical, w = physical_derivation()
    d = operator.derive()
    canonical = operator.first_order().subs({d["u"]: u, d["c"]: c, d["kbar"]: k})
    residual = physical*w-sp.diff(w, u)-w*canonical
    # The independently factored positive square roots need radical
    # simplification, not just rational polynomial factorization.
    assert all(sp.simplify(sp.factor(value)) == 0 for value in residual)


def test_physical_dictionary_includes_every_coefficient_derivative():
    u, c, k, _, w = physical_derivation()
    d = operator.derive()
    actual = operator.physical_map().subs({d["u"]: u, d["c"]: c, d["kbar"]: k})
    assert all(sp.simplify(value) == 0 for value in w-actual)
    assert sp.simplify(w.det()-1/((2*((1+u**2)**6/8))*(2/(c*(1+u**2)**6)))) == 0
    # At a punctured slice a field-only rescaling is not a Cauchy map.
    assert w[1, 0].subs({u: sp.Rational(1, 100), c: 3}) != 0
    assert w[1, 2].subs({u: sp.Rational(1, 100), c: 3}) != 0


@pytest.mark.parametrize("sgn", (-1, 1))
def test_log_radius_generator_by_differentiating_the_state_definitions(sgn):
    r = sp.Symbol("r", positive=True)
    l, lp, q, qp = sp.symbols("l lp q qp", real=True)
    aa, cc, dd, ee, ff, bb = sp.symbols("A C d E f B", real=True)
    mu = sp.sqrt(39)/2
    state = sp.Matrix([l, lp, q/sp.sqrt(r), (sgn*r*qp-q/2)/(mu*sp.sqrt(r))])
    vec = sp.Matrix([lp, -aa*l-cc*q-sgn*r*dd*qp,
                     qp, -ee*l-sgn*r*ff*lp-bb*q])
    dv = sgn*r*state.jacobian([l, lp, q, qp])*vec+r*state.diff(r)
    inv = operator.outer_map(r, sgn)
    expected = sp.Matrix([
        [0, sgn*r, 0, 0],
        [-sgn*r*aa, 0, -sgn*r**sp.Rational(3, 2)*(cc+dd/2),
         -sgn*r**sp.Rational(3, 2)*dd*mu],
        [0, 0, 0, mu],
        [-r**sp.Rational(3, 2)*ee/mu, -sgn*r**sp.Rational(5, 2)*ff/mu,
         -(r*r*bb-sp.Rational(1, 4))/mu, 0],
    ])
    assert all(sp.simplify(v) == 0 for v in dv.jacobian([l, lp, q, qp])*inv-expected)


def test_inner_poschl_teller_equation_without_production_transform():
    z = sp.Symbol("z", real=True)
    psi = sp.Function("psi")(z)
    x = sp.sinh(z)/sp.sqrt(8)
    q = sp.sqrt(sp.cosh(z))*psi
    dx = sp.diff(x, z)
    physical = sp.diff(sp.diff(q, z)/dx, z)/dx+80*q/(1+8*x*x)
    normalized = sp.simplify(physical*dx**2/sp.sqrt(sp.cosh(z)))
    expected = sp.diff(psi, z, 2)+(sp.Rational(39, 4)+3/(4*sp.cosh(z)**2))*psi
    assert sp.simplify(normalized-expected) == 0


def test_jost_hypergeometric_parameters_follow_from_ode():
    y = sp.Symbol("y", positive=True)
    mu = sp.sqrt(39)/2
    f0, f1, f2 = sp.symbols("F F_y F_yy")
    yp = -2*y*(1-y)
    ypp = sp.diff(yp, y)*yp
    transformed = sp.expand((yp**2*f2+(ypp+2*sp.I*mu*yp)*f1+3*y*(1-y)*f0)/(4*y*(1-y)))
    a, b, c = -sp.Rational(1, 2), sp.Rational(3, 2), 1-sp.I*mu
    hypergeometric = y*(1-y)*f2+(c-(a+b+1)*y)*f1-a*b*f0
    assert sp.factor(transformed-hypergeometric) == 0
    assert (c-a-b).is_integer is False
    assert c.is_integer is False


def test_gamma_reflection_connection_and_norm_identity():
    # DLMF 15.10.21, 5.5.1 and 5.5.3: for real mu>0,
    # |Gamma(1-i mu)Gamma(-i mu)|^2=pi^2/sinh(pi mu)^2.
    # Recurrence cancels the |1/2-i mu| factors in the denominator,
    # leaving pi^2/cosh(pi mu)^2. Hence |A|^2=coth(pi mu)^2.
    mu = sp.Symbol("mu", positive=True)
    aa_squared = sp.cosh(sp.pi*mu)**2/sp.sinh(sp.pi*mu)**2
    bb = sp.I/sp.sinh(sp.pi*mu)
    assert sp.simplify(aa_squared-bb*sp.conjugate(bb)-1) == 0
    assert sp.simplify((sp.pi/sp.sin(sp.pi*sp.I*mu))/(-sp.pi)-bb) == 0
    assert sp.gamma(-sp.Rational(1, 2))*sp.gamma(sp.Rational(3, 2)) == -sp.pi


def test_pole_subtracted_mass_has_path_dependent_bounded_limits():
    u, delta = sp.symbols("u delta", real=True, positive=True)
    c, den = 2+delta, 1+u**2
    mass = 16*(1-u**2)*(8+c*den**12)/(c*den**14*(c*den**4-2))
    remainder = mass-80/(delta+8*u**2)
    center_first = sp.limit(remainder.subs(u, 0), delta, 0, dir="+")
    outer_first = sp.limit(sp.cancel(remainder.subs(delta, 0)), u, 0, dir="+")
    assert center_first == -32
    assert outer_first == -141
    assert center_first != outer_first  # bounded does not mean jointly analytic
    d = operator.derive()
    analytic = d["b_analytic"]
    assert sp.factor(analytic-d["V_HH"]+d["mass"]) == 0
    assert sp.limit(analytic.subs(d["u"], 0), d["c"], 2) == d["kbar"]**2-sp.Rational(18, 5)


@pytest.mark.parametrize("bad", [0, -1, True, 0.000001, sp.Rational(1, 10**8)])
def test_finite_parent_domain_rejects_singular_or_outside_delta(bad):
    with pytest.raises((TypeError, ValueError)):
        exact.parameters(bad)


@pytest.mark.parametrize("bad", [0, 3, True, 1.0, sp.oo, sp.nan])
def test_fixed_momentum_domain_cannot_be_silently_broadened(bad):
    with pytest.raises((TypeError, ValueError)):
        exact.parameters(sp.Rational(1, 10**12), bad)


def test_convergence_bound_is_not_small_on_its_entire_parameter_box():
    assert 40000*sp.real_root(sp.Rational(1, 10**9), 3) == 40
    assert 40000*sp.real_root(sp.Rational(1, 10**21), 3) == sp.Rational(1, 250)


def test_positive_frequency_labels_reverse_across_the_bounce():
    # For exp(-i integral omega du), omega=mu/|u|:
    # u<0 has +i mu log(-u), u>0 has -i mu log(u).
    r, mu = sp.symbols("r mu", positive=True)
    incoming = sp.I*mu*sp.log(r)
    outgoing = -sp.I*mu*sp.log(r)
    assert -sp.diff(incoming, r) == -sp.I*mu/r
    assert sp.diff(outgoing, r) == -sp.I*mu/r


def test_g_metric_jet_observes_all_four_homogeneous_physical_components():
    u, c, k, physical, _ = physical_derivation()
    first_two = sp.eye(4)[:2, :]
    g_second = physical[1, :]
    g_third = g_second.diff(u)+g_second*physical
    jet_map = first_two.col_join(g_second).col_join(g_third)
    spring_over_kg = physical[1, 2]
    assert sp.factor(jet_map.det()-spring_over_kg**2) == 0
    for sign in (-1, 1):
        value = spring_over_kg.subs({u: sp.Rational(sign, 100), c: 2, k: 1})
        assert value.is_positive is True  # punctured limit, not the singular center


def test_actual_g_probe_projection_from_literal_action():
    d = operator.derive()
    u, c, k, _, w = physical_derivation()
    coeff = source.coefficients()
    projection = ((1+u**2)**6/2)*w[0, :]
    substitutions = {d["u"]: u, d["c"]: c, d["kbar"]: k}
    assert sp.simplify(projection[0]-coeff["j_light"].subs(substitutions)) == 0
    assert sp.simplify(projection[2]-coeff["j_heavy"].subs(substitutions)) == 0
    assert projection[1] == projection[3] == 0


def test_g_only_source_loading_is_a_variational_identity():
    # Independent generic physical action, with neither an f source nor a
    # special solution ansatz in its variational derivative.
    u = sp.Symbol("u", real=True)
    kg, kf, gg, gf, spring, a = [sp.Function(name)(u) for name in
                                ("Kg", "Kf", "Gg", "Gf", "U", "a")]
    k = sp.Symbol("k", positive=True)
    g, f = sp.Function("g")(u), sp.Function("f")(u)
    sigma = sp.Function("sigma")(u)
    lagrangian = (kg*sp.diff(g, u)**2+kf*sp.diff(f, u)**2
                  -k**2*(gg*g**2+gf*f**2)-spring*(g-f)**2+a**3*sigma*g/2)

    def euler(field):
        return sp.diff(sp.diff(lagrangian, sp.diff(field, u)), u)-sp.diff(lagrangian, field)

    reconstructed_g = f+(sp.diff(kf*sp.diff(f, u), u)+gf*k**2*f)/spring
    reconstructed_sigma = 4/a**3*(sp.diff(kg*sp.diff(reconstructed_g, u), u)
                                       +gg*k**2*reconstructed_g+spring*(reconstructed_g-f))
    assert sp.simplify(euler(f).subs(g, reconstructed_g).doit()) == 0
    g_euler = euler(g).subs(g, reconstructed_g).doit().subs(sigma, reconstructed_sigma)
    assert sp.expand(g_euler) == 0


def test_complex_power_basis_is_norm_preserving_not_a_cauchy_symplectic_claim():
    unitary = sp.Matrix([[1, -sp.I], [1, sp.I]])/sp.sqrt(2)
    assert unitary == connection.unitary_power_coordinates()
    assert unitary.conjugate().T*unitary == sp.eye(2)
    aa, ac, bb, bc, phase = sp.symbols("A Ac B Bc phase", nonzero=True)
    left_from_right = sp.Matrix([[bb, ac*phase], [aa/phase, bc]])
    candidate = sp.Matrix([[-bc, ac*phase], [aa/phase, -bb]])
    identity_factor = aa*ac-bb*bc
    assert sp.expand(candidate.det()+identity_factor) == 0
    assert all(sp.expand(v) == 0 for v in candidate*left_from_right-identity_factor*sp.eye(2))


def test_finite_prepared_light_bound_retains_source_loading_error():
    delta = sp.Rational(1, 10**21)
    target, source_l1 = sp.Integer(2), sp.Integer(3)
    transfer = sp.Rational(1, 250)
    loading = sp.Rational(9, 10**15)
    assert bounds.prepared_error(delta, target, source_l1) == 2*transfer+(2+transfer)*loading
    assert bounds.prepared_error(delta, target, source_l1) > 2*transfer


def test_two_log_phase_subsequences_have_different_limits():
    mu, scale, integer = sp.symbols("mu scale n", positive=True)
    # Both eps sequences tend to zero. The associated transmission phases
    # are exactly +1 and -1; B need not be nonzero for this distinction.
    log_eps_one = sp.log(scale)-sp.pi*integer/mu
    log_eps_two = sp.log(scale)-sp.pi*(integer+sp.Rational(1, 2))/mu
    exponent_one = sp.expand(2*sp.I*mu*(log_eps_one-sp.log(scale)))
    exponent_two = sp.expand(2*sp.I*mu*(log_eps_two-sp.log(scale)))
    assert sp.simplify(sp.exp(exponent_two-exponent_one)) == -1


@cache
def differentiated_analytic_coefficients():
    d = operator.derive()
    keys = ("A", "C", "d_cross", "f_cross", "E", "b_analytic")
    return d, {key: sp.factor(sp.diff(d[key], d["c"])) for key in keys}


@pytest.mark.parametrize("u0", [Fraction(-1, 10), Fraction(0), Fraction(1, 100)])
@pytest.mark.parametrize("c0", [Fraction(2), Fraction(4004001, 2000000)])
def test_fraction_mixed_jets_contain_independently_differentiated_coefficients(u0, c0):
    interval = independent.Interval
    values = independent.profiles(interval(u0), interval(c0), interval(2))
    d, differentiated = differentiated_analytic_coefficients()
    substitutions = {d["u"]: sp.Rational(u0), d["c"]: sp.Rational(c0), d["kbar"]: 2}
    for key, expression in differentiated.items():
        actual = sp.simplify(expression.subs(substitutions))
        lower, upper = values[key].get((0, 1)).pair()
        assert (actual-sp.Rational(lower)).is_nonnegative is True
        assert (sp.Rational(upper)-actual).is_nonnegative is True


@pytest.mark.parametrize("value", [Fraction(0), Fraction(2), Fraction(201, 100), Fraction(1, 10**18)])
def test_fraction_square_root_enclosure_is_outward_by_integer_squaring(value):
    enclosure = independent.Interval(value).sqrt()
    assert enclosure.lo >= 0
    assert enclosure.lo**2 <= value <= enclosure.hi**2


def test_closed_analytic_box_is_continuously_enclosed_and_excludes_the_pole():
    result = independent.analytic_boxes()
    assert result["closed_u_box"] == (Fraction(-1, 10), Fraction(1, 10))
    assert result["closed_c_box"] == (Fraction(2), Fraction(201, 100))
    assert result["u_subintervals"] == 40
    for key, cap in result["analytic_c_derivative_abs_upper"].items():
        assert cap < (1 if key.startswith("j_") else 8)
    assert "B" not in result["analytic_c_derivative_abs_upper"]
    assert "mass" not in result["analytic_c_derivative_abs_upper"]
