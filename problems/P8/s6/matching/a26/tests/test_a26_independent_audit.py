"""Separately authored covariant audit; not generated from the report."""

from fractions import Fraction as Q
from functools import cache

import pytest
import sympy as sp


def _literal_coefficients(x, a, g, printed=False):
    f = sp.Rational(1, 2)-g*(1+x)
    ax = -a*(1+x)
    fx = -g
    a3 = 2*(ax-2*fx)*(ax*x-2*f)/(x*(3*ax*x-4*f))
    disputed = 4*(3*f+16*x*fx*ax**2) if printed else 4*(3*f+16*x*fx)*ax**2
    a4 = (-16*x*ax**3+disputed-x*x*f*a3*a3
          -(16*x*x*fx-12*x*f)*a3*ax-16*fx*(3*f+4*x*fx)*ax
          +8*f*(x*fx-f)*a3+48*f*fx*fx)/(8*(f-x*ax)**2)
    a5 = ((4*fx-2*ax+x*a3)*(-2*ax*ax-3*x*ax*a3+4*fx*ax+4*f*a3)
          /(8*(f-x*ax)**2))
    return f, ax, a3, a4, a5


@cache
def _flat_jet_euler():
    # Full symmetric-H differentiation precedes diagonal evaluation.
    v = sp.Matrix(sp.symbols("p q s", real=True))
    hh = sp.symbols("h00 h11 h22 h01 h02 h12", real=True)
    h = sp.Matrix([[hh[0], hh[3], hh[4]], [hh[3], hh[1], hh[5]],
                   [hh[4], hh[5], hh[2]]])
    eta = sp.diag(-1, 1, 1)
    raised = eta*v
    x = (v.T*eta*v)[0]
    lag = ((raised.T*h*raised)[0]*sp.trace(eta*h)
           -(raised.T*h*eta*h*raised)[0])/x
    diagonal = sp.symbols("d0:3", real=True)
    at = dict(zip(hh, (*diagonal, 0, 0, 0), strict=True))
    euler = -sum(diagonal[i]*sp.diff(lag, v[i], 2).subs(at) for i in range(3))
    euler += sum(diagonal[i]**2*sp.diff(lag, hh[i], v[i], 2).subs(at) for i in range(3))
    euler += sum(diagonal[i]*diagonal[j]*sp.diff(lag, hh[k], v[i], v[j]).subs(at)
                 for i, j, k in ((0, 1, 3), (0, 2, 4), (1, 2, 5)))
    return v, diagonal, x, sp.factor(euler)


@pytest.mark.parametrize("printed", [False, True])
def test_full_source_coefficients_have_the_same_leading_residues(printed):
    x, a, g = sp.symbols("X a g", real=True)
    f, _, a3, a4, a5 = _literal_coefficients(x, a, g, printed)
    r = 2*g-a
    for coefficient, expected in ((a3, r), (a4, -r), (a5, -r*r/(2*f.subs(x, 0)))):
        assert sp.factor(sp.cancel(x*coefficient).subs(x, 0)-expected) == 0


def test_printed_ia_parenthesis_difference_is_regular_not_silently_repaired():
    x, a, g = sp.symbols("X a g", real=True)
    f, ax, _, standard, _ = _literal_coefficients(x, a, g)
    printed = _literal_coefficients(x, a, g, printed=True)[3]
    difference = sp.factor(printed-standard)
    assert sp.factor(difference-3*f*(1-ax*ax)/(2*(f-x*ax)**2)) == 0
    assert sp.factor(sp.cancel(x*difference).subs(x, 0)) == 0


@pytest.mark.parametrize("gradient,diagonal", [((1, 0, 0), (1, 1, 1)),
                                               ((2, 1, 0), (2, 3, 5)),
                                               ((3, 1, 1), (-1, 2, 3)),
                                               ((2, 0, 0), (2, 0, 4))])
def test_generic_covariant_euler_including_offdiagonal_variations(gradient, diagonal):
    v, dd, x, euler = _flat_jet_euler()
    assert sp.factor(x*euler+2*sp.prod(dd)) == 0
    at = dict(zip(v, gradient, strict=True)) | dict(zip(dd, diagonal, strict=True))
    expected = -2*Q(diagonal[0]*diagonal[1]*diagonal[2],
                    -gradient[0]**2+gradient[1]**2+gradient[2]**2)
    assert euler.subs(at) == expected


def test_directional_euler_is_not_a_linear_vacuum_operator():
    v, dd, _, euler = _flat_jet_euler()
    outputs = []
    for p, factor in ((1, 1), (2, 1), (3, 2)):
        at = dict(zip(v, (p, 0, 0), strict=True)) | dict.fromkeys(dd, factor)
        outputs.append(euler.subs(at))
    assert outputs == [2, sp.Rational(1, 2), sp.Rational(16, 9)]
    assert outputs[2]-outputs[0]-outputs[1] == -sp.Rational(13, 18)


def test_homogeneous_rank_one_cancellation_prevents_a_false_pole_only_argument():
    v, dd, _, euler = _flat_jet_euler()
    assert euler.subs(dict(zip(v, (1, 0, 0), strict=True)) | dict(zip(dd, (1, 0, 0), strict=True))) == 0


@pytest.mark.parametrize("a,g", [(sp.Rational(1, 3), sp.Rational(1, 5)),
                                 (-2, -sp.Rational(1, 4)),
                                 (sp.Rational(2, 5), sp.Rational(1, 5))])
def test_full_four_dimensional_small_field_density_not_individual_poles(a, g):
    eps = sp.Symbol("eps", positive=True)
    eta = sp.diag(-1, 1, 1, 1)
    v = sp.Matrix([3, 1, 1, 0])
    h = sp.Matrix([[1, sp.Rational(1, 3), 0, 0], [sp.Rational(1, 3), 2, 0, 0],
                   [0, 0, 3, sp.Rational(1, 2)], [0, 0, sp.Rational(1, 2), 4]])
    raised = eta*v
    x = (v.T*eta*v)[0]
    vhv = (raised.T*h*raised)[0]
    box = sp.trace(eta*h)
    l1 = sp.trace(eta*h*eta*h)
    l3, l4, l5 = vhv*box, (raised.T*h*eta*h*raised)[0], vhv**2
    _, ax, a3, a4, a5 = _literal_coefficients(eps**2*x, a, g)
    literal = ax*eps**2*(l1-box**2)+eps**4*(a3*l3+a4*l4)+eps**6*a5*l5
    leading = sp.cancel(literal/eps**2).subs(eps, 0)
    assert sp.factor(leading-(-a*(l1-box**2)+(2*g-a)*(l3-l4)/x)) == 0


def test_zero_residue_is_not_wrongly_assigned_the_generic_quadratic_obstruction():
    x = sp.Symbol("X", real=True)
    _, _, *aa = _literal_coefficients(x, sp.Rational(2, 5), sp.Rational(1, 5))
    assert all(sp.cancel(coefficient).subs(x, 0).is_finite for coefficient in aa)


def test_published_bounce_jet_and_analytic_nonzero_residue_numerator():
    u = sp.Symbol("u", real=True)
    d = 1+u*u
    sigmoid = 1/(1+sp.exp(-u))
    scale = sigmoid*d**sp.Rational(1, 6)+(1-sigmoid)*d**sp.Rational(1, 20)
    h = sp.diff(scale, u)/scale
    hp0 = sp.diff(h, u).subs(u, 0)
    b = u*u*sp.tanh(u+sp.Rational(1, 10))+sp.tanh(u)
    numerator = 9*d*h-2*b
    assert hp0 == sp.Rational(13, 60)
    assert sp.diff(numerator, u).subs(u, 0) == -sp.Rational(1, 20)
    assert sp.diff(numerator, u).subs(u, 0)/(4*hp0) == -sp.Rational(3, 52)


def test_hypothetical_zero_gradient_sign_does_not_replace_rolling_tensor_sign():
    # cosh(x)<=sum (x²/2)^n for |x|<sqrt(2), hence an exact rational bound.
    cosh_upper = 1/(1-Q(1, 200))
    g_lower = Q(2, 3)/cosh_upper**2
    assert g_lower == Q(39601, 60000) > Q(1, 2)
    x, a, g = sp.symbols("X a g", real=True)
    f, ax, *_ = _literal_coefficients(x, a, g)
    assert f.subs(x, -1) == sp.Rational(1, 2)
    assert (f-x*ax).subs(x, -1) == sp.Rational(1, 2)
    assert f.subs(x, 0) == sp.Rational(1, 2)-g


def test_fractional_kinetic_profile_has_c2_not_c3_leading_power():
    alpha = sp.Rational(509, 250)
    assert 2 < alpha < 3
    coefficient = -alpha*(alpha-1)*(alpha-2)/4
    assert coefficient < 0 and alpha-3 == -sp.Rational(241, 250)
    assert 1 < alpha-1 < 2


def test_exact_entropy_field_shift_retains_cancellations_and_does_not_imply_instability():
    q, qp, b, bp, bpp, h, x, box, xi, z, cross = sp.symbols("q qp b bp bpp h X box xi z cross", nonzero=True)
    w = -(q*bpp+(3*h*q+qp)*bp)/b
    original = q*x*bp**2-q*(bp**2*x+2*bp*cross+z)+w*(b**2-(b+xi)**2)
    # Remove only the displayed cross derivative, adding its exact IBP bulk term.
    integrated = sp.expand(original+2*q*bp*cross+2*xi*((qp*bp+q*bpp)*x+q*bp*box))
    expected = -q*z-w*xi**2+2*xi*((qp*bp+q*bpp)*(x+1)+q*bp*(box+3*h))
    assert sp.expand(integrated-expected) == 0
    assert sp.expand(integrated.subs({xi: 0, z: 0})) == 0


@pytest.mark.parametrize("a,g", [(sp.Rational(1, 3), sp.Rational(1, 5)),
                                 (-2, -sp.Rational(1, 4)),
                                 (sp.Rational(2, 5), sp.Rational(1, 5)), (1, 1)])
def test_primary_residue_interface_matches_literal_independent_limits(a, g):
    from p8_a26_vacuum import coefficients

    x = sp.Symbol("X", real=True)
    f, _, *aa = _literal_coefficients(x, a, g, printed=True)
    actual = coefficients.residues(a, g)
    assert actual["F2_at_zero"] == f.subs(x, 0)
    for name, coefficient in zip(("A3", "A4", "A5"), aa, strict=True):
        assert actual[name] == sp.cancel(x*coefficient).subs(x, 0)
    assert actual["generic_obstruction"] == (2*g-a != 0)


def test_primary_euler_does_not_hide_scalar_clock_or_offdiagonal_factors():
    from p8_a26_vacuum import jets

    v, dd, _, euler = _flat_jet_euler()
    for gradient, diagonal in (((1, 0, 0), (1, 1, 1)), ((2, 1, 0), (2, 3, 5)),
                              ((3, 1, 1), (-1, 2, 3)), ((2, 0, 0), (2, 0, 4))):
        at = dict(zip(v, gradient, strict=True)) | dict(zip(dd, diagonal, strict=True))
        assert jets.evaluate_euler(gradient, diagonal) == euler.subs(at)


@pytest.mark.parametrize("bad", [True, 0.1, sp.Float("0.1"), sp.oo, sp.nan, sp.I])
def test_exact_numeric_replay_rejects_inexact_nonfinite_and_nonreal_data(bad):
    from p8_a26_vacuum import coefficients, jets

    with pytest.raises((TypeError, ValueError)):
        coefficients.residues(bad, sp.Rational(1, 5))
    with pytest.raises((TypeError, ValueError)):
        jets.evaluate_euler((bad, 0, 0), (1, 1, 1))


@pytest.mark.parametrize("gradient", [(0, 0, 0), (1, 1, 0), (1, 2, 0)])
def test_null_and_spacelike_directions_are_not_in_this_timelike_replay(gradient):
    from p8_a26_vacuum import jets

    with pytest.raises(ValueError):
        jets.evaluate_euler(gradient, (1, 1, 1))


def test_degenerate_einstein_and_wrong_size_data_are_not_generic_domain_verdicts():
    from p8_a26_vacuum import coefficients, jets

    with pytest.raises(ValueError):
        coefficients.residues(1, sp.Rational(1, 2))
    with pytest.raises(ValueError):
        jets.evaluate_euler((1, 0), (1, 1, 1))

