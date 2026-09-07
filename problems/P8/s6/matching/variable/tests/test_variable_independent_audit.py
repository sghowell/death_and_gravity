"""Separately derived full-background and canonical-action audit."""

from fractions import Fraction as Q
from functools import cache
from itertools import combinations

import pytest
import sympy as sp


@cache
def _background():
    u = sp.Symbol("u", real=True)
    c = sp.Symbol("c", positive=True)
    d = 1+u*u
    a, y = d**2, 2/d**4
    h = 4*u/d
    hp = sp.diff(h, u)
    null = 2*(y**3/c-1)*hp
    chi2 = 1/(100*a**6)
    kinetic = null-chi2
    b1 = y**3*hp/(c*(c-y))
    b4 = 3*h*h/(2*c*c)-b1/y**3
    b0 = (3*h*h-null/2)/2-3*b1*y
    return u, c, a, y, h, null, chi2, kinetic, b0, b1, b4


@cache
def _equations():
    u, c, _, y, h, null, chi2, kinetic, b0, b1, b4 = _background()
    hp = sp.diff(h, u)
    rho_g, rho_f = 2*(b0+3*b1*y), 2*(b1+b4*y**3)/y**3
    interaction_g = 2*b1*(y-c)
    interaction_f = 2*b1*(c-y)/(c*y**3)
    force_u = sp.diff(b0, u)+3*y*sp.diff(b1, u)+c*(sp.diff(b1, u)+y**3*sp.diff(b4, u))
    raw = (3*h*h-null/2-rho_g,
           -(2*hp+3*h*h)-null/2-interaction_g+rho_g,
           3*h*h/c**2-rho_f,
           (2*hp-3*h*h)/c**2-interaction_f+rho_f,
           sp.diff(chi2, u)/2+3*h*chi2,
           sp.diff(kinetic, u)/2+3*h*kinetic+2*force_u,
           6*b1*(-y*h/c-h)-2*(sp.diff(b1, u)+y**3*sp.diff(b4, u)))
    return tuple(sp.factor(item) for item in raw)


def test_both_einstein_metrics_and_both_canonical_scalars_satisfy_full_equations():
    assert _equations() == (0,)*7


@pytest.mark.parametrize("lapse,expected", [(1, (178, -32, 4, sp.Rational(5599, 100))),
                                            (4, (-26, 4, -sp.Rational(1, 2), sp.Rational(799, 100)))])
def test_actual_center_data_use_the_second_metrics_proper_clock(lapse, expected):
    u, c, _, _, _, _, _, kinetic, b0, b1, b4 = _background()
    assert tuple(value.subs({u: 0, c: lapse}) for value in (b0, b1, b4, kinetic)) == expected


def test_uniform_positive_clock_box_is_exact_and_does_not_assert_full_health():
    assert Q(101, 100)**12 < Q(8, 7)
    assert Q(39600, 10201) > 3
    assert 2*(Q(7, 4)-1)*3-Q(1, 100) == Q(449, 100) > 0
    assert Q(101, 100)**4 < 2  # y>1, so the c=1 chart also avoids c=y.


def test_clock_is_interaction_sourced_while_the_extra_matter_is_truly_free():
    u, c, a, _, h, _, _, kinetic, *_ = _background()
    free_clock_error = sp.factor(sp.diff(kinetic, u)/2+3*h*kinetic)
    assert free_clock_error.subs({u: sp.Rational(1, 10), c: 3}) != 0
    chi_prime = 1/(10*a**3)
    assert sp.factor(sp.diff(chi_prime, u)+3*h*chi_prime) == 0


def test_literal_shift_block_and_kinetic_schur_complement():
    ng, ag, y, c = sp.symbols("Ng ag y c", positive=True)
    shift = sp.Symbol("shift", real=True)
    beta = sp.symbols("b0:5", real=True)
    block = sp.Matrix([[c*c-y*y*ag*ag*shift*shift/ng**2, -y*y*ag*ag*shift/ng**2],
                       [y*y*shift, y*y]])
    trace_root = sp.sqrt(sp.trace(block)+2*c*y)
    root = (block+c*y*sp.eye(2))/trace_root
    assert sp.simplify(root*root-block) == sp.zeros(2)
    es = (1, trace_root+2*y, c*y+2*y*trace_root+y*y,
          2*c*y*y+y*y*trace_root, c*y**3)
    lag = -2*ng*ag**3*sum(b*e for b, e in zip(beta, es, strict=True))
    p = 2*(beta[1]+2*beta[2]*y+beta[3]*y*y)
    actual = sp.factor(sp.diff(lag, shift, 2).subs(shift, 0)/2)
    assert sp.factor(actual-ag**5*y*y*p/(2*ng*(c+y))) == 0
    a, b, mass, v, w, sg, sf = sp.symbols("A B D v w sg sf", nonzero=True)
    kinetic = (a*(v-sg)**2+b*(w-sf)**2+mass*(sf-sg)**2)/2
    stationary = sp.solve([sp.diff(kinetic, sg), sp.diff(kinetic, sf)], (sg, sf))
    assert sp.factor(kinetic.subs(stationary)-(v-w)**2/(2*(1/mass+1/a+1/b))) == 0


def test_literal_traceless_coframe_potential_has_the_correct_polarization_norm():
    ng, ag, c, y = sp.symbols("Ng ag c y", positive=True)
    relative = sp.Symbol("relative", real=True)
    beta = sp.symbols("b0:5", real=True)
    roots = (c, y*sp.exp(relative/2), y*sp.exp(-relative/2), y)
    elementary = [sum(sp.prod(term) for term in combinations(roots, degree)) for degree in range(5)]
    lag = -2*ng*ag**3*sum(b*e for b, e in zip(beta, elementary, strict=True))
    quadratic = sp.diff(lag, relative, 2).subs(relative, 0)/2
    mu = 2*y*(beta[1]+beta[2]*(c+y)+beta[3]*c*y)
    assert sp.factor(quadratic+ng*ag**3*mu/4) == 0
    assert sp.factor(mu.subs({beta[2]: 0, beta[3]: 0})-2*y*beta[1]) == 0


def test_time_dependent_common_relative_map_retains_rotation_and_normalization():
    k1, k2 = sp.symbols("K1 K2", positive=True)
    k1p, k2p = sp.symbols("K1p K2p", real=True)
    total = k1+k2
    transform = sp.Matrix([[1/sp.sqrt(total), -sp.sqrt(k2/(k1*total))],
                           [1/sp.sqrt(total), sp.sqrt(k1/(k2*total))]])
    metric = sp.diag(k1, k2)
    derivative = sp.diff(transform, k1)*k1p+sp.diff(transform, k2)*k2p
    theta_s = (k1p+k2p)/(2*total)
    theta_r = (k1p/k1+k2p/k2)/2-theta_s
    omega = sp.sqrt(k1*k2)*(k2p/k2-k1p/k1)/(2*total)
    cross = sp.simplify(transform.T*metric*derivative)
    assert sp.simplify(transform.T*metric*transform-sp.eye(2)) == sp.zeros(2)
    assert sp.simplify(cross-sp.Matrix([[-theta_s, -2*omega], [0, -theta_r]])) == sp.zeros(2)
    assert sp.simplify(derivative.T*metric*derivative-cross.T*cross) == sp.zeros(2)


def test_c4_tensor_gradient_and_vector_sign_controls_are_not_gap_proofs():
    y, c, p = sp.Integer(2), sp.Integer(4), sp.Integer(8)
    k1, k2 = sp.S.One, y**3/c
    assert (k1+k2) == 3
    assert (1+c*y)/(k1+k2) == 3
    assert y*p*(1+c/y**3) == 24
    assert c*c/y**2 == 4
    assert 2*(8-c)/(c+8)*3 == 2  # The relative normalization pump at u=0.
    assert p/(c+y) > 0
    assert sp.Rational(-64, 1)/(1+y) < 0  # Different c=1 actual solution.


def test_nominal_mass_and_true_canonical_diagonal_share_nonvanishing_second_adiabatic_ratio():
    u, c, a, y, h, *_ = _background()
    mass = 2*y*sp.diff(h, u)*(y**3+c)/(c*(c-y))
    ksum = a**3*(1+y**3/c)
    theta = sp.factor(sp.diff(ksum, u)/(2*ksum))
    nr = -sp.diff(theta, u)+theta**2
    omega2 = (sp.Rational(3, 2)*sp.diff(y, u)/y)**2*c*y**3/(c+y**3)**2
    actual_diagonal = mass-nr-4*omega2
    expected = -c*(7*c**2+146*c-240)/(8*(c+8)**2)
    ratio = sp.factor(sp.diff(mass, u, 2).subs(u, 0)/mass.subs(u, 0)**2)
    assert sp.factor(ratio-expected) == 0
    assert sp.limit(ratio, c, 2, dir="+") == -sp.Rational(1, 5)
    actual_ratio = sp.factor(sp.diff(actual_diagonal, u, 2).subs(u, 0)/actual_diagonal.subs(u, 0)**2)
    assert sp.limit(actual_ratio, c, 2, dir="+") == -sp.Rational(1, 5)
    assert sp.diff(actual_diagonal, u).subs(u, 0) == 0


def test_inner_operator_and_constant_source_have_a_finite_nonadiabatic_defect():
    x = sp.Symbol("x", real=True)
    delta = sp.Symbol("delta", positive=True)
    u, c, _, y, h, *_ = _background()
    mass = 2*y*sp.diff(h, u)*(y**3+c)/(c*(c-y))
    limit = sp.limit(delta*mass.subs({u: sp.sqrt(delta)*x, c: 2+delta}), delta, 0, dir="+")
    assert sp.factor(limit-80/(1+8*x*x)) == 0
    algebraic = (1+8*x*x)/80
    particular = (1+8*x*x)/96
    assert sp.simplify(sp.diff(algebraic, x, 2)+limit*algebraic-1) == sp.Rational(1, 5)
    assert sp.simplify(sp.diff(particular, x, 2)+limit*particular-1) == 0
    assert sp.cancel(particular/algebraic) == sp.Rational(5, 6)


@pytest.mark.parametrize("time,lapse", [(-sp.Rational(1, 10), 1), (0, 1), (0, 4),
                                        (sp.Rational(1, 10), 3), (sp.Rational(1, 20), sp.Rational(201, 100))])
def test_primary_actual_family_matches_independent_physical_clock_values(time, lapse):
    from p8_variable_beta import background

    u, c, a, y, h, null, chi2, kinetic, b0, b1, b4 = _background()
    at = {u: time, c: lapse}
    expected = {"a": a, "b": a*y, "y": y, "h": h, "h_u": sp.diff(h, u),
                "h_f": -h/c, "D_f_h_f": -sp.diff(h, u)/c**2,
                "nbar": null, "kbar": kinetic, "b0": b0, "b1": b1, "b4": b4, "P": 2*b1}
    actual = background.evaluate(time, lapse)
    for key, value in expected.items():
        assert actual[key] == sp.cancel(value.subs(at))
    assert actual["chi_speed"]**2 == chi2.subs(at)
    assert actual["kbar"] > 0


def test_primary_canonical_data_retain_actual_derivatives_and_mass_not_just_center_numbers():
    from p8_variable_beta import canonical

    u, c, a, y, h, *_ = _background()
    primary = canonical.derive()
    total = a**3*(1+y**3/c)
    theta = sp.factor(sp.diff(total, u)/(2*total))
    expected = {"theta_sum": theta, "theta_relative": -theta,
                "N_sum": sp.diff(theta, u)+theta**2,
                "N_relative": -sp.diff(theta, u)+theta**2,
                "mass_squared": 2*y*sp.diff(h, u)*(y**3+c)/(c*(c-y)),
                "c_light_squared": (1+c*y)/(1+y**3/c)}
    for key, value in expected.items():
        actual = primary[key].subs({primary["u"]: u, primary["c"]: c})
        assert sp.factor(actual-value) == 0
    assert sp.factor(primary["f_sum_squared"]-2*primary["K_sum"]) == 0
    assert sp.factor(primary["f_relative_squared"]-2*primary["K_relative"]) == 0


@pytest.mark.parametrize("lapse", [0, 2, sp.Rational(3, 2), 5, True, 0.25, sp.oo, sp.nan])
def test_outside_family_or_inexact_lapse_is_not_silently_admitted(lapse):
    from p8_variable_beta import background

    with pytest.raises((TypeError, ValueError)):
        background.evaluate(0, lapse)


@pytest.mark.parametrize("time", [-sp.Rational(11, 100), sp.Rational(11, 100), True, 0.1, sp.oo])
def test_finite_local_window_and_exact_clock_guards(time):
    from p8_variable_beta import background

    with pytest.raises((TypeError, ValueError)):
        background.evaluate(time, 3)

