"""Root audit from anisotropic coframes and continuous rational inequalities.

The reference action is expanded over subsets of four independent
coframes, rather than using the candidate's isotropic generating function.
The interval conclusion uses the written continuous proof and the bound
polynomial, not the exact fixture evaluations alone.
"""

from fractions import Fraction
from functools import cache
from itertools import combinations

import pytest
import sympy as sp
from p8_variable_vacuum import bounds, model, obstruction

F = Fraction


@cache
def anisotropic_action():
    gs = sp.symbols("Ng Ax Ay Az", positive=True)
    fs = sp.symbols("Nf Bx By Bz", positive=True)
    root = sp.Symbol("root", positive=True)
    phi = sp.Symbol("phi", real=True)
    beta = tuple(sp.Function(f"b{index}")(phi) for index in range(5))
    density = 0
    for degree in range(5):
        for selected in combinations(range(4), degree):
            density -= 2*beta[degree]*sp.prod(fs[index] if index in selected else gs[index] for index in range(4))
    at = {fs[index]: root*gs[index] for index in range(4)}
    g_eq = [sp.factor(sp.diff(density, gs[index]).subs(at)/sp.prod(gs[j] for j in range(4) if j != index)) for index in range(4)]
    f_eq = [sp.factor((sp.diff(density, fs[index])/sp.prod(fs[j] for j in range(4) if j != index)).subs(at)) for index in range(4)]
    U = beta[0]+3*root*beta[1]+3*root**2*beta[2]+root**3*beta[3]
    V = beta[1]+3*root*beta[2]+3*root**2*beta[3]+root**3*beta[4]
    clock = sp.factor((-sp.diff(density, phi)/sp.prod(gs)).subs(at))
    return {"gs": gs, "fs": fs, "r": root, "phi": phi, "beta": beta,
            "density": density, "g_equations": g_eq, "f_equations": f_eq,
            "U": U, "V": V, "clock": clock}


@pytest.mark.parametrize("index", range(4))
def test_each_independently_varied_g_and_f_coframe_equation(index):
    d = anisotropic_action()
    assert sp.expand(d["g_equations"][index]+2*d["U"]) == 0
    assert sp.expand(d["f_equations"][index]+2*d["V"]/d["r"]**3) == 0
    primary = model.derive()
    mapping = dict(zip(primary["betas"], d["beta"]))
    mapping[primary["r"]] = d["r"]
    assert sp.expand(primary["U"].subs(mapping)-d["U"]) == 0
    assert sp.expand(primary["V"].subs(mapping)-d["V"]) == 0


def test_clock_variation_has_four_dimensional_not_lapse_binomial_factors():
    d = anisotropic_action()
    root, phi, beta = d["r"], d["phi"], d["beta"]
    expected = 2*sum(sp.binomial(4, j)*root**j*sp.diff(beta[j], phi) for j in range(5))
    assert sp.expand(d["clock"]-expected) == 0
    wrong = 2*sp.diff(d["U"], phi)
    assert sp.expand(d["clock"]-wrong) != 0
    primary = model.derive()
    derivatives = {sp.diff(primary["betas"][j], primary["phi"]): sp.diff(beta[j], phi) for j in range(5)}
    assert sp.expand(primary["clock_equation"].subs(derivatives).subs(primary["r"], root)-expected) == 0


def test_dimensionful_clock_slope_is_not_a_rolling_vacuum_stress():
    M, tau, kinetic, slope = sp.symbols("M tau k slope", positive=True)
    phi_u = M*sp.sqrt(kinetic)
    beta_u = M**2*slope/tau**2
    assert sp.simplify(beta_u/phi_u-M*slope/(tau**2*sp.sqrt(kinetic))) == 0
    d = obstruction.evaluate(0, 4, 2)
    # These are labels of fixed coefficient values. Constant vacuum
    # fields have zero gradients, not the listed rolling kinetic density.
    assert d["kbar"] == sp.Rational(799, 100)
    assert d["g_flat_residual"] == -2


def test_generic_f_root_residual_identity_from_undivided_equations():
    hp, h2, y, c, root = sp.symbols("hp h2 y c root", positive=True)
    n = 2*(y**3/c-1)*hp
    b1 = y**3*hp/(c*(c-y))
    b4 = 3*h2/(2*c*c)-b1/y**3
    b0 = 3*h2/2-n/4-3*b1*y
    U, V = b0+3*b1*root, b1+b4*root**3
    denominator = root**2+root*y+y*y
    identity = U+3*y**3*V/denominator-3*h2/2+n/4-9*h2*y**3*root**3/(2*c*c*denominator)
    assert sp.factor(identity) == 0
    assert sp.expand(denominator-3*root*y-(root-y)**2) == 0
    theta = 3*h2*(c-y)/(2*c*hp)
    assert sp.factor(b4+b1*(1-theta)/y**3) == 0
    assert sp.factor(b1+b4*y**3/(1-theta)) == 0


def test_continuous_primitive_and_branch_bounds_without_root_sampling():
    dmax = F(101, 100)
    h2_upper, hp_lower, y3_lower = F(4, 25), F(3), F(7)
    assert 4*F(99, 100)/dmax**2 > hp_lower
    assert 8/dmax**12 > y3_lower
    assert 2/dmax**4 > 1
    theta_upper = F(3, 2)*h2_upper/hp_lower
    assert theta_upper == F(2, 25) < 1
    root_cube_upper = 8/(1-theta_upper)
    assert root_cube_upper == F(200, 23) < F(9, 4)**3
    # c1: r<=y, n/4>9. Positive branch: r<9/4 and
    # n/4>(3/2)(7/c-1). These hold on the full coefficient interval.
    c1_bound = F(3, 2)*h2_upper*(1+2**4)-F(3, 2)*(y3_lower-1)
    assert c1_bound == -F(123, 25)
    c = sp.Symbol("c", positive=True)
    positive = sp.Rational(F(3, 2)*h2_upper)*(c*c+4*sp.Rational(9, 4)**2)-sp.Rational(3, 2)*(7*c-c*c)
    wanted = (87*c*c-525*c+243)/50
    assert sp.expand(positive-wanted) == 0
    assert sp.factor(-sp.Rational(459, 50)-wanted-3*(c-2)*(117-29*c)/50) == 0


def test_continuous_Bernstein_conversion_is_independent_of_supplied_list():
    z, c = sp.symbols("z c", real=True)
    q = sp.Rational(3, 2)*sp.Rational(4, 25)*(c*c+4*sp.Rational(9, 4)**2)-sp.Rational(3, 2)*(7*c-c*c)
    margin = sp.Poly(sp.expand((-sp.Rational(459, 50)-q).subs(c, 2+2*z)), z)
    power = [F(margin.nth(i)) for i in range(3)]
    converted = [power[0], power[0]+power[1]/2, sum(power)]
    assert converted == [F(0), F(177, 50), F(3, 25)]
    assert converted == [F(value) for value in bounds.derive()["endpoint_margin_Bernstein"]]
    assert -F(459, 50)/16 == -F(459, 800)


@pytest.mark.parametrize("c", [1, sp.Rational(201, 100), 3, 4])
def test_center_f_and_clock_solution_still_fails_g_for_independent_root(c):
    result = obstruction.evaluate(0, c, 2)
    assert result["r_star_cubed"] == 8
    assert result["f_flat_residual"] == result["clock_profile_residual"] == 0
    assert result["g_flat_residual"] == 2-16/sp.Rational(c) < 0
    if c != 2:
        wrong = obstruction.evaluate(0, c, c)
        assert wrong["f_flat_residual"] != 0


def test_no_vacuum_claim_from_metric_equations_alone_or_changed_action():
    d = anisotropic_action()
    values = {d["r"]: 2, **dict(zip(d["beta"], [-6, 1, 0, 0, -sp.Rational(1, 8)]))}
    assert d["U"].subs(values) == d["V"].subs(values) == 0
    slopes = {sp.diff(beta, d["phi"]): int(index == 0) for index, beta in enumerate(d["beta"])}
    assert d["clock"].subs(slopes).subs(values) == 2
    at_center = obstruction.evaluate(0, 4, 2)
    assert at_center["g_flat_residual"]+2 == 0
    # A constant b0 shift alters only the g volume term -2*b0*Ng*a^3.
    ng, a = sp.symbols("Ng a", positive=True)
    extra_density = -4*ng*a**3
    assert sp.diff(extra_density, ng)/a**3 == -4


def test_strict_domain_is_checked_after_cache_warmup():
    obstruction.evaluate(0, 4, 2)
    for args in [(0, 2, 2), (sp.Rational(101, 1000), 4, 2),
                 (0, 4, 0), (0, 4, True), (0.0, 4, 2)]:
        with pytest.raises((TypeError, ValueError)):
            obstruction.evaluate(*args)
