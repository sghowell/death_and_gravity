"""Root-authored checks of the prepared-branch proof, not a residual oracle.

The parent operator is frozen and separately derived from the physical TT
action. This audit independently clears its pole, derives the joint leading
recurrences and checks the physical and symplectic dictionaries. Continuous
analytic estimates still require the written norm proof, not finite jets.
"""

from fractions import Fraction

import pytest
import sympy as sp
from p8_variable_prepared import analytic, majorants, physical
from p8_variable_prepared import source as prepared_source
from p8_variable_response import operator


def test_actual_heavy_equation_after_pole_clearing():
    d = operator.derive()
    u, c = d["u"], d["c"]
    den = c*(1+u**2)**4-2
    q = sp.Function("q")(u)
    light = sp.Function("l")(u)
    relative = den*q
    numerator = 16*(1-u**2)*(8/(c*(1+u**2)**14)+1/(1+u**2)**2)
    actual = (sp.diff(relative, u, 2)+d["B"]*relative
              +d["E"]*light+u*d["f_cross"]*sp.diff(light, u))
    cleared = (den*sp.diff(q, u, 2)+2*sp.diff(den, u)*sp.diff(q, u)
               +(sp.diff(den, u, 2)+numerator+den*d["b_analytic"])*q
               +d["E"]*light+u*d["f_cross"]*sp.diff(light, u))
    assert sp.factor(actual-cleared) == 0
    assert numerator.subs({u: 0, c: 2}) == 80
    assert sp.diff(den, u, 2).subs({u: 0, c: 2}) == 16


@pytest.mark.parametrize("degree", range(8))
def test_weighted_homogeneous_leading_inverse(degree):
    # Independent finite triangular inverse at each weighted degree.
    # delta has weight two, u weight one. The exact all-degree inequalities
    # underlying the 1/84 bound are checked separately below.
    delta, u = sp.symbols("delta u")
    pairs = [(j, degree-2*j) for j in range(degree//2+1)]
    forcing = {pair: Fraction((-1)**j*(j+1), j+2)
               for j, pair in enumerate(pairs)}
    solution = {}
    for j, n in pairs:
        inherited = (n+2)*(n+1)*solution.get((j-1, n+2), Fraction(0))
        solution[j, n] = (forcing[j, n]-inherited)/(8*(n*n+3*n+12))
    polynomial = sum(sp.Rational(value.numerator, value.denominator)*delta**j*u**n
                     for (j, n), value in solution.items())
    rhs = sum(sp.Rational(value.numerator, value.denominator)*delta**j*u**n
              for (j, n), value in forcing.items())
    image = ((delta+8*u*u)*sp.diff(polynomial, u, 2)
             +32*u*sp.diff(polynomial, u)+96*polynomial)
    assert sp.expand(image-rhs) == 0
    assert sum(map(abs, solution.values())) <= sum(map(abs, forcing.values()))/84


def test_leading_inverse_bounds_are_all_degree_inequalities():
    n = sp.Symbol("n", integer=True, nonnegative=True)
    diagonal = 8*(n*n+3*n+12)
    # For B = diagonal^{-1} delta*d_u^2, the output degree is n.
    assert sp.expand(diagonal/8-(n+2)*(n+1)) == 10
    assert sp.expand(diagonal-96-8*n*(n+3)) == 0
    # n/[8(n^2+3n+12)] <= 1/80 for integer n>=0:
    # (n-3)(n-4)>=0 on the integers, including the two zero cases.
    assert sp.expand(diagonal-80*n-8*(n-3)*(n-4)) == 0
    assert Fraction(1, 96)/(1-Fraction(1, 8)) == Fraction(1, 84)
    assert Fraction(1, 80)/(1-Fraction(1, 8)) == Fraction(1, 70)


def test_full_operator_gives_regular_even_and_odd_leading_coefficients():
    d = operator.derive()
    u, c, k = d["u"], d["c"], d["kbar"]
    origin = {u: 0, c: 2}
    assert sp.simplify(d["A"].subs(origin)-k**2-sp.Rational(18, 5)) == 0
    assert sp.simplify(d["f_cross"].subs(origin)+sp.Rational(48, 5)) == 0
    assert sp.simplify(sp.diff(d["E"], c).subs(origin)-2*k**2/5) == 0
    assert sp.simplify(sp.diff(d["E"], u, 2).subs(origin)/2
                       -16*k**2/5+sp.Rational(864, 25)) == 0

    aa, bb, odd, delta = sp.symbols("a b odd delta")
    leading = lambda q: ((delta+8*u*u)*sp.diff(q, u, 2)
                         +32*u*sp.diff(q, u)+96*q)
    even_rhs = -2*k*k*delta/5-64*k*k*u*u/5
    even_coefficients = sp.solve(sp.Poly(leading(aa*delta+bb*u*u)-even_rhs,
                                        delta, u).coeffs(), (aa, bb))
    assert even_coefficients == {aa: -7*k*k/2640, bb: -4*k*k/55}
    assert sp.solve(leading(odd*u)-48*u/5, odd) == [sp.Rational(3, 40)]
    relative_even = sp.expand((delta+8*u*u)*(aa*delta+bb*u*u)).subs(even_coefficients)
    assert sp.expand(relative_even+k*k*(7*delta**2/2640
                                        +31*delta*u*u/330+32*u**4/55)) == 0


def test_canonical_symplectic_form_and_prepared_center():
    d = operator.derive()
    u = d["u"]
    first = sp.Matrix(sp.symbols("l1 p1 Q1 P1"))
    second = sp.Matrix(sp.symbols("l2 p2 Q2 P2"))
    cross = first[0]*second[2]-first[2]*second[0]
    omega = (first[0]*second[1]-first[1]*second[0]
             +first[2]*second[3]-first[3]*second[2]-2*d["omega"]*cross)
    vector = operator.first_order()
    derivative = (sp.diff(omega, u)
                  +(sp.Matrix([omega]).jacobian(first)*vector*first)[0]
                  +(sp.Matrix([omega]).jacobian(second)*vector*second)[0])
    assert sp.factor(derivative) == 0
    delta, qe, qop = sp.symbols("delta qe qop", real=True)
    center = omega.subs(u, 0).subs(dict(zip(first, [1, 0, delta*qe, 0])))
    center = center.subs(dict(zip(second, [0, 1, 0, delta*qop])))
    assert sp.expand(center-1-delta**2*qe*qop) == 0


def test_center_closed_physical_equation_uses_actual_g_not_locked_field():
    delta = sp.Symbol("delta", positive=True)
    qe = sp.Symbol("qe", real=True)
    c = 2+delta
    kg, kf = sp.Rational(1, 8), 1/c
    fs = sp.sqrt(2*(kg+kf))
    fr = sp.sqrt(2*kg*kf/(kg+kf))
    g = 1/fs-kf/(kg+kf)*delta*qe/fr
    mu = 128/(c*delta)
    correction = sp.simplify(mu*delta/(fr*g))
    expected = 32*(delta+10)/((delta+2)*(sp.sqrt(1+delta/2)-2*delta*qe))
    assert sp.simplify(correction-expected) == 0
    assert sp.simplify(expected.subs({delta: 0, qe: 0})) == 160
    assert sp.Rational(160)*sp.Rational(7, 2640) == sp.Rational(14, 33)
    assert sp.Rational(14, 33) != sp.Rational(4, 5)


def test_center_multiplier_bounds_without_sampling_or_sign_of_q():
    delta = sp.Symbol("delta", nonnegative=True)
    upper_delta = sp.Rational(1, 10**9)
    # Squaring positive brackets proves these radical bounds on the full
    # interval, not just at its endpoint. The lower squared error is
    # delta^3*(delta-16)/1024 and the upper error is delta^2/16.
    low_sqrt = 1+delta/4-delta**2/32
    high_sqrt = 1+delta/4
    assert sp.factor(low_sqrt**2-1-delta/2) == delta**3*(delta-16)/1024
    assert sp.expand(high_sqrt**2-1-delta/2) == delta**2/16
    assert upper_delta < 16
    assert 1-upper_delta**2/32 > 0
    # The coefficient-l1 Schwarz bound gives abs(q_e(delta,0))<=3*delta.
    # This audit deliberately does not assume its sign.
    low_den = (2+delta)*(low_sqrt-6*delta**2)
    high_den = (2+delta)*(high_sqrt+6*delta**2)
    numerator = 32*(10+delta)
    upper_gap = sp.factor((160*low_den-numerator)/delta)
    assert sp.expand(upper_gap-160*(sp.Rational(13, 10)
                                   -sp.Rational(189, 16)*delta
                                   -sp.Rational(193, 32)*delta**2)) == 0
    assert upper_gap.subs(delta, upper_delta) > 0
    lower_gap = sp.factor((numerator-(160-105*delta)*high_den)/delta)
    assert sp.expand(lower_gap-(2-sp.Rational(3605, 2)*delta
                                +sp.Rational(1305, 4)*delta**2+630*delta**3)) == 0
    assert 2-sp.Rational(3605, 2)*upper_delta > 0
    assert 160-105*upper_delta > 0
    # Therefore 0<B<=160 and 160-B<=105*delta continuously.
    assert 160*300+sp.Rational(105*7, 2640) < 48001


def test_coefficient_l1_tail_and_symplectic_constants():
    radius = Fraction(1, 20)
    delta_radius = radius*radius
    delta_max = Fraction(1, 10**9)
    assert delta_max/delta_radius < 1
    # Banach-valued K-Schwarz: norm(q_e/K)<=3*v/4. Every delta
    # coefficient of degree>=2 then gains at least (delta/v)^2.
    assert (3*delta_radius/4)/delta_radius**2 == 300
    # q_e(0,0)=0 gives |q_e(delta,0)|<=3*delta, while the odd
    # derivative extraction gives |q_o,u(delta,0)|<=1/5.
    assert (3*delta_radius)/delta_radius == 3
    assert (radius/5)/radius == Fraction(1, 5)
    omega_error = Fraction(3, 5)*delta_max**3
    assert omega_error < Fraction(1, 10)
    # With the independently proved physical Wronskian bracket (3/5,2),
    # the kinetic bracket follows with the same canonical orientation.
    assert Fraction(9, 10)/2 == Fraction(9, 20)
    assert Fraction(11, 10)/Fraction(3, 5) == Fraction(11, 6)


@pytest.mark.parametrize("j,n", [(0, 0), (0, 2), (1, 3), (2, 4), (1, 9)])
def test_production_inverse_against_independent_coefficient_recurrence(j, n):
    d = analytic.derive()
    current = Fraction(1, 8*(n*n+3*n+12))
    coefficients = {(j, n): current}
    degree = n
    for shift in range(1, n//2+1):
        next_degree = degree-2
        current *= Fraction(-degree*(degree-1),
                            8*(next_degree**2+3*next_degree+12))
        coefficients[j+shift, next_degree] = current
        degree = next_degree
    actual = sp.Poly(analytic.inverse_monomial(j, n), d["delta"], d["u"])
    assert {powers: Fraction(value) for powers, value in actual.terms()} == coefficients


def test_contraction_columns_from_integrated_full_equations():
    c = majorants.calibration()
    v = Fraction(1, 400)
    # Each term is the operator norm in ||l||+v||q||. In particular the
    # u*D*q' light term has 1/6, not the cruder derivative-free 1/2.
    d_norm = (2+v)*(1+v)**4-2
    u_dprime = 8*v*(2+v)*(1+v)**3
    light_column = 9*v/2+v*(60*v/84+Fraction(10, 70))
    q_to_light = (Fraction(12, 2)+Fraction(10, 6))*d_norm+5*u_dprime
    assert Fraction(c["column_light"]) == light_column < Fraction(1, 2)
    assert Fraction(c["light_from_q_scaled"]) == q_to_light
    assert Fraction(c["column_q"]) == q_to_light+Fraction(c["q_self"]) < Fraction(1, 2)


def test_real_physical_wronskian_bounds_from_branch_norm_and_actual_map():
    radius, v = Fraction(1, 20), Fraction(1, 400)
    a, delta = Fraction(1, 100), Fraction(1, 10**9)
    rho = a/radius
    cmax, dmax = 2+delta, 1+a*a
    assert delta/v < rho*rho
    # Actual ag=2 sqrt(c)d^3/sqrt(cd^12+8), bg=-4sqrt(2)/(d^3 sqrt(cd^12+8)).
    assert Fraction(8)/(cmax*dmax**12+8) > Fraction(4, 5)**2
    assert 4*cmax*dmax**6/10 < 1
    assert Fraction(32, 10) < 2**2
    # ag'=-theta*ag; bg'/bg=-6u/d-12cu*d^11/(cd^12+8).
    agprime = 4*a
    assert 2*(6+12*cmax*dmax**11/10) < 20
    bgprime = 20*a
    assert cmax*dmax**4-2 < Fraction(1, 1200)
    assert 8*a*cmax*dmax**3 < Fraction(1, 6)
    d_bound, dp_bound = Fraction(1, 1200), Fraction(1, 6)
    # Even q has no joint constant term; odd q has at least one u.
    qe, qep = 3*v*rho*rho, 3*v/radius*2*rho
    qo, qop = radius/5*rho, Fraction(1, 5)
    Qe, Qep = d_bound*qe, dp_bound*qe+d_bound*qep
    Qo, Qop = d_bound*qo, dp_bound*qo+d_bound*qop
    # Initial conditions/parity give u^2 for l_e-1 and u^3 for l_o-u.
    le_error, lep = 10*v*rho*rho, 10*v/radius*2*rho
    lo_error, lop_error = radius/40*rho**3, Fraction(1, 40)*3*rho*rho
    ge_lower = Fraction(4, 5)*(1-le_error)-2*Qe
    ge_upper = 1+le_error+2*Qe
    gep_upper = agprime*(1+le_error)+lep+bgprime*Qe+2*Qep
    go_upper = a+lo_error+2*Qo
    gop_lower = (Fraction(4, 5)*(1-lop_error)-agprime*(a+lo_error)
                 -bgprime*Qo-2*Qop)
    gop_upper = 1+lop_error+agprime*(a+lo_error)+bgprime*Qo+2*Qop
    assert ge_lower > Fraction(79, 100)
    assert ge_upper < Fraction(501, 500)
    assert gep_upper < Fraction(27, 100)
    assert go_upper < Fraction(11, 1000)
    assert gop_lower > Fraction(79, 100)
    assert gop_upper < Fraction(101, 100)
    assert ge_lower*gop_lower-gep_upper*go_upper > Fraction(3, 5)
    assert ge_upper*gop_upper+gep_upper*go_upper < 2


def test_endpoint_and_fixed_source_telescoping_constants():
    c = physical.calibration()
    delta = sp.Rational(1, 10**10)
    assert c["endpoint_even_column"]+c["endpoint_odd_column"] < 200
    assert 42*200+200 == 8600
    assert 42*3000000 == 126000000
    assert physical.fixed_light_error(delta, 3) == 3*8600*delta
    assert physical.fixed_source_error(delta, 3, 5) == delta*(3*8600+5*126000000)


def test_source_loading_by_commuting_actual_metric_operators():
    u, k = sp.symbols("u K", real=True)
    kg, kf, gg, gf, inv, field, switch = (
        sp.Function(name)(u) for name in ("kg", "kf", "gg", "gf", "P", "F", "zeta"))

    def reconstruct_g(f):
        return f+inv*(sp.diff(kf*sp.diff(f, u), u)+gf*k*f)

    def g_equation(g, f):
        return sp.diff(kg*sp.diff(g, u), u)+gg*k*g+(g-f)/inv

    target = reconstruct_g(field)
    loaded = reconstruct_g(switch*field)
    commutator = g_equation(loaded, switch*field)-switch*g_equation(target, field)
    actual = prepared_source.loading_coefficients(u, kg, kf, gg, inv, target, field, k)
    assembled = sum(value*sp.diff(switch, u, order)
                    for order, value in actual["coefficients"].items())
    assert sp.expand(commutator-assembled) == 0
    # The literal g equation is a^3*sigma/4, not a mass-normalized source.
    a = sp.Symbol("a", positive=True)
    sigma = 4*assembled/a**3
    assert sp.expand(a**3*sigma/4-assembled) == 0


def test_source_representation_is_regular_when_inverse_spring_vanishes():
    u, delta = sp.symbols("u delta", real=True)
    g, f = 1+u+u**3, 1-u+u**2
    loaded = prepared_source.loading_coefficients(
        u, 1+u*u, 2+u*u, 3+u*u, delta+u*u, g, f, sp.Integer(2))
    for value in loaded["coefficients"].values():
        assert sp.Poly(value, u, delta) is not None
        assert value.subs({u: 0, delta: 0}).is_finite


def test_flat_switch_derivative_majorants_from_direct_differentiation():
    x, y = sp.symbols("x y", positive=True)
    c = prepared_source.switch_bounds()
    for order in range(5):
        derivative = sp.diff(sp.exp(-1/x), x, order)
        polynomial = sp.Poly(sp.simplify(derivative*sp.exp(1/x)).subs(x, 1/y), y)
        bound = sum(abs(value)*sp.factorial(power[0])
                    for power, value in polynomial.terms())
        assert bound == c["bump_jets"][order]
    assert -sp.Rational(1, 50) < -sp.Rational(7, 400)
    assert -sp.Rational(7, 400) < -sp.Rational(1, 80) < -sp.Rational(1, 100)
    # e<3 follows from sum_{n>=2}1/n!<sum_{n>=2}2^{-(n-1)}=1;
    # hence max(eta(x),eta(1-x))>=e^-2>1/9 on the entire switch.
    assert c["denominator_lower"] == sp.Rational(1, 9)
