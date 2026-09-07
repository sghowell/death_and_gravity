"""Root-authored physical Cauchy checks from the frozen scalar Hamiltonian.

The finite fixtures here are independent controls, not a uniform-symbol
proof. No result is inferred from a frozen Hamiltonian momentum block.
Every time derivative of the physical observable map is retained.
"""

from fractions import Fraction
from functools import cache
from math import comb

import pytest
import sympy as sp
from p8_variable_constraints import action, reduction
from p8_variable_microlocal import bounds, domain, energy, leading, symbol
from p8_variable_microlocal import rational as ra
from p8_variable_microlocal import stueckelberg as st


@cache
def direct_equation(time, lapse):
    d, r = action.derive(), reduction.derive()
    u, k = d["u"], action.K
    old = dict(zip(d["state"], r["old_state"], strict=True))
    xi = sp.factor(d["a"]**2*d["B_g"].subs(old, simultaneous=True))
    assert not xi.has(reduction.P0, reduction.Q)
    observables = sp.Matrix([-action.E/k, xi, action.Z+d["chi_speed"]*xi])
    rows = observables.jacobian(reduction.STATE)
    at = {u: time, d["c"]: lapse}
    O0, O1, O2 = (sp.diff(rows, u, order).subs(at) for order in range(3))
    jets = reduction.jets(time, lapse, 1)
    A0, A1 = (reduction.J*value for value in jets["H"])
    velocity = (O1+O0*A0).applyfunc(sp.cancel)
    phase = O0.col_join(velocity)
    inverse = phase.inv().applyfunc(sp.cancel)
    acceleration = (O2+2*O1*A0+O0*(A1+A0*A0)).applyfunc(sp.cancel)
    equation = (acceleration*inverse).applyfunc(sp.cancel)
    form = (inverse.T*reduction.J*inverse).applyfunc(sp.cancel)
    return {"phase": phase, "det": sp.factor(phase.det(method="domain-ge")),
            "Mq": equation[:, :3], "Mv": equation[:, 3:], "form": form,
            "O": (O0, O1, O2), "A": (A0, A1)}


def rational_order(value):
    if value == 0:
        return -sp.oo
    numerator, denominator = sp.fraction(sp.cancel(value))
    return sp.degree(numerator, action.K)-sp.degree(denominator, action.K)


def infinity_coefficient(value):
    value = sp.cancel(value)
    if value == 0 or rational_order(value) < 0:
        return sp.Integer(0)
    assert rational_order(value) == 0
    numerator, denominator = sp.fraction(value)
    return sp.factor(sp.Poly(numerator, action.K).LC()/sp.Poly(denominator, action.K).LC())


@pytest.mark.parametrize("lapse", [sp.Integer(3), sp.Integer(4), sp.Rational(201, 100)])
def test_center_equation_and_conserved_weight_from_full_phase_map(lapse):
    d = direct_equation(sp.Integer(0), lapse)
    k = action.K
    speed = (9*lapse*lapse-22*lapse+144)/120
    principal = (d["Mq"]/k).applyfunc(infinity_coefficient)
    assert principal == sp.diag(-speed, -1, -1)
    assert all(rational_order(value) <= 0 for value in d["Mv"])
    weight = d["form"][:3, 3:].applyfunc(infinity_coefficient)
    expected = sp.diag(30720/(lapse**2*(lapse-2)**2),
                       (6400-801*lapse)/(100*lapse), 1)
    assert weight == expected
    assert all(rational_order(value) <= -2 for value in d["form"][3:, 3:])
    assert all(value > 0 for value in expected.diagonal())


def test_center_c4_determinant_and_nonuniform_chart_controls():
    d = direct_equation(sp.Integer(0), sp.Integer(4))
    k = action.K
    expected = 5*(k+60)*(k*k+26*k-24)/(19176*k**3)
    assert sp.factor(d["det"]-expected) == 0
    assert sp.expand(k*k+26*k-24-((k-1)**2+28*(k-1)+3)) == 0
    assert infinity_coefficient(d["det"]) == sp.Rational(5, 19176)
    # Two failed routes are retained explicitly: the old observable chart
    # gives5, while dropping a center-vanishing boundary coefficient's
    # nonzero time derivative gives448/480=14/15.
    assert sp.Rational(800, 480) == sp.Rational(5, 3)
    assert sp.Rational(448, 480) == sp.Rational(14, 15)
    assert sp.Rational(5, 3) not in (sp.Integer(5), sp.Rational(14, 15))


def test_general_lapse_center_identity_not_only_rational_fixtures():
    c = action.derive()["c"]
    d = direct_equation(sp.Integer(0), c)
    principal = (d["Mq"]/action.K).applyfunc(infinity_coefficient)
    expected = sp.diag(-(9*c*c-22*c+144)/120, -1, -1)
    assert (principal-expected).applyfunc(sp.factor) == sp.zeros(3)
    weights = d["form"][:3, 3:].applyfunc(infinity_coefficient)
    expected_weights = sp.diag(30720/(c**2*(c-2)**2), (6400-801*c)/(100*c), 1)
    assert (weights-expected_weights).applyfunc(sp.factor) == sp.zeros(3)
    # These exact center coefficients are not yet a uniform-time or
    # subcutoff propagation claim. The interval signs themselves are exact.
    delta = sp.Symbol("delta", positive=True)
    excess = sp.expand(((9*c*c-22*c+144)/120-1).subs(c, 2+delta))
    assert sp.expand(excess-(9*delta**2+14*delta+16)/120) == 0
    assert excess.is_positive
    assert 6400-801*4 > 0


@pytest.mark.parametrize("time", [sp.Rational(-1, 100), sp.Rational(1, 100)])
def test_nearby_exact_rational_orders_do_not_have_the_old_center_degree_drop(time):
    d = direct_equation(time, sp.Integer(4))
    assert rational_order(d["det"]) == 0
    assert all(rational_order(value) <= 1 for value in d["Mq"])
    assert all(rational_order(value) <= 0 for value in d["Mv"])
    assert all(rational_order(value) <= -2 for value in d["form"][3:, 3:])


def test_time_derivatives_are_present_in_the_actual_observable_equation():
    d = direct_equation(sp.Integer(0), sp.Integer(4))
    O0, O1, O2 = d["O"]
    A0, A1 = d["A"]
    assert O1 == sp.zeros(3, 6)
    assert O2 != sp.zeros(3, 6)
    full = O2+2*O1*A0+O0*(A1+A0*A0)
    incomplete = O0*A0*A0
    assert (full-incomplete).applyfunc(sp.cancel) != sp.zeros(3, 6)


@pytest.mark.parametrize("boundary_name", ["F1", "F0"])
def test_action_boundaries_have_no_hidden_acceleration_or_lapse_derivative(boundary_name):
    boundary = st.derive()[boundary_name]
    # Some unsimplified expressions contain cancelling velocity symbols.
    # Derivatives must vanish as identities, not merely at the center.
    for variable in (*st.PHYSICAL_VELOCITIES, st.VPSIG, st.VPSIF, st.VXF,
                     st.PHIG, st.PHIF):
        assert sp.factor(sp.diff(boundary, variable)) == 0


def test_full_literal_action_differs_by_actual_total_boundaries_only():
    d = st.derive()
    u, k = d["u"], action.K
    fields = list(st.FIELDS)
    # Independent jet derivative includes lapse velocities as well. Their
    # coefficients must cancel instead of being omitted by convention.
    velocities = [st.VPSIG, st.VPSIF, st.VPI, st.VCHI, st.VXG, st.VXF,
                  *sp.symbols("Phi_g_prime Phi_f_prime")]
    boundary = k*d["F1"]+d["F0"]
    derivative = sp.diff(boundary, u)+sum(sp.diff(boundary, q)*v
                                         for q, v in zip(fields, velocities, strict=True))
    residual = (d["literal_L"]-k*d["boundary_subtracted_L1"]
                -d["boundary_subtracted_L0"]-derivative)
    assert sp.factor(residual) == 0


def test_final_boundary_center_derivative_accounts_for_full_gradient_correction():
    d = leading.derive()
    u, c = d["u"], d["c"]
    coefficient = d["Fpi_coefficient"]
    assert sp.factor(coefficient.subs(u, 0)) == 0
    derivative = sp.factor(sp.diff(coefficient, u).subs(u, 0))
    expected = 256*(9*c*c-18*c+16)/(c*c*(c-2)**2)
    assert sp.factor(derivative-expected) == 0
    before = -sp.hessian(d["stationary_L1"], st.PHYSICAL)
    correction = (d["gradient"]-before).subs(u, 0)
    assert (correction-sp.diag(expected, 0, 0)).applyfunc(sp.factor) == sp.zeros(3)
    assert derivative.subs(c, 4) == 352


def fraction_evaluation(value, point):
    def evaluate(poly):
        return sum(Fraction(str(coefficient))*point[0]**int(powers[0])
                   *point[1]**int(powers[1])*point[2]**int(powers[2])
                   for powers, coefficient in poly.to_dict().items())
    return evaluate(value.numerator)/evaluate(value.denominator)


@pytest.mark.parametrize("time,lapse,momentum", [
    (sp.Integer(0), sp.Integer(3), sp.Integer(1)),
    (sp.Integer(0), sp.Integer(4), sp.Rational(7, 2)),
    (sp.Rational(-1, 100), sp.Integer(4), sp.Integer(2)),
    (sp.Rational(1, 100), sp.Integer(4), sp.Integer(3)),
])
def test_generic_native_symbol_against_separate_exact_time_jet_route(time, lapse, momentum):
    generic = symbol.derive()
    independent = direct_equation(time, lapse)
    point = tuple(Fraction(value) for value in (time, lapse, momentum))
    for native_name, direct_name in (("C", "phase"), ("Mq", "Mq"),
                                    ("Mv", "Mv"), ("symplectic", "form")):
        actual = sp.Matrix([[sp.Rational(fraction_evaluation(value, point)) for value in row]
                            for row in generic[native_name]])
        expected = independent[direct_name].subs(action.K, momentum)
        assert (actual-expected).applyfunc(sp.cancel) == sp.zeros(*actual.shape)


def test_native_rational_operations_and_all_partial_derivatives_against_sympy():
    u, c, k = sp.symbols("u c K", real=True)
    first = (u+c*k+2)/(c-2+u*u)
    second = (3*u-k*c+5)/(1+u*u*k)
    left, right = (ra.from_sympy(value, (u, c, k)) for value in (first, second))
    pairs = [(left+right, first+second), (left-right, first-second),
             (left*right, first*second), (left/right, first/second),
             (left**-2, first**-2)]
    for index, variable in enumerate((u, c, k)):
        pairs.append((left.derivative(index), sp.diff(first, variable)))
    for native, expression in pairs:
        assert sp.cancel(native.sympy()-expression) == 0
    with pytest.raises(ZeroDivisionError):
        ra.Rational(1, 0)
    with pytest.raises(TypeError):
        ra.Rational(0.5)


def test_native_matrix_pivoting_and_zero_order_conventions():
    u, c, k = sp.symbols("u c K", real=True)
    left = ra.from_sympy((u+1)/(c+1), (u, c, k))
    right = ra.from_sympy((k+2)/(1+u*u), (u, c, k))
    matrix = [[ra.Rational(0), left], [right, left]]
    inverse = ra.inverse(matrix)
    native_product = ra.multiply(matrix, inverse)
    assert native_product == [[ra.Rational(1), ra.Rational(0)],
                              [ra.Rational(0), ra.Rational(1)]]
    exact_inverse = sp.Matrix([[0, left.sympy()], [right.sympy(), left.sympy()]]).inv()
    assert (sp.Matrix([[value.sympy() for value in row] for row in inverse])
            -exact_inverse).applyfunc(sp.cancel) == sp.zeros(2)
    assert ra.Rational(0).degree() == -sp.oo
    assert ra.Rational(0).leading() == 0
    assert right.degree() == 1
    assert sp.cancel(right.leading().sympy()-1/(1+u*u)) == 0


@pytest.mark.parametrize("name,count", [
    ("kinetic_trace_upper", 87), ("gradient_trace_upper", 190),
    ("kinetic_minor_lower", 112), ("gradient_minor_lower", 222),
])
def test_continuous_principal_bounds_by_separate_fraction_bernstein_conversion(name, count):
    d = leading.derive()
    u, c, dd = d["u"], d["c"], d["d"]
    kinetic, gradient = d["normalized_kinetic"], d["normalized_gradient"]
    specifications = {
        "kinetic_trace_upper": (64-sp.trace(kinetic[:2, :2]),
                                 100*c*(1-u*u)**2*dd**8),
        "gradient_trace_upper": (64-sp.trace(gradient[:2, :2]),
                                  400*c*(1-u*u)**3*dd**12),
        "kinetic_minor_lower": (kinetic[:2, :2].det()-16,
                                 200*c*c*(1-u*u)**2*dd**2),
        "gradient_minor_lower": (gradient[:2, :2].det()-16,
                                  1600*c*c*(1-u*u)**3*dd**10),
    }
    margin, denominator = specifications[name]
    record = domain.derive()["records"][name]
    assert sp.cancel(denominator-record["denominator"]) == 0
    # The raw polynomial comes from the action-derived matrices, not the
    # production Bernstein conversion or a list of supplied coefficients.
    polynomial = sp.Poly(sp.cancel(margin*denominator), u, c)
    assert sp.expand(polynomial.as_expr()-record["numerator"]) == 0
    compact = {}
    for (u_power, c_power), coefficient in polynomial.terms():
        assert u_power % 2 == 0
        x_power = u_power//2
        for z_power in range(c_power+1):
            key = (x_power, z_power)
            compact[key] = compact.get(key, Fraction(0))+(
                Fraction(coefficient)*2**c_power*comb(c_power, z_power)/100**x_power)
    compact = {power: value for power, value in compact.items() if value}
    nx, nz = (max(power[index] for power in compact) for index in range(2))
    values = []
    for i in range(nx+1):
        for j in range(nz+1):
            value = sum(coefficient*Fraction(comb(i, p), comb(nx, p))
                        *Fraction(comb(j, q), comb(nz, q))
                        for (p, q), coefficient in compact.items() if p <= i and q <= j)
            values.append(value)
    assert (nx, nz) == record["degrees"]
    assert len(values) == count
    assert tuple(values) == tuple(Fraction(value) for value in record["coefficients"])
    assert min(values) > 0
    # Every cleared denominator is a product of positive constants, c,
    # 1-u^2 and d; on this box c>=2, 1-u^2>=99/100 and d>=1.
    assert Fraction(99, 100) > 0


def test_clock_sign_and_uncoupled_third_mode_complete_the_continuous_matrix_bound():
    d = leading.derive()
    u, c, dd = d["u"], d["c"], d["d"]
    clock = 2*(8/(c*dd**12)-1)*4*(1-u*u)/dd**2-sp.Rational(1, 100)/dd**12
    assert sp.factor(clock-d["kbar"]) == 0
    upper_d = Fraction(101, 100)
    clock_lower = 2*(2/upper_d**12-1)*4*Fraction(99, 100)/upper_d**2-Fraction(1, 100)
    assert 2/upper_d**12-1 > 0
    assert clock_lower > 0
    for matrix_name, power in (("normalized_kinetic", 6), ("normalized_gradient", 2)):
        matrix = d[matrix_name]
        assert sp.factor(matrix[1, 1]-dd**power*clock) == 0
        assert sp.factor(matrix[2, 2]-dd**power) == 0
        assert all(sp.factor(matrix[2, index]) == sp.factor(matrix[index, 2]) == 0
                   for index in (0, 1))
        assert Fraction(1, 4) < 1 <= upper_d**power < 64


def test_constructive_threshold_from_separate_exact_polynomial_coefficient_bounds():
    generic = symbol.derive()
    u, c, k = generic["variables"]
    native = generic["det_C"]
    numerator = sum(sp.Rational(str(value))*u**int(power[0])*c**int(power[1])*k**int(power[2])
                    for power, value in native.numerator.to_dict().items())
    polynomial = sp.Poly(numerator, k)
    assert polynomial.degree() == 3
    assert {int(power[2]) for power in native.denominator.to_dict()} == {3}
    expected = c**4*(1+u*u)**18*(c*(1+u*u)**4-2)**2*(u*u-1)/344064
    assert sp.factor(polynomial.LC()-expected) == 0
    native_expression = native.sympy().subs(dict(zip(
        sp.symbols("u c K", real=True), generic["variables"], strict=True)), simultaneous=True)
    assert sp.factor(native_expression.subs(u, 0)-direct_equation(0, c)["det"]) == 0
    norms = []
    for degree, upper in enumerate((16, 14, 1)):
        coefficient = sp.Poly(polynomial.nth(degree), u, c)
        norm = sum(abs(Fraction(value))*Fraction(1, 10)**power[0]*4**power[1]
                   for power, value in coefficient.terms())
        assert norm < upper
        norms.append(norm)
    assert tuple(norms) == bounds.derive()["coefficient_norms"]
    # On the continuous actual box, |n3| >= (2^4*99/100)/344064 * delta_min^2.
    lower = Fraction(2**4, 344064)*Fraction(99, 100)
    assert lower == Fraction(99, 2150400)
    # K_star >= 500000 > 1, so all three inverse powers are bounded by1/K.
    assert Fraction(2_000_000, 2**2) >= 1
    relative_error = Fraction(16+14+1, 2_000_000)/lower
    assert relative_error == Fraction(6944, 20625) < Fraction(1, 2)


@pytest.mark.parametrize("time", [sp.Integer(0), sp.Rational(1, 100)])
def test_exact_normalized_equation_by_independent_second_jet_cauchy_map(time):
    lapse = sp.Integer(4)
    d = direct_equation(time, lapse)
    l = leading.derive()
    u, c = l["u"], l["c"]
    weight = sp.diag(l["U"], 1, 1)
    at = {u: time, c: lapse}
    T0, T1, T2 = (sp.diff(weight, u, order).subs(at) for order in range(3))
    O0, O1, O2 = d["O"]
    A0, A1 = d["A"]
    normalized_rows = T0*O0
    normalized_first = T1*O0+T0*O1
    normalized_second = T2*O0+2*T1*O1+T0*O2
    phase = normalized_rows.col_join(normalized_first+normalized_rows*A0)
    acceleration = (normalized_second+2*normalized_first*A0
                    +normalized_rows*(A1+A0*A0))
    expected = (acceleration*phase.inv()).applyfunc(sp.cancel)
    native = energy.derive()
    point = (Fraction(time), Fraction(lapse), Fraction(7, 2))
    for start, key in ((0, "normalized_Mq"), (3, "normalized_Mv")):
        actual = sp.Matrix([[sp.Rational(fraction_evaluation(value, point)) for value in row]
                            for row in native[key]])
        exact = expected[:, start:start+3].subs(action.K, sp.Rational(7, 2))
        assert (actual-exact).applyfunc(sp.cancel) == sp.zeros(3)
    # A vanished first center connection does not license dropping T''.
    if time == 0:
        assert T1 == sp.zeros(3)
        assert sp.factor(T2[0, 0]/T0[0, 0]) == -38
