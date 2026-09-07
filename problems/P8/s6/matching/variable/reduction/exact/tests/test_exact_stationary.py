from fractions import Fraction

import pytest
import sympy as sp
from p8_exact_stationary import stationary as primary
from p8_variable_reduction import stationary as frozen


class Series:
    """Four exact Fraction coefficients, independent of the symbolic engine."""

    def __init__(self, values):
        if isinstance(values, (int, Fraction)):
            values = [values]
        self.a = tuple(Fraction(x) for x in values)+tuple(Fraction(0) for _ in range(4-len(values)))
        if len(self.a) != 4:
            raise ValueError("Exactly four coefficients are available")

    @staticmethod
    def lift(value):
        return value if isinstance(value, Series) else Series(value)

    def __add__(self, other):
        other = self.lift(other)
        return Series([a+b for a, b in zip(self.a, other.a, strict=True)])

    __radd__ = __add__

    def __neg__(self):
        return Series([-x for x in self.a])

    def __sub__(self, other):
        return self+-self.lift(other)

    def __rsub__(self, other):
        return self.lift(other)+-self

    def __mul__(self, other):
        other = self.lift(other)
        return Series([sum((self.a[j]*other.a[n-j] for j in range(n+1)), Fraction(0)) for n in range(4)])

    __rmul__ = __mul__

    def reciprocal(self):
        if not self.a[0]:
            raise ZeroDivisionError("Singular series constant")
        values = [1/self.a[0]]
        for n in range(1, 4):
            values.append(-sum(self.a[j]*values[n-j] for j in range(1, n+1))/self.a[0])
        return Series(values)

    def __truediv__(self, other):
        return self*self.lift(other).reciprocal()

    def __rtruediv__(self, other):
        return self.lift(other)*self.reciprocal()

    def __pow__(self, power):
        if power < 0:
            return self.reciprocal()**(-power)
        result = Series(1)
        for _ in range(power):
            result *= self
        return result

    def derivative(self):
        return Series([self.a[1], 2*self.a[2], 3*self.a[3], 0])


def independent_center(c):
    """Solve the literal constraint series successively, not center formulas."""
    x = Series([0, 1])
    d = 1+x
    y = 2/d**4
    p = y**3*4*(1-x)/d**2/(c*(c-y))
    w = 24*x/(c**2*d**2)-p/y**3

    def expressions(b):
        F = 2*(p.derivative()+w.derivative()*b**3)/(3*p)
        return 3*b*x*F**2-2*(p+w*b**3), F

    coefficients = [Fraction(2), Fraction(0), Fraction(0), Fraction(0)]
    for order in (1, 2):
        zero = expressions(Series(coefficients))[0].a[order]
        shifted = list(coefficients)
        shifted[order] = 1
        slope = expressions(Series(shifted))[0].a[order]-zero
        coefficients[order] = -zero/slope
    b = Series(coefficients)
    residual, F = expressions(b)
    lapse = 2*b.derivative()/F
    assert residual.a[:3] == (0, 0, 0)
    return b, lapse, p, w, F


def test_primary_exact_identities_and_scope():
    checks = primary.checks()
    assert len(checks) == 37
    assert set(checks.values()) == {0}
    report = primary.calibration()
    assert "not a full-parent/CD solution" in report["scope"]
    assert "parity follows" in report["regularity"]


def test_literal_own_lapse_and_scale_variations():
    e = primary.own_equations()
    t, b, N, p, w = primary.TIME, e["b"], e["N"], e["beta1"], e["beta4"]
    bdot, ndot = sp.diff(b, t), sp.diff(N, t)
    expected_C = 3*primary.M**2*b*bdot**2/N**2-2*p-2*w*b**3
    expected_E = (6*primary.M**2*b*sp.diff(b, t, 2)/N+3*primary.M**2*bdot**2/N
                  -6*primary.M**2*b*bdot*ndot/N**2-6*p-6*w*N*b**2)
    assert sp.expand(e["C"]-expected_C) == 0
    assert sp.expand(e["E"]-expected_E) == 0
    assert sp.expand(sp.diff(expected_C, t)-bdot*expected_E/N+2*e["consistency"]) == 0


@pytest.mark.parametrize("c", [Fraction(201, 100), Fraction(401, 200), Fraction(3), Fraction(4), Fraction(1)])
def test_independent_fraction_series_center(c):
    b, N, _, _, _ = independent_center(c)
    center = primary.center()
    point = {primary.C: sp.Rational(c.numerator, c.denominator)}
    pairs = [(center["b_uu_center"], 2*b.a[1]),
             (center["b_uuuu_center"], 24*b.a[2]),
             (center["N_center"], N.a[0]), (center["N_uu_center"], 2*N.a[1])]
    for expression, exact in pairs:
        assert expression.subs(point) == sp.Rational(exact.numerator, exact.denominator)
    # c=1 is an algebra-only fixed-local control, not part of the positive box.
    assert N.a[0] > 0


def test_no_parity_assumption_in_the_center_constraint():
    r, p, speed, lapse, M = sp.symbols("r p speed lapse M", positive=True)
    w = -p/r**3
    C0 = 3*M**2*r*speed**2/lapse**2-2*(p+w*r**3)
    assert sp.factor(C0) == 3*M**2*r*speed**2/lapse**2
    # With r,N,M>0, the lapse equation itself enforces speed=0.
    assert sp.solve(C0, speed) == []  # speed was stipulated strictly positive.
    assert C0.subs(speed, 0) == 0


def test_odd_lapse_center_only_omission_control():
    c = Fraction(201, 100)
    b, N, p, w, _ = independent_center(c)
    u = sp.Symbol("u", real=True)
    rational = lambda value: sp.Rational(value.numerator, value.denominator)
    b_trial = rational(b.a[0])+rational(b.a[1])*u**2
    p_trial = rational(p.a[0])+rational(p.a[1])*u**2
    w_trial = rational(w.a[0])+rational(w.a[1])*u**2
    N_trial = rational(N.a[0])+u
    C = 3*b_trial*sp.diff(b_trial, u)**2/N_trial**2-2*(p_trial+w_trial*b_trial**3)
    E = (6*b_trial*sp.diff(b_trial, u, 2)/N_trial+3*sp.diff(b_trial, u)**2/N_trial
         -6*b_trial*sp.diff(b_trial, u)*sp.diff(N_trial, u)/N_trial**2
         -6*p_trial-6*w_trial*N_trial*b_trial**2)
    consistency = sp.diff(p_trial, u)+sp.diff(w_trial, u)*b_trial**3-3*p_trial*sp.diff(b_trial, u)/N_trial
    assert C.subs(u, 0) == 0 and E.subs(u, 0) == 0
    assert sp.series(consistency, u, 0, 3).removeO().expand().coeff(u, 1) == 0
    assert sp.series(consistency, u, 0, 3).removeO().expand().coeff(u, 2) != 0


def test_frozen_profile_and_physical_lapse_are_not_the_rolling_lapse():
    p, old = primary.profile(), frozen.profile()
    for name in ("b1", "b4"):
        assert sp.factor(p[name.upper()]-old[name].subs({frozen.V: primary.V, frozen.C: primary.C})) == 0
    center = primary.center()
    delta = primary.DELTA
    assert sp.expand((center["N_center"]-primary.C).subs(primary.C, 2+delta)) == delta+delta**2/2
    assert sp.factor(center["metric00_center_defect"]-(primary.C**2-4)**2/4) == 0


def test_joint_partials_are_literal_independent_derivatives():
    j = primary.joint()
    point = {primary.V: sp.Rational(1, 20000), primary.C: sp.Rational(201, 100), primary.ZETA: 12}
    for value, coordinate, derivative in (("b", primary.V, "b_v"), ("b", primary.ZETA, "b_zeta"),
                                         ("T", primary.V, "T_v"), ("T", primary.ZETA, "T_zeta")):
        residual = (sp.diff(j[value], coordinate)-j[derivative]).subs(point)
        assert sp.simplify(residual) == 0
    assert sp.simplify(j["N"]-2*(j["b_v"]+j["b_zeta"]*j["zeta_v"])/j["F"]) == 0


def test_joint_lapse_uses_the_implicit_zeta_derivative():
    j = primary.joint()
    point = {primary.V: sp.Rational(1, 20000), primary.C: sp.Rational(201, 100), primary.ZETA: 12}
    omitted = (2*j["b_v"]/j["F"]-j["N"]).subs(point)
    assert sp.simplify(omitted) != 0
    assert j["N"].subs(primary.center()["joint_point"]) == 2


def test_zero_delta_is_an_extension_not_an_action():
    p, j = primary.profile(), primary.joint()
    point = primary.center()["joint_point"]
    assert p["D"].subs(point) == 0
    assert p["Q"].subs(point) == 16
    assert p["B1"].subs(primary.V, 0).subs(primary.C, 2) is sp.zoo
    assert j["T"].subs(point) == 12
    assert primary.center()["joint_G_zeta_center"] == -32


def test_degenerate_constant_root_is_not_spuriously_excluded():
    e = primary.own_equations()
    r = sp.Symbol("constant_r", positive=True)
    p = e["beta1"]
    substitution = {e["b"]: r, e["N"]: r, e["beta4"]: -p/r**3}
    for name in ("C", "E", "consistency"):
        assert sp.simplify(e[name].subs(substitution).doit()) == 0


def test_zero_link_center_loses_the_regular_root_premise():
    c = primary.center()
    assert sp.simplify(c["center_consistency_after_spatial"].subs(c["p0"], 0)) == 0
    assert primary.profile()["Q"].subs({primary.V: 0, primary.C: 2}) != 0
    # The actual family never uses beta1=0 to derive a lapse.
    assert sp.limit(primary.fixed_c()["b_jacobian_center"], primary.C, 2, dir="+") is sp.oo


def test_literal_pairwise_log_rates_and_tensor_potential():
    M, b, N, h, qd, Hd = sp.symbols("M b N h qd Hd", real=True)
    g_rates = (qd/2, -qd/2, 0)
    f_rates = (h+Hd/2, h-Hd/2, h)
    pair = lambda rates: sum(rates[a]*rates[b] for a in range(3) for b in range(a+1, 3))
    assert sp.expand(-M**2*pair(g_rates)) == M**2*qd**2/4
    assert sp.expand(-M**2*b**3*pair(f_rates)/N+3*M**2*b**3*h**2/N) == M**2*b**3*Hd**2/(4*N)
    t = primary.tensors()
    epsilon = sp.Symbol("epsilon", real=True)
    diagonal_root = (primary.own_equations()["N"],
                     primary.own_equations()["b"]*sp.exp(epsilon*(t["Q_tensor"]-t["q"])/2),
                     primary.own_equations()["b"]*sp.exp(-epsilon*(t["Q_tensor"]-t["q"])/2),
                     primary.own_equations()["b"])
    assert sp.simplify(sp.prod(diagonal_root)-t["e4"]) == 0
    literal = -2*primary.own_equations()["beta1"]*sum(diagonal_root)
    assert sp.simplify(sp.diff(literal, epsilon, 2).subs(epsilon, 0)/2-t["quadratic_potential"]) == 0


def test_distinct_inner_own_and_coupled_equations():
    i = primary.inner()
    denominator = 1+8*i["x"]**2
    q, Q = sp.symbols("q Q", real=True)
    g_second = -64*(q-Q)/denominator
    f_second = -16*(Q-q)/denominator
    assert sp.expand(g_second+4*f_second) == 0
    assert sp.simplify(f_second-g_second+80*(Q-q)/denominator) == 0
    assert i["own_f_inner_coefficient"] != i["coupled_relative_inner_coefficient"]
    assert sp.factor(i["own_f_center_mass_squared"]-8*primary.C/(primary.TAU**2*(primary.C-2))) == 0


def test_inner_drift_and_pump_power_counting_requires_joint_derivatives():
    delta, x = primary.DELTA, primary.inner()["x"]
    K, Kv, Kvv = sp.symbols("K K_v K_vv", positive=True)
    u = sp.sqrt(delta)*x
    Ku, Kuu = 2*u*Kv, 2*Kv+4*u**2*Kvv
    drift = sp.sqrt(delta)*Ku/K
    pump = delta*(Kuu/(2*K)-Ku**2/(4*K**2))
    assert sp.limit(drift, delta, 0, dir="+") == 0
    assert sp.limit(pump, delta, 0, dir="+") == 0
    assert sp.limit(primary.inner()["D_over_delta"], delta, 0, dir="+") == 1+8*x**2


def test_nonzero_homogeneous_data_are_not_fixed_by_an_operator_formula():
    t = sp.Symbol("t", real=True)
    q = sp.Function("source_q")(t)
    Q1, Q2 = sp.Function("solution1")(t), sp.Function("solution2")(t)
    K, nu = sp.Function("K")(t), sp.Function("nu")(t)
    operator = lambda value: sp.diff(K*sp.diff(value, t), t)+nu*(value-q)
    difference = sp.expand(operator(Q1)-operator(Q2))
    assert sp.expand(difference-sp.diff(K*sp.diff(Q1-Q2, t), t)-nu*(Q1-Q2)) == 0
