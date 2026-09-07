"""Independent Fraction coframes/Taylor jets and ordinary time variations."""

from fractions import Fraction as F

import pytest
import sympy as sp
from p8_sixth_reduction import tensor as t
from p8_variable_reduction import stationary as fourth


class Jet:
    """Four-coefficient Fraction power series, independent of the primary jets."""

    def __init__(self, values=0):
        if isinstance(values, Jet):
            self.a = values.a
        elif isinstance(values, (tuple, list)):
            self.a = tuple(map(F, values))+tuple(F(0) for _ in range(4-len(values)))
        else:
            self.a = (F(values), F(0), F(0), F(0))

    def __add__(self, other):
        return Jet(tuple(x+y for x, y in zip(self.a, Jet(other).a, strict=True)))

    __radd__ = __add__

    def __neg__(self):
        return Jet(tuple(-x for x in self.a))

    def __sub__(self, other):
        return self+-Jet(other)

    def __rsub__(self, other):
        return Jet(other)+-self

    def __mul__(self, other):
        other = Jet(other)
        return Jet(tuple(sum(self.a[k]*other.a[n-k] for k in range(n+1)) for n in range(4)))

    __rmul__ = __mul__

    def inverse(self):
        result = [1/self.a[0]]
        for n in range(1, 4):
            result.append(-sum(self.a[k]*result[n-k] for k in range(1, n+1))/self.a[0])
        return Jet(result)

    def __truediv__(self, other):
        return self*Jet(other).inverse()

    def __rtruediv__(self, other):
        return Jet(other)*self.inverse()

    def power(self, exponent, root0=None):
        exponent = F(exponent)
        if root0 is None:
            assert exponent.denominator == 1
            root0 = self.a[0]**exponent.numerator
        z, power, choose, result = self/self.a[0]-1, Jet(1), F(1), Jet(1)
        for n in range(1, 4):
            power *= z
            choose *= (exponent-n+1)/n
            result += choose*power
        return root0*result

    def derivative(self):
        return Jet(tuple((n+1)*self.a[n+1] for n in range(3)))


def _fraction_coframes(r2, A, Ap, ell, ellp, ellpp, qp, qpp, qppp):
    rates = (qp/2, -qp/2, F(0))
    ratesp = (qpp/2, -qpp/2, F(0))
    ratespp = (qppp/2, -qppp/2, F(0))
    # Direct conformal diagonal Ricci components in the physical-g mixed frame.
    ricci = [-3*ellp-qp**2/2]
    riccip = [-3*ellpp-qp*qpp]
    for h, hp, hpp in zip(rates, ratesp, ratespp, strict=True):
        ricci.append(-ellp-hp-2*ell**2-2*ell*h)
        riccip.append(-ellpp-hpp-4*ell*ellp-2*ellp*h-2*ell*hp)
    scalar, scalarp = sum(ricci), sum(riccip)
    P = [value-scalar/6 for value in ricci]
    Pp = [value-scalarp/6 for value in riccip]
    h = [A*value for value in P]
    hp = [Ap*value+A*derivative for value, derivative in zip(P, Pp, strict=True)]
    eps = Jet((0, 1))
    roots = [(1+eps*value).power(F(1, 2), root0=F(1)) for value in h]
    volume_over_lapse = roots[1]*roots[2]*roots[3]/roots[0]
    rates_f = [ell+rates[i]+eps*hp[i+1]/(2*(1+eps*h[i+1])) for i in range(3)]
    EH = -r2*volume_over_lapse*sum(rates_f[i]*rates_f[j] for i in range(3) for j in range(i+1, 3))
    boundary = r2*volume_over_lapse*sum(rates_f)
    potential = -2*r2/A*(sum(roots)-roots[0]*roots[1]*roots[2]*roots[3])
    return {"h": h, "hp": hp, "EH2": EH.a[2], "boundary2": boundary.a[2], "potential3": potential.a[3]}


def _fraction_center(c):
    c, v = F(c), Jet((0, 1))
    d = 1+v
    y = 2*d.power(-4)
    h2, hp = 16*v/d.power(2), 4*(1-v)/d.power(2)
    b1 = y.power(3)*hp/(c*(c-y))
    b4 = 3*h2/(2*c*c)-b1/y.power(3)
    r3 = -b1/b4
    r2 = r3.power(F(2, 3), root0=4)
    kappa = -1/(4*b4)
    A = 4*kappa/r2
    ell = F(2, 3)*r3.derivative()/r3
    A0, A2, A4 = A.a[0], 2*A.a[1], 24*A.a[2]
    lp, lppp = ell.a[0], 6*ell.a[1]
    return {"A0": A0, "A2": A2, "A4": A4,
            "E6": -A0**2/2,
            "E4": -(8*A0**2*lp+7*A0*A2)/2,
            "E2": -(A2**2+22*A0*A2*lp+26*A0**2*lp**2+2*A0**2*lppp+A0*A4)/2,
            "fourth_E2": 2*kappa.a[1]}


def test_primary_identities_and_continuous_margin():
    assert len(t.checks()) == 30
    assert set(t.checks().values()) == {0}
    assert t.calibration()["sixth_small_delta_low_absolute_floor_over_M2"] == sp.Rational(793, 400)


@pytest.mark.parametrize("values", [
    (F(4), F(1, 2), F(2, 3), F(-1, 3), F(4, 5), F(-2, 7), F(3, 5), F(-5, 4), F(7, 3)),
    (F(9, 4), F(3, 7), F(-1, 2), F(2, 5), F(-3, 4), F(5, 6), F(-2, 3), F(7, 5), F(-4, 7)),
    (F(4), F(1, 2), F(0), F(0), F(-6), F(0), F(1, 3), F(2, 5), F(4, 7)),
])
def test_literal_fraction_metric_EH_and_potential_expansion(values):
    independent, primary = _fraction_coframes(*values), t.derive()
    symbols = (t.R2, t.INV[0], t.INV[1], t.LOG[0], t.LOG[1], t.LOG[2], t.Q[1], t.Q[2], t.Q[3])
    point = dict(zip(symbols, values, strict=True))
    assert primary["EH_second_coefficient"].subs(point) == independent["EH2"]
    assert primary["potential_third_coefficient"].subs(point) == independent["potential3"]
    assert primary["covariant_EH_boundary_second"].subs(point) == independent["boundary2"]
    for actual, expected in zip(primary["h_mixed_diagonal"], independent["h"], strict=True):
        assert actual.subs(point) == expected
    for actual, expected in zip(primary["h_time_derivative"], independent["hp"], strict=True):
        assert actual.subs(point) == expected


def test_literal_covariant_diagonal_EH_boundary_sign():
    time = sp.Symbol("time", real=True)
    lapse = sp.Function("N", positive=True)(time)
    scales = [sp.Function(f"a{i}", positive=True)(time) for i in range(3)]
    metric = sp.diag(lapse**2, *(-a**2 for a in scales))
    inverse = metric.inv()

    def partial(expression, coordinate):
        return sp.diff(expression, time) if coordinate == 0 else sp.S.Zero

    gamma = [[[sp.simplify(sum(inverse[a, d]*(partial(metric[d, c], b)+partial(metric[d, b], c)
                                            -partial(metric[b, c], d))/2 for d in range(4)))
               for c in range(4)] for b in range(4)] for a in range(4)]
    ricci = sp.Matrix(4, 4, lambda b, c: sp.simplify(sum(
        partial(gamma[a][c][b], a)-partial(gamma[a][a][b], c)
        +sum(gamma[a][a][d]*gamma[d][c][b]-gamma[a][c][d]*gamma[d][a][b] for d in range(4))
        for a in range(4))))
    scalar = sp.simplify(sp.trace(inverse*ricci))
    volume = sp.prod(scales)
    rates = [sp.diff(a, time)/a for a in scales]
    covariant = -lapse*volume*scalar/2
    ADM = -volume/lapse*sum(rates[i]*rates[j] for i in range(3) for j in range(i+1, 3))
    boundary = volume/lapse*sum(rates)
    assert sp.factor(covariant-ADM-sp.diff(boundary, time)) == 0
    assert sp.factor(covariant-ADM+sp.diff(boundary, time)) != 0


@pytest.mark.parametrize("c", [F(201, 100), F(5, 2), F(3), F(4)])
def test_fraction_actual_clock_profile_and_Euler_center(c):
    independent, primary = _fraction_center(c), t.center()
    point = {t.C: c}
    for order in (2, 4, 6):
        assert primary["sixth_E_dimensionless"][order].subs(point) == independent[f"E{order}"]
    assert primary["fourth_E_dimensionless"][2].subs(point) == independent["fourth_E2"]
    jets = t.center_jets()["substitution"]
    for order, name in ((0, "A0"), (2, "A2"), (4, "A4")):
        assert jets[t.INV[order]].subs(point) == independent[name]


def test_original_stationary_profile_bridge_and_constant_ratio_high_symbol():
    old, new = fourth.center_jets(), t.center_jets()
    point = new["substitution"]
    assert sp.factor(point[t.INV[0]]-2*old["own_f_inverse_parameter_bar_u0_u2"][0]) == 0
    assert sp.factor(point[t.INV[2]]-2*old["own_f_inverse_parameter_bar_u0_u2"][1]) == 0
    assert sp.factor(point[t.LOG[1]]-old["r_squared_u0_u2_u4"][1]/8) == 0
    calibration = fourth.constant_ratio_calibration()
    high = t.calibration()["sixth_highest_symbol"].subs({t.R2: fourth.R**2, t.INV[0]: fourth.M**2*fourth.R/fourth.BETA1})
    assert sp.factor(high-calibration["six_derivative_flat_curvature_coefficient"]) == 0


def test_arbitrary_composed_clock_chain_and_boundary():
    time = sp.Symbol("time", real=True)
    theta = time+time**2
    r, A = 2+theta**2, sp.Rational(1, 2)+theta+theta**2
    q = time**2-time**3+time**4
    ell = sp.diff(r, time)/r
    point = {t.R2: r**2}
    for symbols, expression in ((t.INV, A), (t.LOG, ell), (t.Q, q)):
        point.update({symbol: sp.diff(expression, time, n) for n, symbol in enumerate(symbols)})
    d = t.derive()
    # Ordinary differentiation of functions, independent of the jet derivation.
    for name in ("B1", "B2", "B3", "time_IBP_boundary"):
        expression = d[name]
        expected = sp.diff(expression.subs(point, simultaneous=True), time)
        actual = t.total_derivative(expression).subs(point, simultaneous=True)
        assert sp.factor((actual-expected).subs(time, sp.Rational(1, 7))) == 0
    residual = (d["raw_quadratic"]-d["normal_quadratic"]).subs(point, simultaneous=True)
    residual -= sp.diff(d["time_IBP_boundary"].subs(point, simultaneous=True), time)
    assert sp.factor(residual.subs(time, sp.Rational(-1, 11))) == 0


def test_center_before_time_derivatives_is_a_failing_control():
    c = t.center()
    incorrectly_frozen_low = -2*c["B_dimensionless"][1]
    assert incorrectly_frozen_low.subs(t.C, 2) == 0
    assert c["sixth_E_dimensionless"][2].subs(t.C, 2) == -2
    assert incorrectly_frozen_low.subs(t.C, 4) == 60
    assert c["sixth_E_dimensionless"][2].subs(t.C, 4) == -104


def test_fourth_sixth_cancellation_must_not_be_omitted():
    c = t.center()
    assert c["fourth_E_dimensionless"][2].subs(t.C, 2) == 2
    assert c["sixth_E_dimensionless"][2].subs(t.C, 2) == -2
    assert c["combined_E_dimensionless"][2].subs(t.C, 2) == 0
    assert c["combined_low_delta_polynomial"].all_coeffs() == [
        -sp.Rational(349, 128), -sp.Rational(143, 32), -sp.Rational(173, 32), sp.Rational(7, 2), 0]
    assert c["sixth_E_dimensionless"][4].subs(t.C, 4) == -22
    assert c["sixth_E_dimensionless"][6].subs(t.C, 4) == -sp.Rational(1, 8)


def test_literal_scale_restoration_and_normalized_floor():
    c = t.center()
    for n in (2, 4, 6):
        assert sp.factor(c["sixth_E_physical"][n]/(t.M**2*t.TAU**(n-2))-c["sixth_E_dimensionless"][n]) == 0
    assert t.calibration()["sixth_low_over_Mstar_squared_limit"] == -sp.Rational(2, 5)
    assert c["sixth_low_absolute_lower_bound"] > sp.Rational(19, 10)


def test_finite_jet_order_is_not_silently_truncated():
    for value in (t.INV[-1], t.LOG[-1], t.Q[-1]):
        with pytest.raises(ValueError):
            t.total_derivative(value)
