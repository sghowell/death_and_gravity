"""Own-f checks, with a separate Fraction Taylor/determinant route."""

from fractions import Fraction as F
from itertools import permutations

import pytest
import sympy as sp
from p8_bimetric import matching as old_matching
from p8_bimetric.background import MF2, MG2, NU
from p8_variable_beta import background as old_background
from p8_variable_reduction import stationary as s


class Jet:
    """Independent four-coefficient Fraction ring; no symbolic algebra."""

    def __init__(self, values=0):
        if isinstance(values, Jet):
            self.a = values.a
        elif isinstance(values, (tuple, list)):
            self.a = tuple(F(value) for value in values)+tuple(F(0) for _ in range(4-len(values)))
        else:
            self.a = (F(values), F(0), F(0), F(0))

    def __add__(self, other):
        other = Jet(other)
        return Jet(tuple(a+b for a, b in zip(self.a, other.a, strict=True)))

    __radd__ = __add__

    def __neg__(self):
        return Jet(tuple(-a for a in self.a))

    def __sub__(self, other):
        return self+-Jet(other)

    def __rsub__(self, other):
        return Jet(other)+-self

    def __mul__(self, other):
        other = Jet(other)
        return Jet(tuple(sum(self.a[j]*other.a[n-j] for j in range(n+1)) for n in range(4)))

    __rmul__ = __mul__

    def inverse(self):
        result = [1/self.a[0]]
        for n in range(1, 4):
            result.append(-sum(self.a[j]*result[n-j] for j in range(1, n+1))/self.a[0])
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
        z = self/self.a[0]-1
        result, power, choose = Jet(1), Jet(1), F(1)
        for n in range(1, 4):
            power = power*z
            choose *= (exponent-n+1)/n
            result += choose*power
        return root0*result

    def derivative(self):
        return Jet(tuple((n+1)*self.a[n+1] for n in range(3)))


def _determinant(matrix):
    answer = Jet(0)
    for p in permutations(range(4)):
        inversions = sum(p[i] > p[j] for i in range(4) for j in range(i+1, 4))
        product = Jet((-1)**inversions)
        for i in range(4):
            product *= matrix[i][p[i]]
        answer += product
    return answer


def _multiply(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def _fraction_profile(c):
    v, c = Jet((0, 1)), F(c)
    d = 1+v
    y = 2*d.power(-4)
    h2, hp = 16*v/d.power(2), 4*(1-v)/d.power(2)
    b1 = y.power(3)*hp/(c*(c-y))
    b4 = 3*h2/(2*c*c)-b1/y.power(3)
    r3 = -b1/b4
    r2 = r3.power(F(2, 3), root0=F(4))
    ell = F(2, 3)*r3.derivative()/r3
    P00 = -2*(hp+ell+2*v*ell.derivative())-h2+v*ell.power(2)
    kappa = -1/(4*b4)
    f2 = 4*kappa*P00
    space = 4*kappa*v*(4/d+ell).power(2)
    k = 2*(y.power(3)/c-1)*hp-1/(100*d.power(12))
    return {"r2": r2, "f2": f2, "space": space, "kappa": kappa,
            "remainder": c*c-r2-f2, "k": k}


def test_all_primary_identities_and_margins():
    assert len(s.checks()) == 56
    assert set(s.checks().values()) == {0}
    assert all(value > 0 for value in s.domain_margins().values())
    assert s.calibration()["formal_kappa"] == s.M**4*s.R**3/(4*s.BETA1)


@pytest.mark.parametrize("b", [F(3, 2), F(-5, 3)])
def test_literal_full_fraction_determinant_potential(b):
    lower = [[F(2), F(1, 2), F(-1, 3), F(2, 5)],
             [F(1, 2), F(-1), F(3, 4), F(1, 7)],
             [F(-1, 3), F(3, 4), F(4), F(-1, 6)],
             [F(2, 5), F(1, 7), F(-1, 6), F(-2)]]
    signs = (1, -1, -1, -1)
    H = [[signs[i]*lower[i][j] for j in range(4)] for i in range(4)]
    H2, H3 = _multiply(H, H), _multiply(_multiply(H, H), H)
    root_trace = Jet((4, sum(H[i][i] for i in range(4))/2,
                      -sum(H2[i][i] for i in range(4))/8,
                      sum(H3[i][i] for i in range(4))/16))
    determinant = _determinant([[Jet((int(i == j), H[i][j])) for j in range(4)] for i in range(4)])
    volume = determinant.power(F(1, 2), root0=F(1))
    literal = -2*b*(root_trace-volume)
    data = s.potential_hessian()
    point = {s.R: 1, s.BETA1: sp.Rational(b.numerator, b.denominator)}
    point.update({data["h"][i, j]: sp.Rational(lower[i][j].numerator, lower[i][j].denominator)
                  for i in range(4) for j in range(i, 4)})
    assert literal.a[1] == 0
    assert sp.factor(data["quadratic"].subs(point)) == literal.a[2]
    assert sp.factor(data["cubic"].subs(point)) == literal.a[3]


def test_full_trace_inverse_not_massless_one_half():
    p = s.potential_hessian()
    E = s.ETA
    true = E-s.ETA*sp.trace(s.ETA*E)/3
    wrong = E-s.ETA*sp.trace(s.ETA*E)/2
    assert true-s.ETA*sp.trace(s.ETA*true) == E
    assert wrong-s.ETA*sp.trace(s.ETA*wrong) != E
    assert sp.factor(p["hessian"].det()) != 0
    assert sp.factor(p["hessian"].det().subs(s.BETA1, 0)) == 0


def test_covariant_EH_sign_by_literal_two_lapse_curvature():
    time = sp.Symbol("time", real=True)
    a, r = (sp.Function(name)(time) for name in ("a", "r"))
    h = sp.diff(a, time)/a
    Hf = sp.diff(a*r, time)/(a*r*r)
    Rg = -6*(sp.diff(h, time)+2*h**2)
    Rf = -6*(sp.diff(Hf, time)/r+2*Hf**2)
    literal = -s.M**2*a**3*r**4*Rf/2
    retained = -s.M**2*a**3*r**2*Rg/2-3*s.M**2*a**3*sp.diff(r, time)**2
    boundary = 3*s.M**2*sp.diff(a**3*r*sp.diff(r, time), time)
    assert sp.factor(literal-retained-boundary) == 0
    assert sp.factor(literal-retained+boundary) != 0


def test_constant_r_recovers_pinned_source_preserving_S63():
    old, new = old_matching.tt_schur(), s.constant_ratio_calibration()
    at = {MG2: new["G"], MF2: new["F"], NU: new["nu"], old_matching.D: new["D"]}
    for old_key, new_key in (("physical_g_kernel", "kernel"), ("kernel_remainder", "kernel_remainder"), ("c_C", "c_C")):
        assert sp.factor(old[old_key].subs(at, simultaneous=True)-new[new_key]) == 0
    assert new["nu"] == 2*s.R*s.BETA1
    # Literal positive-canonical stress, not a sign assigned to an abstract j.
    amplitude = sp.Symbol("amplitude", real=True)
    gradient = sp.Matrix(sp.symbols("chi0:4", real=True))
    metric = sp.diag(1, -1-amplitude/sp.sqrt(2), -1+amplitude/sp.sqrt(2), -1)
    matter = sp.sqrt(-metric.det())*(gradient.T*metric.inv()*gradient)[0]/2
    Pi = (gradient[1]**2-gradient[2]**2)/sp.sqrt(2)
    assert sp.simplify(sp.diff(matter, amplitude).subs(amplitude, 0)-Pi/2) == 0
    assert "j=+2 Pi_TT" in new["source_scope"]


def test_constant_clock_gradient_conformal_boundary_is_not_zero():
    a = s.conformal_action()
    at = {a["gradient"][i]: int(i == 0) for i in range(4)}
    for key in ("ricci", "hessian"):
        at.update(dict.fromkeys(a[key], 0))
    assert sp.expand(a["Qhat"].subs(at)) == -12*a["s"]**2*a["t"]
    assert sp.expand(a["div_V"].subs(at)) == 3*a["s"]**2*a["t"]
    assert sp.expand(a["compact_ibp"].subs(at)) == 4*a["kappa_phi"]*a["s"]**3


@pytest.mark.parametrize("c", [F(201, 100), F(5, 2), F(3), F(4)])
def test_independent_fraction_profile_jets(c):
    p, data = _fraction_profile(c), s.center_jets()
    point = {s.C: sp.Rational(c.numerator, c.denominator)}
    pairs = (("r2", "r_squared_u0_u2_u4"), ("f2", "f2_00_u0_u2_u4"),
             ("space", "f2_space_over_a2_u0_u2_u4"), ("kappa", "kappa_bar_u0_u2_u4"),
             ("remainder", "metric00_remainder_u0_u2_u4"))
    for independent, primary in pairs:
        for n, factor in enumerate((1, 2, 24)):
            assert data[primary][n].subs(point) == factor*p[independent].a[n]
    rem, k, f2 = p["remainder"].a, p["k"].a, p["f2"].a
    assert data["f2_00_varphi2"].subs(point) == 2*f2[1]/k[0]
    assert data["remainder00_varphi4"].subs(point) == 24*rem[2]/k[0]**2-8*rem[1]*k[1]/k[0]**3


def test_profile_matches_unchanged_variable_coefficients_continuously():
    p, old = s.profile(), old_background.derive()
    at = {s.C: old_background.C, s.V: old_background.U**2}
    for name in ("y", "b1", "b4", "kbar"):
        assert sp.factor(p[name].subs(at)-old[name]) == 0
    assert sp.factor(p["r_cubed"]+p["b1"]/p["b4"]) == 0
    assert sp.factor(p["kappa_bar"]-p["r_cubed"]/(4*p["b1"])) == 0


def test_continuous_fourth_jet_floor_and_moving_clock_term():
    j = s.center_jets()
    assert j["remainder_u4_minus_10752_cleared"].all_coeffs() == [26832, 87552, 132288, 22272, 0]
    assert j["f2_u2_above_60_for_delta_le_1over100_margin"] > 5
    assert j["f2_00_varphi2_limit"] == sp.Rational(6400, 2399)
    assert j["remainder00_varphi4_limit"] == sp.Rational(107520000, 5755201)
    naive = j["metric00_remainder_u0_u2_u4"][2]/j["kbar_center"]**2
    assert sp.factor((j["remainder00_varphi4"]-naive).subs(s.C, 3)) > 0


def test_small_center_value_does_not_erase_derivative_defect():
    j = s.center_jets()
    assert j["metric00_remainder_u0_u2_u4"][0] == (s.C-2)**2
    assert j["metric00_remainder_u0_u2_u4"][1].subs(s.C, 2) == 0
    assert j["metric00_remainder_u0_u2_u4"][2].subs(s.C, 2) == 10752
    assert j["f2_00_u0_u2_u4"][0].subs(s.C, 2) == 0
    assert j["f2_00_u0_u2_u4"][1].subs(s.C, 2) == 64
    assert j["own_f_inverse_parameter_bar_u0_u2"][0] == s.C*(s.C-2)/32
    assert j["own_f_inverse_parameter_bar_u0_u2"][1].subs(s.C, 2) == 1


def test_stationary_order_count_requires_next_operator():
    # Generic scalar algebra of the stationary expansion: q4 drops out of
    # S6 because H0*q2+J2=0, but the q2 kinetic/cubic terms do not.
    e, H0, J2, K2, T0, q4 = sp.symbols("e H0 J2 K2 T0 q4", nonzero=True)
    q2 = -J2/H0
    q = e**2*q2+e**4*q4
    action = H0*q**2/2+T0*q**3/6+e**2*(J2*q+K2*q**2/2)
    assert sp.expand(action).coeff(e, 4) == -J2**2/(2*H0)
    assert sp.factor(sp.expand(action).coeff(e, 6)-K2*q2**2/2-T0*q2**3/6) == 0
    assert sp.expand(action).coeff(e, 6) != 0
