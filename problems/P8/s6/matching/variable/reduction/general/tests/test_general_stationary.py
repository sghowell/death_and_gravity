"""General-beta tests, including literal Fraction square-root/principal minors."""

from fractions import Fraction as F
from itertools import combinations, permutations

import pytest
import sympy as sp
from p8_general_reduction import general as g


class Jet:
    """Independent truncated Fraction ring through the cubic coefficient."""

    def __init__(self, value=0):
        if isinstance(value, Jet):
            self.a = value.a
        elif isinstance(value, tuple):
            self.a = tuple(F(x) for x in value)+tuple(F(0) for _ in range(4-len(value)))
        else:
            self.a = (F(value), F(0), F(0), F(0))

    def __add__(self, other):
        other = Jet(other)
        return Jet(tuple(x+y for x, y in zip(self.a, other.a, strict=True)))

    __radd__ = __add__

    def __neg__(self):
        return Jet(tuple(-x for x in self.a))

    def __sub__(self, other):
        return self+-Jet(other)

    def __rsub__(self, other):
        return Jet(other)+-self

    def __mul__(self, other):
        other = Jet(other)
        return Jet(tuple(sum(self.a[j]*other.a[n-j] for j in range(n+1)) for n in range(4)))

    __rmul__ = __mul__


def _multiply(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(len(right)))
             for j in range(len(right[0]))] for i in range(len(left))]


def _determinant(matrix):
    if not matrix:
        return Jet(1)
    answer = Jet(0)
    for perm in permutations(range(len(matrix))):
        parity = sum(perm[i] > perm[j] for i in range(len(perm)) for j in range(i+1, len(perm)))
        term = Jet((-1)**parity)
        for i, j in enumerate(perm):
            term *= matrix[i][j]
        answer += term
    return answer


def _literal_fraction_potential(ratio, beta):
    lower = [[F(2), F(1, 2), F(-1, 3), F(2, 5)],
             [F(1, 2), F(-1), F(3, 4), F(1, 7)],
             [F(-1, 3), F(3, 4), F(4), F(-1, 6)],
             [F(2, 5), F(1, 7), F(-1, 6), F(-2)]]
    signs = (1, -1, -1, -1)
    h = [[signs[i]*lower[i][j] for j in range(4)] for i in range(4)]
    h2, h3 = _multiply(h, h), _multiply(_multiply(h, h), h)
    root = [[Jet((ratio*int(i == j), ratio*h[i][j]/2,
                  -ratio*h2[i][j]/8, ratio*h3[i][j]/16))
             for j in range(4)] for i in range(4)]
    squared = _multiply(root, root)
    for i in range(4):
        for j in range(4):
            assert squared[i][j].a == (ratio*ratio*int(i == j), ratio*ratio*h[i][j], 0, 0)
    elementary = []
    for n in range(5):
        en = Jet(0)
        for subset in combinations(range(4), n):
            en += _determinant([[root[i][j] for j in subset] for i in subset])
        elementary.append(en)
    potential = -2*sum(beta[n]*elementary[n] for n in range(5))
    return lower, elementary, potential


def test_all_exact_identities_and_the_original_nineteen():
    checks = g.checks()
    assert len(checks) == 45
    assert sum(key.startswith("core") for key in checks) == 19
    assert set(checks.values()) == {0}


@pytest.mark.parametrize("ratio,b1,b2,b3", [
    (F(2), F(3), F(1, 2), F(1, 4)),
    (F(2), F(-3), F(-1, 2), F(-1, 4)),
    (F(3, 2), F(2), F(-1, 3), F(2, 5)),
    (F(1), F(1), F(-2, 3), F(1, 3)),
])
def test_literal_general_square_root_principal_minor_series(ratio, b1, b2, b3):
    b4 = -(b1+3*ratio*b2+3*ratio*ratio*b3)/ratio**3
    beta = (F(7), b1, b2, b3, b4)
    lower, elementary, literal = _literal_fraction_potential(ratio, beta)
    primary = g.potential()
    at = {g.RATIO: sp.Rational(ratio)}
    at.update({primary["h"][i, j]: sp.Rational(lower[i][j]) for i in range(4) for j in range(i, 4)})
    for n, series in enumerate(elementary):
        assert tuple(sp.factor(value.subs(at)) for value in primary["elementary_series"][n]) == series.a
    at.update(dict(zip(g.BETA, (sp.Rational(value) for value in beta), strict=True)))
    assert tuple(sp.factor(value.subs(at)) for value in primary["literal_potential_coefficients"]) == literal.a
    data = g.polynomial_at(beta, ratio)
    at.update({g.P: data["P"], g.Z: data["Z"]})
    assert literal.a[1] == 0
    assert sp.factor(primary["quadratic"].subs(at)) == literal.a[2]
    assert sp.factor(primary["cubic"].subs(at)) == literal.a[3]


def test_root_normalization_and_not_the_full_vacuum_equation():
    data = g.potential()
    assert sp.factor(data["F_root"].subs(data["root_rule"])) == 0
    assert sp.diff(data["F_root"], g.BETA[0]) == 0
    assert sp.diff(data["stationary_potential_polynomial"], g.BETA[0]) == 1
    assert "not necessarily a vacuum" in g.calibration()["root_scope"]


def test_full_ten_component_hessian_and_one_third_trace_inverse():
    p, data = g.potential(), g.stationary()
    assert p["hessian"].shape == (10, 10)
    assert p["hessian_determinant"] == 3*g.P**10*g.RATIO**10/16
    trace = sp.trace(g.ETA*data["B"])
    wrong = g.MF2*g.RATIO/g.P*(data["B"]-g.ETA*trace/2)
    assert sp.expand(sp.trace(g.ETA*wrong)-sp.trace(g.ETA*data["relative_h2"])) != 0
    assert p["hessian"].subs(g.P, 0) == sp.zeros(10)


@pytest.mark.parametrize("sign", [1, -1])
def test_regular_positive_and_negative_P_are_algebraic_not_health_claims(sign):
    beta = (7, 3*sign, F(sign, 2), F(sign, 4), F(-9*sign, 8))
    data = g.inverse_at(beta, 2, mg2=2, mf2=5)
    assert data["P"] == 6*sign
    assert data["Z"] == 8*sign
    assert data["kappa"] == sign*sp.Rational(25, 3)
    assert data["leading_Planck_squared"] == 22
    assert data["additive_inverse_coefficient"] == sign*sp.Rational(20, 3)
    assert data["relative_inverse_coefficient"] == sign*sp.Rational(5, 3)
    assert "not a parent-health verdict" in data["scope"]
    assert data["kappa"] != data["MF2"]*data["r"]**3/(4*data["P"])


def test_double_root_keeps_polynomial_cubic_and_rejects_only_the_inverse():
    beta = (0, 1, F(-2, 3), F(1, 3), 0)
    data = g.polynomial_at(beta, 1)
    assert data["F_root"] == data["P"] == data["hessian_determinant_on_root"] == 0
    assert data["Z"] == -sp.Rational(1, 3)
    p = g.potential()
    at = dict(zip(g.BETA, beta, strict=True))
    assert sp.factor(p["F_root"].subs(at)) == (g.RATIO-1)**2
    # H=diag(1,1,1,0), with covariant h=eta H, gives a nonzero cubic.
    matrix = g.ETA*sp.diag(1, 1, 1, 0)
    point = {g.RATIO: 1, g.P: 0, g.Z: -sp.Rational(1, 3)}
    point.update({p["h"][i, j]: matrix[i, j] for i in range(4) for j in range(i, 4)})
    assert sp.factor(p["cubic"].subs(point)) == -sp.Rational(1, 12)
    with pytest.raises(ValueError, match="P=0"):
        g.inverse_at(beta, 1)


def test_endpoint_only_disconnected_case_is_not_algebraically_integrated_out():
    data = g.polynomial_at((3, 0, 0, 0, 0), F(7, 3))
    assert data["F_root"] == data["P"] == data["Z"] == 0
    with pytest.raises(ValueError, match="P=0"):
        g.inverse_at((3, 0, 0, 0, 0), F(7, 3))


def test_all_variable_kappa_and_r_scalar_terms_are_retained():
    data = g.leading_action()
    assert sp.diff(data["F2"], data["X"]) == data["F2_X"]
    assert data["A1"] == 2*data["F2_X"]
    assert data["A2"] == -data["A1"]
    assert data["A3"] == data["A4"] == data["A5"] == 0
    assert data["K"].subs({g.KP: 0, g.KPP: 0}) == 0
    assert data["f4"].subs({g.KP: 0, g.KPP: 0, g.KPPP: 0}) == 0
    assert sp.expand(data["f4"]).has(g.KPPP)
    assert "retained" in g.leading_action.__doc__
    assert "free chi unchanged" in data["matter"]


def test_unequal_Einstein_conformal_kinetic_sign():
    time = sp.Symbol("time", real=True)
    a, r = sp.Function("a")(time), sp.Function("r")(time)
    hubble = sp.diff(a, time)/a
    hf = sp.diff(a*r, time)/(a*r*r)
    rg = -6*(sp.diff(hubble, time)+2*hubble**2)
    rf = -6*(sp.diff(hf, time)/r+2*hf**2)
    literal = -g.MF2*a**3*r**4*rf/2
    retained = -g.MF2*a**3*r*r*rg/2-3*g.MF2*a**3*sp.diff(r, time)**2
    boundary = 3*g.MF2*sp.diff(a**3*r*sp.diff(r, time), time)
    assert sp.factor(literal-retained-boundary) == 0
    assert sp.factor(literal-retained+boundary) != 0


def test_smooth_deformation_preserves_lower_action_but_changes_cubic():
    freedom, data = g.cubic_freedom(), g.potential()
    for key in ("F_root", "P_polynomial", "stationary_potential_polynomial"):
        assert sp.expand(data[key].subs(freedom["rule"], simultaneous=True)-data[key]) == 0
    beta = (7, 3, F(1, 2), F(1, 4), F(-9, 8))
    eta, ratio = F(3, 7), F(2)
    delta = (3*ratio*eta, -2*eta, eta/ratio, 0, -eta/ratio**3)
    other = tuple(x+y for x, y in zip(beta, delta, strict=True))
    before, after = g.inverse_at(beta, ratio, 2, 5), g.inverse_at(other, ratio, 2, 5)
    for key in ("P", "kappa", "relative_inverse_coefficient", "leading_Planck_squared",
                "stationary_potential_polynomial"):
        assert before[key] == after[key]
    assert after["Z"]-before["Z"] == eta
    assert freedom["cubic_shift"] != 0


def test_full_table_is_not_only_its_s_zero_projection():
    data = g.curvature_linear_cubic()
    assert len(data["operators"]) == len(data["coefficients"]) == 10
    assert len(data["retained_complement"]) == 8
    assert set(data["operators"])-set(data["retained_complement"]) == {"X_G_vv", "R_X_squared"}
    assert sp.expand(data["raw_over_N"]-data["table_over_N"]) == 0
    assert sp.diff(data["coefficients"]["G_H_squared"], g.Z) == 8*g.S**2
    assert sp.diff(data["coefficients"]["G_v_Hv"], g.Z) != 0
    assert "not full S6" in data["scope"]


def test_vanishing_center_coefficient_cannot_be_dropped_before_IBP():
    data = g.curvature_linear_cubic()
    coefficient = data["N"]*data["coefficients"]["X_G_H"]
    assert coefficient.subs(g.S, 0) == 0
    # All derivatives of its other coefficients multiply s at this point.
    derivative_at_center = sp.diff(coefficient, g.S).subs(g.S, 0)*g.SPHI
    assert sp.factor(derivative_at_center+8*data["N"]*(g.P+g.Z)*g.SPHI**2) == 0
    assert derivative_at_center != 0
    assert sp.diff(coefficient.subs(g.S, 0), g.S) == 0  # deliberately invalid projection order


def test_raw_center_Z_cancellation_and_sign_do_not_define_a_full_invariant():
    data = g.curvature_linear_cubic()
    center = data["s_zero_coefficient_X_G_vv"]
    assert sp.diff(center, g.Z) == 0
    assert sp.factor(center-g.MF2**3*g.RATIO**4*g.SPHI**2/(2*g.P**2)) == 0
    assert sp.factor(data["H_free_projection_Xi"].subs(g.S, 0)-center*data["X"]) == 0


def test_nonboundary_deformation_uses_compact_variation_not_only_density():
    data = g.nonboundary_control()
    assert sp.expand(data["literal_conformal_variation"]-data["expected_variation"]) == 0
    assert sp.factor(sp.diff(data["literal_conformal_variation"], data["sigma"])
                     +2*data["delta_L6"]) == 0
    scale = sp.Symbol("constant_conformal_scale", positive=True)
    scaled_density = scale**4*data["delta_L6"].subs(data["Lambda"], data["Lambda"]/scale**2)
    assert sp.factor(sp.diff(scaled_density, scale).subs(scale, 1)+2*data["delta_L6"]) == 0
    nonzero = {g.RATIO: 1, g.MF2: 1, g.P: 1, data["eta"]: 1, data["Lambda"]: -3}
    assert data["delta_L6"].subs(nonzero) == -1
    assert data["covariant_metric_Euler_trace"].subs(nonzero) == 1
    assert "off-shell physical g" in data["scope"]


@pytest.mark.parametrize("bad", [True, False, 1.0, "1", sp.oo, sp.nan, sp.sqrt(2), sp.Symbol("x")])
def test_exact_rational_interface_rejects_inexact_or_unresolved_values(bad):
    with pytest.raises(TypeError):
        g.polynomial_at((0, bad, 0, 0, -1), 1)
    with pytest.raises(TypeError):
        g.polynomial_at((0, 1, 0, 0, -1), bad)


@pytest.mark.parametrize("bad_beta", [(), (0, 1), (0, 1, 0, 0, -1, 3), "beta"])
def test_beta_length_and_container_guard(bad_beta):
    with pytest.raises(ValueError):
        g.polynomial_at(bad_beta, 1)


@pytest.mark.parametrize("ratio", [0, -1])
def test_positive_root_chart_guard(ratio):
    with pytest.raises(ValueError):
        g.polynomial_at((0, 1, 0, 0, -1), ratio)


@pytest.mark.parametrize("mg2,mf2", [(0, 1), (1, 0), (-1, 1), (1, -1)])
def test_positive_Einstein_inverse_guard(mg2, mf2):
    with pytest.raises(ValueError):
        g.inverse_at((0, 1, 0, 0, -1), 1, mg2, mf2)


def test_off_root_data_is_not_a_stationary_inverse():
    data = g.polynomial_at((0, 1, 0, 0, -1), 2)
    assert data["F_root"] == -7
    with pytest.raises(ValueError, match="not a stationary"):
        g.inverse_at((0, 1, 0, 0, -1), 2)
