import pytest
import sympy as sp
from p8_m1 import nonlinear, robustness


def test_exact_normal_form_charts_and_covariant_deformation():
    assert all(value == 0 for value in robustness.audit_checks().values())


def test_positive_diagonal_buffer_does_not_cure_mismatch():
    data = robustness.controls()
    K, G = data["bounce"], data["mismatch_gradient"]
    assert K[0, 0] > 0 and K.det() > 0 and G[0, 0] > 0 and G.det() > 0
    assert K[0, 0]-G[0, 0] == sp.Rational(1, 100)
    nf = robustness.normal_form(K, G)
    assert nf["sigma"] == 0 and nf["r2"] == sp.Rational(1, 1199)
    assert data["mismatch_speeds"][1] > 1


@pytest.mark.parametrize("a,b", [(0, 0), (sp.Rational(1, 4), 0), (-sp.Rational(1, 4), 0),
                               (sp.Rational(1, 4), sp.Rational(1, 100)), (0, sp.Rational(1, 10))])
def test_normal_form_health_and_causality_controls(a, b):
    K = sp.eye(2)
    G = sp.Matrix([[1-a, -b], [-b, 1]])
    nf = robustness.normal_form(K, G)
    high = 1-nf["sigma"]/2+sp.sqrt(nf["sigma"]**2/4+nf["r2"])
    if b != 0 or a < 0:
        assert high > 1
    else:
        assert high == 1


def test_matter_diagonal_changes_are_outside_the_normal_form_contract():
    K = sp.eye(2)
    G = sp.Matrix([[sp.Rational(9, 10), sp.Rational(1, 20)],
                   [sp.Rational(1, 20), sp.Rational(9, 10)]])
    assert sorted(G.eigenvals()) == [sp.Rational(17, 20), sp.Rational(19, 20)]
    with pytest.raises(ValueError, match="Equal"):
        robustness.normal_form(K, G)


@pytest.mark.parametrize("K,G", [(sp.eye(3), sp.eye(3)), (sp.Matrix([[1, 1], [0, 1]]), sp.eye(2)),
                                (sp.zeros(2), sp.zeros(2)), (sp.ones(2), sp.ones(2))])
def test_bad_or_degenerate_matrix_inputs_rejected(K, G):
    with pytest.raises(ValueError):
        robustness.normal_form(K, G)


def test_covariant_deformation_sign_and_kinetic_boundary():
    d = nonlinear.d
    epsilon = robustness.deformation()["epsilon"]
    speed = 1/(1+epsilon/d**2)
    assert speed.subs({nonlinear.u: 0, epsilon: 1}) == sp.Rational(1, 2)
    assert speed.subs({nonlinear.u: 0, epsilon: -sp.Rational(1, 2)}) == 2
    assert (1+epsilon/d**2).subs({nonlinear.u: 0, epsilon: -1}) == 0
    assert sp.limit(speed, nonlinear.u, sp.oo) == 1


def test_covariant_deformation_has_quantitative_derivative_tail_bounds():
    data = robustness.deformation_tail_bounds()
    assert all(value == 0 for value in data["residuals"].values())
    assert [row["decay_power_n"] for row in data["rows"]] == [3, 2, 1]
    assert [row["coefficient_bound_per_abs_epsilon"] for row in data["rows"]] == ["2", "1/2", "1/32"]
    for row in data["rows"]:
        assert sp.Rational(row["first_derivative_squared_bound_per_epsilon_squared"]) > 0
        assert sp.Rational(row["second_derivative_bound_per_abs_epsilon"]) > 0
        assert row["X_and_mixed_derivatives"] == "0"
    import json
    old = json.loads(nonlinear.WITNESS.read_text())
    assert sp.Rational(old["canonical_kinetic_ge_1_over_4_when_d_ge"]) == sp.Rational(10481, 1680)
