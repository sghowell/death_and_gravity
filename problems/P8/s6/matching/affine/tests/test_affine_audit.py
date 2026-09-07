"""Separately authored source/signature/IBP audit; no primary lift import."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_affine import audit as a


@pytest.mark.parametrize("name", tuple(a.identities()))
def test_exact_audit_identities(name):
    assert a.identities()[name] == 0


def test_all_continuous_proof_controls():
    assert all(a.checks().values())


def test_five_target_relations_and_dependent_second_coefficient():
    d, p = a.parameters(), a.principal()["simplified"]
    assert tuple(a.principal_checks()) == (
        "target_f", "target_alpha1", "target_alpha3", "target_alpha4", "target_alpha5"
    )
    assert p["alpha2"] == 0
    assert sp.simplify(p["f"] - d["w"] / (2 * a.H)) == 0
    assert sp.simplify(p["alpha3"] - 1 / (a.H * a.X)) == 0
    assert sp.simplify(p["alpha4"] - (11 * a.X - 4 * a.H + 4)
                       / (4 * a.H * a.X * d["w"])) == 0
    assert sp.simplify(p["alpha5"] + 1 / (a.H * a.X * d["w"])) == 0


def test_direct_bridge_to_frozen_cd_target_not_primary_affine_dictionary():
    from p8 import matter, rational_candidates

    source = rational_candidates.specification("CD_matter")
    clock = rational_candidates.j.t
    replace = {a.H: (1 + clock**2)**3, a.Y: rational_candidates.X}
    expected = a.target()["repository"]
    assert sp.cancel(expected["F2"].subs(replace) - source["F2"]) == 0
    assert sp.cancel(expected["A1"].subs(replace) - source["A1"]) == 0
    assert sp.cancel(expected["A3"].subs(replace) - source["A3"]) == 0
    a2, a4, a5 = matter.ia_completion(
        source["F2"], sp.diff(source["F2"], rational_candidates.X),
        source["A1"], source["A3"], rational_candidates.X,
    )
    for key, value in (("A2", a2), ("A4", a4), ("A5", a5)):
        assert sp.cancel(expected[key].subs(replace) - value) == 0


@pytest.mark.parametrize("number", range(1, 7))
def test_signature_from_dense_traceful_hessians(number):
    v = tuple(Fraction(number + j, j + 1) for j in range(4))
    hessian = tuple(tuple(Fraction((i + 1) * (j + 1) - number, i + j + 1)
                          for j in range(4)) for i in range(4))
    repo = a.contractions(v, hessian)
    paper = a.contractions(v, hessian, paper=True)
    for name, sign in (("X", -1), ("Box", -1), ("vHv", 1),
                       ("L1", 1), ("L2", 1), ("L3", -1), ("L4", -1), ("L5", 1)):
        assert paper[name] == sign * repo[name]


def test_wrong_signature_for_a3_is_not_an_equivalent_cd_coefficient():
    actual = a.principal()["simplified"]["alpha3"].subs({a.H: 1, a.X: -1})
    required_repo = a.target()["repository"]["A3"].subs({a.H: 1, a.Y: 1})
    assert -actual == required_repo == 1
    assert actual != required_repo


@pytest.mark.parametrize("paper", (False, True))
@pytest.mark.parametrize("number", range(1, 5))
def test_literal_product_rule_and_source_signature(paper, number):
    vector = (Fraction(5, 2), 1, Fraction(-2, 3), number)
    hessian = ((2, 1, -1, 0), (1, 4, 2, 1), (-1, 2, -3, 5), (0, 1, 5, 7))
    output = a.divergence_fixture(
        vector, hessian, Fraction(number, 3), Fraction(4, 5), Fraction(-5, 7), paper=paper
    )
    assert output["residual"] == 0
    assert output["direct"] == output["expanded"]


def test_printed_q2_failure_is_off_trajectory_and_exact():
    d = a.lower_order()
    error = d["printed_error"].subs(d["q"], 0)
    assert sp.simplify(error.subs({a.U: 1, a.X: -sp.Rational(9, 10)})) != 0
    # Testing only q=0 on x=-1 misses the printed error.
    assert sp.simplify(error.subs(a.X, -1)) == 0
    assert d["lower_order_claim_conditional_on_literal_Q2"] is True


def test_pure_pR_cross_coefficient_from_literal_gradient_contractions():
    v = (Fraction(3, 2), Fraction(1, 5), -1, Fraction(2, 7))
    hessian = ((3, 1, -2, 4), (1, -3, 2, 1), (-2, 2, 5, -1), (4, 1, -1, 2))
    tensor = a.contractions(v, hessian, paper=True)
    pu, px, p = Fraction(4, 5), Fraction(-2, 3), Fraction(3, 2)
    metric = (-1, 1, 1, 1)
    raised = tuple(metric[i] * v[i] for i in range(4))
    gradp = tuple(pu * v[i] + 2 * px * sum(hessian[i][j] * raised[j] for j in range(4))
                  for i in range(4))
    literal = Fraction(3, 2) * sum(metric[i] * gradp[i]**2 for i in range(4)) / p
    expanded = (Fraction(3, 2) * pu**2 * tensor["X"] / p
                + 6 * pu * px * tensor["vHv"] / p
                + 6 * px**2 * tensor["L4"] / p)
    assert literal == expanded
    d = a.palatini_control()
    fixture = {d["p"]: sp.Rational(3, 2), d["p_u"]: sp.Rational(4, 5),
               d["p_x"]: -sp.Rational(2, 3), a.X: -1}
    assert sp.simplify((d["Q2_corrected"] - d["cross_coefficient"]).subs(fixture)) == 0
    assert sp.simplify((d["Q2_printed"] - d["cross_coefficient"]).subs(fixture)) != 0


def test_initial_parent_cubic_and_scalar_have_no_spurious_clock_boundary_shift():
    d = a.lower_order()
    zero_q = {d["q"]: 0}
    cubic = d["F3"].subs(zero_q).doit().subs(a.X, -1)
    scalar = (d["F2_parent"] - d["F_repository"]).subs(zero_q).doit().subs(a.X, -1)
    assert sp.simplify(cubic) == 0
    assert sp.simplify(scalar) == 0


@pytest.mark.parametrize("u", (0, 1, -1, Fraction(3, 2), Fraction(-3, 2), 1000, -1000))
@pytest.mark.parametrize("x", (Fraction(-11, 10), -1, Fraction(-9, 10)))
def test_exact_domain_fixtures_only_supplement_continuous_proof(u, x):
    d = a.domain_at(u, x)
    assert Fraction(9, 20) <= d["f"] <= Fraction(11, 20)
    assert Fraction(9, 40) <= d["p_squared"] <= Fraction(11, 40)
    assert d["quotient_factor"] >= Fraction(4, 5)
    assert d["denominator"] == 1
    assert 0 < d["q_abs_bound"] < Fraction(1, 32)


def test_closed_tube_endpoint_rejects_falsely_strict_quotient_bound():
    d = a.domain_at(0, Fraction(-9, 10))
    assert d["quotient_factor"] == Fraction(4, 5)
    assert not d["quotient_factor"] > Fraction(4, 5)


@pytest.mark.parametrize("u,x", ((True, -1), (0, False), (1.0, -1),
                                (0, -1.0), (sp.oo, -1), (0, sp.nan)))
def test_exact_domain_api_rejects_inexact_or_nonfinite_inputs(u, x):
    with pytest.raises(TypeError):
        a.domain_at(u, x)


@pytest.mark.parametrize("x", (Fraction(-6, 5), Fraction(-4, 5), 0, 1))
def test_domain_api_does_not_admit_open_extension_as_closed_certificate(x):
    with pytest.raises(ValueError):
        a.domain_at(0, x)


def test_tensor_input_guards():
    identity = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
    with pytest.raises(ValueError):
        a.contractions((1, 2, 3), identity)
    with pytest.raises(ValueError):
        a.contractions((1, 2, 3, 4), ((1,),))
    with pytest.raises(ValueError):
        a.contractions((1, 2, 3, 4), ((1, 1, 0, 0), (0, 1, 0, 0),
                                    (0, 0, 1, 0), (0, 0, 0, 1)))
    with pytest.raises(TypeError):
        a.contractions((True, 0, 0, 0), identity)
    with pytest.raises(TypeError):
        a.contractions((1, 2, 3, 4), identity, paper=1)


def test_audit_does_not_claim_connection_rank_or_uv_from_normal_form():
    d = a.calibration()
    assert d["source_Q2_correction_required"] is True
    assert d["connection_rank_or_health_proved_here"] is False
    assert d["UV_or_cutoff_claim"] is False
