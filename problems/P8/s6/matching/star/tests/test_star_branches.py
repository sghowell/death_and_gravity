import pytest
import sympy as sp
from p8_star import background, branches, independent


def test_star_dynamic_clock_and_sum_identities():
    assert all(value == 0 for value in background.checks().values())


def test_star_all_rational_star_sizes_match_independent_chain_rule():
    for fixture in independent.dynamic_fixtures():
        actual = background.dynamic_reconstruction(
            fixture["Gs"], fixture["Rs"], fixture["cs"], fixture["H"], fixture["nh"],
            central_G=fixture["G_u"])
        for key in ("K", "Kprime", "Hprime", "Rprimes"):
            assert actual[key] == fixture[key]
        assert actual["scaled_H_prime"].is_nonpositive is True


def test_star_flat_GR_limit_has_only_central_kinetic_weight():
    actual = background.dynamic_reconstruction([], [], [], -7, 6, central_G=3)
    assert actual["K"] == 3 and actual["Kprime"] == 0 and actual["Hprime"] == -1


def test_star_bounce_is_excluded_without_dividing_H_or_link_signs():
    actual = background.dynamic_reconstruction([1, 2], [2, 3], [7, "1/3"], 0, 5)
    assert actual["Kprime"] == 0 and actual["Hprime"] < 0


def test_star_endpoint_only_is_not_a_genuine_algebraic_root_branch():
    assert branches.classify([7, 0, 0, 0, -5]) == {
        "kind": "endpoint_only", "positive_roots": (), "identically_zero": True}
    assert branches.classify([0, 0, 0, -2, 0]) == {
        "kind": "genuine", "positive_roots": (), "identically_zero": False}


def test_star_linear_double_and_two_distinct_positive_roots():
    assert branches.classify([0, 1, "-1/2", 0, 0])["positive_roots"] == ((1, 1),)
    assert branches.classify([0, 1, -1, 1, 0])["positive_roots"] == ((1, 2),)
    assert branches.classify([0, 2, "-3/2", 1, 0])["positive_roots"] == ((1, 1), (2, 1))
    assert branches.classify([0, 1, 0, 1, 0])["positive_roots"] == ()


def test_star_irrational_positive_roots_are_not_dropped_or_float_approximated():
    assert branches.classify([0, -2, 0, 1, 0])["positive_roots"] == ((sp.sqrt(2), 1),)


def test_star_algebraic_fixture_is_an_actual_solution_and_invalidates_wrong_K():
    control = branches.algebraic_control()
    assert all(value == 0 for value in control["residuals"].values())
    assert all(value == 0 for value in branches.checks().values())
    assert control["false_K_all_defect"].is_negative is True
    for fixture in independent.actual_branch_fixtures():
        assert fixture["wrong_all_K_undivided_defect"] < 0


def test_star_genuinely_disconnected_center_can_bounce_only_outside_theorem():
    control = branches.endpoint_control()
    assert control["central_Euler"] == control["leaf_Euler"] == control["scalar_Euler"] == 0
    assert control["H_u_prime_at_zero"] == 2
    assert independent.disconnected_fixtures()[1]["H_u_prime"] == 2
    with pytest.raises(ValueError):
        background.dynamic_reconstruction([], [], [], 0, 0)


def test_star_compact_logarithmic_rate_envelope_is_not_a_solution_estimate():
    assert background.compact_comparison_bound(["1/2", 7], [8, "3/4"]) == 8
    assert background.compact_comparison_bound([], []) == 0
    with pytest.raises(ValueError):
        background.compact_comparison_bound([1], [])
    with pytest.raises(ValueError):
        background.compact_comparison_bound([-1], [1])


@pytest.mark.parametrize("kwargs", [
    {"Gs": [-1], "Rs": [1], "cs": [1]},
    {"Gs": [0], "Rs": [1], "cs": [1]},
    {"Gs": [1], "Rs": [0], "cs": [1]},
    {"Gs": [1], "Rs": [1], "cs": [-1]},
    {"Gs": [1, 2], "Rs": [1], "cs": [1]},
    {"Gs": [1], "Rs": [1], "cs": [1], "nh": -1},
    {"Gs": [1], "Rs": [1], "cs": [1], "central_G": -1},
])
def test_star_dynamic_public_domain_guards(kwargs):
    arguments = {"H": 0, "nh": 0} | kwargs
    with pytest.raises((TypeError, ValueError)):
        background.dynamic_reconstruction(**arguments)
