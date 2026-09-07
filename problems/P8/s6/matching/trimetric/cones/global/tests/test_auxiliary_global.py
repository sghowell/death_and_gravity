import pytest
import sympy as sp
from p8_trimetric_global import (
    background,
    guards,
    independent,
    monotonic,
    obstruction,
    vacuum,
)


@pytest.mark.parametrize("module", (background, vacuum, obstruction, guards, monotonic))
def test_exact_source_background_and_scope_identities(module):
    assert all(sp.simplify(value) == 0 for value in module.checks().values())


def test_one_nonzero_link_is_enough_without_adding_a_disconnected_clock_to_K():
    d = monotonic.derive(1)
    assert d["K"] == d["links"][0]["G"]/d["links"][0]["R"]**2
    assert all(sp.simplify(value) == 0 for value in monotonic.checks(1).values())


@pytest.mark.parametrize("bad", (0, 3, True, 1.0))
def test_empty_or_inexact_nonzero_link_index_sets_are_rejected(bad):
    with pytest.raises(ValueError):
        monotonic.derive(bad)


def test_normalizer_is_positive_without_a_link_sign_or_TT_inverse_assumption():
    d = monotonic.derive()
    assert d["K"].is_positive is True
    assert all(link["p"] not in d["K"].free_symbols for link in d["links"])
    point = {d["links"][0]["p"]: 1, d["links"][1]["p"]: -1,
             d["links"][0]["R"]: 1, d["links"][1]["R"]: 1}
    assert d["K"].subs(point).is_positive is True


def test_independent_fraction_source_dynamics_and_actual_solution_fixtures():
    r = independent.checks()
    assert len(r["coefficientwise_identities"]) == 5
    assert len(r["off_shell_source_dynamics_normalization_fixtures"]) == 5
    assert len(r["actual_positive_NEC_rolling_fixtures"]) == 4
    assert len(r["same_action_mixed_weight_fixtures"]) == 3
    assert {fixture["q"] for fixture in r["actual_positive_NEC_rolling_fixtures"]} == {"1", "-1"}


def test_bounce_equation_does_not_divide_by_H_u():
    d = monotonic.derive()
    at_bounce = d["combined_residual"].subs({d["H_u"]: 0, d["H_u_prime"]: 1})
    assert sp.simplify(at_bounce-2*d["K"]-d["n_h"]) == 0
    assert at_bounce.is_positive is True


def test_degenerate_crossing_is_not_decided_only_by_acceleration_at_the_middle():
    t = sp.Symbol("T", real=True)
    h = t**3
    assert sp.diff(h, t).subs(t, 0) == 0
    required_null = -2*sp.diff(h, t)  # K=1 kinematic diagnostic, not a parent solution.
    assert required_null == -6*t**2
    assert required_null.subs(t, 1) < 0


def test_generic_residual_budget_is_a_threshold_not_an_actual_correction_bound():
    assert monotonic.bounce_residual_budget(1, 1, 1)["combined_residual_lower_bound"] == 3
    assert monotonic.bounce_residual_budget(1, 1, 1, sp.Rational(1, 10), 5)["combined_residual_lower_bound"] == sp.Rational(5, 2)
    assert monotonic.bounce_residual_budget(1, 1, 1, 1, 5)["status"] == "INCONCLUSIVE"


@pytest.mark.parametrize("bad", (True, 0.1, float("inf"), sp.Float("0.1"), sp.oo, sp.I))
def test_exact_budget_interfaces_reject_inexact_and_nonreal_inputs(bad):
    with pytest.raises(TypeError):
        monotonic.bounce_residual_budget(bad, 1)
    with pytest.raises(TypeError):
        obstruction.negative_link_three_slice(bad, 1)


def test_explicit_actual_physical_CD_endpoint_dictionary_and_strictness():
    half = monotonic.cd_endpoint_test(sp.Rational(1, 2), 1)
    assert half["necessary_normalized_endpoint_error"] == sp.Rational(8, 5)
    assert half["status"] == "EXCLUDED_BY_ACTUAL_PARENT_HUBBLE_MONOTONICITY"
    assert monotonic.cd_endpoint_test(1, 0)["necessary_normalized_endpoint_error"] == 2
    assert monotonic.cd_endpoint_test(sp.Rational(1, 2), sp.Rational(8, 5))["status"] == "INCONCLUSIVE"


def test_complementary_negative_link_three_slice_test_uses_actual_cap_ratio():
    assert obstruction.negative_link_three_slice(sp.Rational(1, 2), 1)["strict_chord_margin"] == 1
    assert obstruction.negative_link_three_slice(sp.Rational(1, 4), 1)["status"] == "INCONCLUSIVE"
    assert obstruction.negative_link_three_slice(sp.Rational(1, 2), 1, sp.Rational(1, 3))["status"] == "INCONCLUSIVE"


@pytest.mark.parametrize("args", ((0, 1), (2, 1), (1, 0), (1, 1, -1), (1, 1, 1)))
def test_three_slice_domain_is_not_relaxed_by_a_failed_strict_test(args):
    with pytest.raises(ValueError):
        obstruction.negative_link_three_slice(*args)


def test_arbitrary_positive_vacuum_ratios_use_actual_h_kinetics():
    c = vacuum.controls()
    assert tuple(c.values()) == (1, -2, -1, 2, 3, 4, sp.Rational(10, 3))


def test_negative_cap_does_not_forge_an_actual_vacuum():
    c = guards.controls()
    assert c["same_example_calibrated_E_e00"] == c["same_example_calibrated_E_v00"] == 0
    assert c["positive_cap_does_not_supply_full_vacuum_E_u00"] == -2


def test_true_mixed_singular_flat_control_is_covered_by_background_not_inverse_theorem():
    assert guards.controls()["singular_TT_hessian"] == 0
    assert all(guards.checks()[f"mixed_singular_full_flat_{key}"] == 0 for key in ("E_e", "E_v", "E_u"))


def test_all_zero_links_have_a_real_but_underdetermined_bouncing_metric_control():
    assert guards.controls()["all_zero_links_arbitrary_metric_H_prime_at_zero"] == 2
    assert all(guards.checks()[f"all_zero_links_full_arbitrary_u_{key}"] == 0 for key in ("E_e", "E_v", "E_u"))


def test_primary_theorem_allows_actual_nonzero_link_deSitter_without_sign_change():
    assert guards.controls()["actual_deSitter_H_u"] == sp.Rational(1, 2)
    assert all(guards.checks()[f"actual_deSitter_nonzero_link_{key}"] == 0 for key in ("E_e", "E_v", "E_u"))
