import pytest
import sympy as sp
from p8_determinant import background, controls, independent, potential


def test_determinant_all_exact_primary_residuals():
    for module in (potential, background, controls):
        assert all(value == 0 for value in module.checks().values())


def test_determinant_full_coframe_and_polynomial_independent_replay():
    report = independent.checks()
    assert report["full_coframe_first_jet_directions"] == 192
    assert report["full_auxiliary_stationary_directions"] == 16
    assert len(report["coefficientwise_identities"]) == 10


def test_determinant_primary_stresses_match_every_independent_matrix_fixture():
    for row in independent.stress_fixtures():
        actual = potential.evaluate(row["beta"], row["a"], row["n"], lambda_=row["lambda"])
        for key in ("S_a", "S_N", "rho", "pressure", "nulls"):
            assert actual[key] == row[key]


def test_determinant_clock_measures_match_independent_dynamic_rates():
    for row in independent.rate_fixtures():
        actual = background.reconstruct(row["Gs"], row["ys"], row["cs"], row["H"], row["nulls"])
        for key in ("K", "Kprime", "Hprime", "yprimes", "weighted_null"):
            assert actual[key] == row[key]
        assert actual["scaled_H_prime"].is_nonpositive is True


def test_determinant_actual_stiff_and_vacuum_controls_have_full_equations():
    stiff, vacuum = controls.proportional_stiff(), controls.proportional_de_sitter()
    assert stiff["cosmological_constants"] == (-8, sp.Rational(1, 2), -sp.Rational(8, 81))
    assert stiff["K"] == 36 and stiff["S_y"] == 2
    assert stiff["scaled_H_prime"].is_negative is True
    assert all(value == 0 for value in stiff["residuals"].values())
    assert all(value == 0 for value in vacuum["residuals"].values())
    assert vacuum["scaled_H_prime"] == 0


def test_determinant_singular_time_sum_is_not_admitted_by_regular_flag():
    control = controls.singular_sum_bianchi_control()
    actual = potential.evaluate(control["beta"], control["a"], control["n"])
    assert (actual["S_a"], actual["S_N"]) == (2, 0)
    assert actual["regular_sum"] is False
    assert actual["pressure"] == (0, 0, 0)
    assert any(actual["rho"])
    assert control["bianchi"] == (0, 0, 0) and len(set(control["Q_i"])) == 3
    assert control["actual_solution"] is control["branch_health_claim"] is False


def test_determinant_zero_spatial_sum_and_disconnected_fields_are_distinct():
    singular = potential.evaluate((1, -1), (1, 1), (1, 2))
    assert singular["S_a"] == 0 and singular["regular_sum"] is False
    assert singular["rho"] == singular["pressure"] == (0, 0)
    disconnected = potential.evaluate((1, 0, 1), (1, 2, 3), (1, 5, 2))
    assert disconnected["active"] == (0, 2)
    assert disconnected["rho"][1] == disconnected["pressure"][1] == 0
    assert background.disconnected_GR(3, 2)["Hprime"] == -sp.Rational(1, 3)


def test_determinant_auxiliary_elimination_does_not_move_the_physical_source():
    actual = controls.auxiliary_map()
    assert actual["B"] == -3*actual["lambda"]/2
    assert actual["p_over_beta"] == actual["lambda"]/2
    assert actual["reduced_action"] == actual["expected_reduced_action"]
    assert actual["matter_stays_on_original_EH_leaves"] is True
    assert actual["matter_on_auxiliary_w"] is False
    assert controls.source_relocation_control() == (-sp.Rational(1, 2), sp.Rational(1, 2),
                                                    sp.Rational(1, 2), sp.Rational(1, 2))


def test_determinant_unweighted_source_and_inverse_ratio_K_are_wrong():
    result = background.reconstruct([1, 2], [1, 3], [1, 2], 0, [1, 2])
    assert result["K"] == 19 and result["K"] != 1+sp.Rational(2, 9)
    assert result["weighted_null"] == 325 and result["weighted_null"] != 3
    assert result["Hprime"] < 0


@pytest.mark.parametrize("bad", [True, 0.5, sp.Float("0.5"), sp.oo, sp.nan, sp.sqrt(2)])
def test_determinant_numeric_interfaces_reject_inexact_nonfinite_or_nonrational_inputs(bad):
    with pytest.raises((TypeError, ValueError)):
        potential.evaluate([bad], [1], [1])
    with pytest.raises((TypeError, ValueError)):
        background.reconstruct([bad], [1], [1], 1, [1])


@pytest.mark.parametrize("bad", [0, -1])
def test_determinant_positive_kinetic_coframe_and_lambda_domains(bad):
    with pytest.raises(ValueError):
        potential.evaluate([1], [bad], [1])
    with pytest.raises(ValueError):
        potential.evaluate([1], [1], [bad])
    with pytest.raises(ValueError):
        potential.evaluate([1], [1], [1], lambda_=bad)
    with pytest.raises(ValueError):
        background.reconstruct([bad], [1], [1], 0, [1])


def test_determinant_normalization_source_and_component_guards():
    for args in (([1], [2], [1], 0, [1]), ([1], [1], [2], 0, [1]),
                 ([1], [1], [1], 0, [-1]), ([], [], [], 0, []), ([1, 2], [1], [1], 0, [1])):
        with pytest.raises(ValueError):
            background.reconstruct(*args)
    with pytest.raises(TypeError):
        background.reconstruct([1], [1], [1], 0, [1], physical=True)
    with pytest.raises(ValueError):
        background.reconstruct([1], [1], [1], 0, [1], physical=1)
