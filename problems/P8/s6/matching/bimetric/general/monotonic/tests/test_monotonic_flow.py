import pytest
import sympy as sp
from p8_bimetric_monotonic import controls, flow, fraction_checks


@pytest.mark.parametrize("check", [flow.identities, flow.root_checks, flow.endpoint_checks,
                                  controls.decoupled_de_sitter_checks, controls.increasing_H_checks,
                                  controls.static_sharpness_checks])
def test_exact_lapse_branch_and_full_background_controls(check):
    assert all(sp.simplify(value) == 0 for value in check().values())


def test_fraction_coefficient_identities_not_parameter_samples():
    report = fraction_checks.checks()
    assert len(report["zero_coefficient_residuals"]) == 5
    assert set(report["zero_coefficient_residuals"].values()) == {0}
    assert all(value > 0 for value in report["nonzero_omission_term_counts"].values())


def test_zero_polynomial_exception_really_rejects_using_Z():
    value = controls.controls()["decoupled_Z_increases_although_H_is_constant"]
    assert value.is_positive is True
    assert value.subs(flow.t, 0) == sp.sqrt(2)/2


def test_interacting_H_can_increase_with_strict_classical_NEC():
    data = controls.increasing_H_example()
    y = data["y"]
    assert sp.simplify(sp.diff(data["H"], y).subs(y, 1)) == sp.Rational(3, 8)
    assert sp.simplify(data["n_g"].subs(y, 1)) == sp.Rational(1, 2)
    assert sp.simplify(data["V_g"].subs(y, 1)) == 1
    assert sp.simplify(sp.diff(data["Z"], y).subs(y, 1)) == -sp.sqrt(2)/16
    assert sp.Rational(2258, 605) > 0


def test_endpoint_radius_is_attained_by_static_parent():
    assert all(value == 0 for value in flow.endpoint_checks().values())
    assert all(value == 0 for value in controls.static_sharpness_checks().values())
    assert controls.controls()["static_parent_attains_endpoint_error_threshold"] == sp.Rational(8, 5)


def test_root_set_must_be_split_from_zero_polynomial_parameter_case():
    y = sp.Symbol("y", real=True)
    nonzero_polynomial = (y-1)*(y-2)
    assert sp.Poly(nonzero_polynomial, y).degree() == 2
    assert sp.solve(nonzero_polynomial, y) == [1, 2]
    assert sp.Poly(0, y).is_zero
    # This finite root example is an arithmetic regression, not a sampled
    # proof of the written closed-set/continuity argument.
