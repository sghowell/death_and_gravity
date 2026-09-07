from fractions import Fraction as Q

import pytest
import sympy as sp
from p8_composite_response import bounds, independent, model, response


def test_full_background_action_and_flux_identities():
    assert all(sp.simplify(value) == 0 for value in model.checks().values())


def test_physical_limiting_normalization_keeps_its_derivatives():
    assert all(sp.simplify(value) == 0 for value in response.checks().values())


def test_physical_probe_source_units_and_rescaling():
    assert all(value == 0 for value in response.source_checks().values())


def test_separate_fraction_polynomials_and_comparison():
    assert len(independent.checks()["coefficientwise_identities"]) == 3
    margins = response.rational_checks()
    assert margins["relative_kernel_coefficient_lower"] == "41/304"
    assert margins["finite_full_minus_locked_lower"] == "1/75"


def test_finite_family_initial_offset_is_not_dropped():
    eps = sp.Symbol("epsilon", nonnegative=True)
    root = model.initial_root(eps)
    assert root.subs(eps, 0) == 2
    assert sp.diff(root, eps).subs(eps, 0) == sp.Rational(9, 4)
    assert sp.expand(root**2-4) == 3*eps**3+9*eps**2+9*eps
    # R0=2 would instead prescribe eta0=9/(1+epsilon)^3.
    assert sp.diff(9/(1+eps)**3, eps).subs(eps, 0) == -27


def test_explicit_finite_parameter_and_amplitude_gates():
    r = bounds.build()
    assert r["epsilon_range"]["negative_binary_exponent"] == 187116
    assert r["sigma_amplitude_range"]["negative_binary_exponent"] == 67162
    assert r["background_state_Lipschitz"] == 26372
    assert r["phase_matrix_norm"] == 33572
    assert all(Q(value) > 0 for value in r["rational_margins"].values())


def test_every_required_box_denominator_is_positive():
    r = bounds.build()
    for key in r["positive_box_coefficients"]:
        assert Q(r["records"][key]["value"]["lower"]) > 0
    assert Q(r["records"]["e"]["value"]["upper"]) < Q(1, 100)


def test_prepared_same_metric_initial_data_and_locked_flux():
    d = model.derive()
    eps = d["epsilon"]
    value, velocity = sp.symbols("value velocity", real=True)
    data = sp.Matrix([eps*value, eps*d["A_g"]*velocity,
                      value, d["A_f"]*velocity, value, d["C_lock"]*velocity])
    right = sp.Matrix(d["matrix"])*data
    assert sp.cancel(right[0]-eps*velocity) == 0
    assert sp.cancel(right[2]-velocity) == 0
    assert sp.cancel(right[4]-velocity) == 0
    assert data[0]-eps*data[2] == 0
    assert sp.cancel(right[0]-eps*right[2]) == 0


def test_large_band_does_not_inherit_positive_kernel_lower_bound():
    assert Q(1, 4)-Q(140, 19)*Q(1, 4) < 0


def test_opposite_probe_sign_does_not_inherit_positive_response_claim():
    assert -Q(1, 60) < 0


def test_zero_probe_zero_retarded_data_has_no_forced_response():
    d = model.derive()
    assert sp.Matrix(d["matrix"])*sp.zeros(6, 1) == sp.zeros(6, 1)


def test_gate_logarithms_are_integer_majorants():
    for value in (Q(0), Q(1, 2), Q(1), Q(3, 2), Q(1450), Q(68152776910250)):
        assert value <= 2**bounds.log2_majorant(value)
    with pytest.raises(ValueError):
        bounds.log2_majorant(-1)
    for value in (True, 0.5, float("nan")):
        with pytest.raises(TypeError):
            bounds.log2_majorant(value)
