"""Independent curved pole calculations, including failed-reference control."""
import pytest
import sympy as sp
from p8_vector_quadratic import (
    bimetric,
    checks,
    geometry,
    invariants,
    kernel,
    response,
    tensors,
    verify,
)


def all_zero(values):
    for value in values.values():
        if isinstance(value, sp.MatrixBase):
            assert value == sp.zeros(value.rows, value.cols)
        else:
            assert value == 0


def test_ordered_frame_derivatives_and_curvature_commutators():
    all_zero(tensors.checks())


def test_independent_coordinate_curvatures():
    all_zero(geometry.checks())


@pytest.mark.parametrize("order", (0, 2, 4))
def test_frozen_flat_Feynman_integral(order):
    assert checks.flat()["flat_frozen_Feynman_pole_order_"+str(order)] == 0


def test_curved_scalar_mass_conformal_Proca():
    all_zero(checks.scalar_conformal()["compact_variation_checks"])


def test_independent_two_derivative_tensor_basis():
    all_zero(checks.general_second())


def test_constant_anisotropic_auxiliary_metric():
    all_zero(checks.constant_anisotropic())


def test_literal_reference_is_withheld_with_nonzero_control():
    data = checks.withheld_control()
    assert data["nonzero_exact_fixture"] == sp.Rational(237, 20)
    assert data["nonzero_compact_variation_residual"] != 0
    assert kernel.pole(4) != invariants.evaluate(4)/(960*tensors.mass**4)


@pytest.mark.parametrize("order", (0, 2, 4))
def test_actual_lapse_chain_rule_and_compact_local_action(order):
    all_zero(response.checks(order))


def test_actual_mass_profiles_replay_frozen_first_jets():
    all_zero(response.actual_mass_checks())


def test_full_basis_and_quadratic_clock_boundary():
    assert (len(invariants.ZERO), len(invariants.SECOND), len(invariants.FOURTH)) == (2, 9, 48)
    assert all(value is True for value in verify.proof_checks().values())


def test_all_certified_arithmetic_is_exact():
    values = [bimetric.local_quadratic_density()]+[kernel.pole(order) for order in (0, 2, 4)]
    assert all(not value.atoms(sp.Float) for value in values)


def test_warm_cache_does_not_accept_nested_boolean_index():
    tensors.plus.covariant((0,), 0, 0)
    with pytest.raises(ValueError):
        tensors.plus.covariant((False,), 0, 0)
    geometry.connection(0, 0, 0)
    with pytest.raises(ValueError):
        geometry.connection(False, 0, 0)


def test_unsupported_inputs_and_scope_controls():
    assert verify.controls()["rejected_inputs"] == 170
