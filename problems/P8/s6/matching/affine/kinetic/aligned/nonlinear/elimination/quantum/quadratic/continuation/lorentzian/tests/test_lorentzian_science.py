"""Independent Lorentzian geometry, poles and physical-time projection."""
import pytest
import sympy as sp
from p8_vector_lorentzian_pole import checks, frame, response, verify, wick
from p8_vector_quadratic import tensors


def zeros(values):
    for value in values.values():
        assert value == (sp.zeros(value.rows, value.cols) if isinstance(value, sp.MatrixBase) else 0)


def test_physical_frame_curvature_and_commutators():
    zeros(checks.frame_checks())


def test_direct_signature_agrees_in_arbitrary_dimension():
    zeros(checks.signature_checks())


@pytest.mark.parametrize("order", (0, 2, 4))
def test_frozen_Feynman_pole_in_Lorentzian_signature(order):
    assert checks.flat()["Lorentzian_flat_Feynman_pole_"+str(order)] == 0


def test_original_physical_curvatures_replayed():
    zeros(checks.physical_geometry())


@pytest.mark.parametrize("order", (2, 4))
def test_actual_lapse_projection_grading_and_measure(order):
    zeros(response.checks(order))


@pytest.mark.parametrize("order", (2, 4))
def test_physical_bare_counterterm_sign(order):
    assert response.finite_counterterm(order) == 2*response.operator(order, 1)


def test_local_conversion_does_not_rotate_the_background_history():
    data = checks.history_control()
    assert data["nonzero_history_difference"] == -4*response.u**2
    assert data["original_half_time_scale_factor"] == sp.Rational(25, 16)
    assert data["changed_half_time_scale_factor"] == sp.Rational(9, 16)
    assert data["fourth_pole_second_time_square_bounce"] == sp.Rational(128, 6561)


def test_warm_cache_native_indices_and_jets():
    frame.covariant(tensors.plus, (0,), 0, 0)
    with pytest.raises(ValueError):
        frame.covariant(tensors.plus, (False,), 0, 0)
    wick.loop_pole(2, 0)
    with pytest.raises(ValueError):
        wick.loop_pole(2, False)


def test_exact_arithmetic_only():
    assert all(not wick.loop_pole(order, jet).atoms(sp.Float) for order in (2, 4) for jet in (0, 1))


def test_unchanged_clock_and_symbolic_dimension_checks():
    assert all(value is True for value in verify.proof_checks().values())


def test_unsupported_inputs_rejected():
    assert verify.controls()["rejected_inputs"] == 204
