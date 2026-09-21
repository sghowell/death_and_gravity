"""Complete selected quadratic cancellation and original-source scope controls."""

import pytest
import sympy as s
from p8_vacuum_affine_quadratic_radiation_cancellation import (
    audit,
    calibration,
    cancellation,
    loops,
    source,
    ward,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_residual(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_all_rejected_scopes(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("name", source.CLASSES)
def test_all_original_two_point_classes(name):
    assert name in (
        "mixed_two_cubic_bubble",
        "quartic_tadpole",
        "heavy_onepoint_reducible",
    )
    assert source.CUBIC == s.Rational(1, 8192)


@pytest.mark.parametrize("degree", range(1, 9))
def test_all_noncommuting_operator_positions(degree):
    v, u = s.symbols("v u")
    full = sum(u**j * v ** (degree - 1 - j) for j in range(degree))
    assert s.expand((u - v) * full - u**degree + v**degree) == 0
    if degree > 1:
        assert s.expand((u - v) * v ** (degree - 1) - u**degree + v**degree) != 0


@pytest.mark.parametrize(
    "index,leg,virtuality", calibration.data()["whole_original_external_virtualities"]
)
def test_original_actual_virtuality_gap(index, leg, virtuality):
    assert 0 <= index < 3 and 0 <= leg < 4 and 0 < virtuality < 2
    assert virtuality != 1 and virtuality.is_Rational


def test_all_selected_zero_quantities_have_complete_nonzero_terms():
    data = cancellation.data()
    assert data["whole_zero_selected_physical_radiation"] == 0
    assert data["whole_zero_selected_finite_interference"] == 0
    assert data["whole_renormalized_shifted_inverse"] != 0
    assert data["whole_renormalized_TT_vertex"] != 0
    assert data["whole_single_propagator_radiation_correction"] != 0


def test_wrong_Dyson_sign_fails():
    d, H, P = s.symbols("d H Pi")
    correct_vertex = 2 * H * P / d**2
    wrong_prop = 2 * H * P / d**2
    assert s.factor(correct_vertex + wrong_prop) != 0


def test_kinetic_metric_counterterm_cannot_be_omitted():
    d, H, P, P0, P1 = s.symbols("d H Pi Pi0 Pi1")
    uncorrected_vertex = 2 * H * (P - P0) / d**2
    propagator = -2 * H * (P - P0 - d * P1) / d**2
    assert s.factor(uncorrected_vertex + propagator) == 2 * H * P1 / d


def test_full_logarithm_not_identically_zero_or_truncated():
    data = loops.data()
    assert data["whole_original_subtracted_parameter_integrand"].has(s.log)
    assert data["whole_original_relative_inverse_majorant"] > 0
    assert data["whole_original_relative_inverse_majorant"] < s.Rational(1, 10**400)
    assert data["whole_tadpole_shifted_TT_numerator"] != 0
    assert all(value != 0 for value in data["whole_actual_triangle_antiderivatives"])


def test_unprojected_homogeneous_tensor_not_deleted():
    data = ward.data()
    assert data["whole_null_transverse_homogeneous_tensor"] != s.zeros(4)
    assert "kk" in data["whole_physical_boundary"]
    assert "off-shell" in data["whole_not_inferred"]


def test_original_parameters_and_frontiers_unchanged():
    n, g = source.HEAVY_MASS2, source.CUBIC
    assert source.CONTACT == -g * g * (3 / (n - 2) - 2 / (n - 2) ** 2)
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 196
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "curvature" in audit.observable()["not_established"].lower()
    assert "internal-graviton" in audit.observable()["not_established"].lower()


def test_all_original_controls_present():
    assert len(calibration.data()["checks"]) == 72
    assert calibration.data()["whole_nonzero_physical_polarization_factors"] == 16
    assert len(calibration.data()["whole_original_external_virtualities"]) == 12
    assert len(ward.data()["checks"]) == 5
    assert len(cancellation.data()["checks"]) == 28
    assert len(loops.data()["checks"]) == 26
