"""Independent polynomial, Gaussian, integral and exact-input regressions."""

import pytest
import sympy as sp
from p8_polynomial_vacuum import (
    audit,
    calibration,
    gaussian,
    model,
    potential,
    powercount,
)


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_identity(name, value):
    if isinstance(value, sp.MatrixBase):
        assert all(v == 0 for v in value), name
    else:
        assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_bound_or_written_proof_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_unsupported_input_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_and_scope_control(name, value):
    assert bool(value), name


def test_independent_rational_parameter_calibration():
    d = calibration.point(2, 1)
    assert d["D"] == 4 and d["heavy_mass_squared"] == 6
    assert d["cubic_squared"] == 256 and d["quartic"] == 160
    assert d["stable_margin"] == 32 and d["kernel_upper"] == sp.Rational(176, 3)
    assert d["tree_forward_b2"] == 8


def test_selected_fixed_parameters_and_complete_square_margin():
    d = calibration.point()
    m = model.data()["actual_parameters"]
    assert d["lambda"] == sp.Rational(1, 10**600)
    assert d["gamma"] == sp.Rational(1024, 10**800)
    assert d["cubic_squared"] == m["cubic_coupling_squared"]
    assert d["stable_margin"] == m["positive_completed_square_quartic_margin"] > 0


@pytest.mark.parametrize("y", (0, 1, 10**197))
def test_exact_Euclidean_kernel_calibrations(y):
    p = calibration.kernel_point(y, -sp.Rational(7, 3))
    d = calibration.point()
    assert d["stable_margin"] / 2 <= p["F"] < d["kernel_upper"]
    assert p["Schur_kernel"] > y + 1
    assert p["log_increment"] > 0


def test_zero_field_remainder_and_log_increment_zero():
    p = calibration.kernel_point(1, 0)
    assert p["log_increment"] == p["rational_one_loop_remainder_upper"] == 0
    assert p["Schur_kernel"] == 2


def test_arbitrarily_many_cubic_vertices_are_not_enumeration_cutoff():
    assert powercount.degree(0, 0, 1000) == -996
    assert powercount.degree(2, 2, 1000) == -1000


def test_no_divergent_heavy_kinetic_or_higher_heavy_monomial():
    rows = powercount.candidates()
    assert {(a, b) for a, b, _ in rows if b} == {(0, 1), (0, 2), (2, 1)}
    assert all(b <= 2 and not (a >= 2 and b >= 2) for a, b, _ in rows)
    assert max(v for (a, b, _), v in rows.items() if b == 2) == 0


def test_literal_nonlocal_quartic_has_mixed_site_terms():
    d = gaussian.data()
    phi1, phi2 = sp.symbols("phi1 phi2", real=True)
    coupling = sp.Symbol("cubic_G", real=True)
    mixed = sp.diff(d["effective_quartic_action"], phi1, 2, phi2, 2)
    assert mixed != 0
    assert mixed.subs(coupling, 0) == 0


def test_one_loop_kernel_uses_same_actual_parameter_calibration():
    d = calibration.point()
    p = potential.data()
    assert p["actual_kernel_upper"] == d["kernel_upper"]
    assert p["actual_rational_Phi6_upper_coefficient"] == d["kernel_upper"] ** 3 / 1728
    assert 0 < p["actual_rational_Phi6_upper_coefficient"] < sp.Rational(1, 10**618)


def test_counts_and_original_scope():
    assert len(audit.residuals()) == 37
    assert (
        sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in audit.residuals().values()
        )
        == 51
    )
    assert len(audit.gates()) == 39
    assert len(audit.controls()) == 14
    assert audit.rejected_inputs() == 87


def test_mass_weighted_one_loop_bound_uses_actual_heavy_scale():
    d = calibration.point()
    p = potential.data()
    expected = (
        (d["stable_margin"] / 2) ** 3
        + (d["cubic_squared"] / d["heavy_mass_squared"]) ** 3 / d["heavy_mass_squared"]
    ) / 432
    assert p["actual_mass_weighted_rational_Phi6_upper_coefficient"] == expected
    assert 0 < expected < sp.Rational(1, 10**815)
