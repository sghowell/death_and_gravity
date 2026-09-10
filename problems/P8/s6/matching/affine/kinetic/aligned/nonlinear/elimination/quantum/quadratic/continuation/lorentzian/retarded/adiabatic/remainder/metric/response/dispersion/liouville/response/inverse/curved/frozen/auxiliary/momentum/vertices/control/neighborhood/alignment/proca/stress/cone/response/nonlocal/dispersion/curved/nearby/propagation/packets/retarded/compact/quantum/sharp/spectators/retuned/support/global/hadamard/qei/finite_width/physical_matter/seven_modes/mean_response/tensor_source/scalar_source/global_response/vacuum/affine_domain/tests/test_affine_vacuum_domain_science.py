"""Independent affine action, domain, exact-input and vacuum regressions."""

import pytest
import sympy as sp
from p8_affine_vacuum_domain import (
    audit,
    bounds,
    calibration,
    covariance,
    family,
    local,
    lower,
)


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_independent_identity(name, value):
    if isinstance(value, sp.MatrixBase):
        assert all(v == 0 for v in value), name
    else:
        assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_exact_domain_bound_and_written_proof_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args",
    calibration.bad_cases(),
    ids=[row[0] for row in calibration.bad_cases()],
)
def test_unsupported_exact_calibration_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_obstruction_and_scope_control(name, value):
    assert bool(value), name


def test_old_polynomial_lift_requires_two_distinct_rank_loci():
    d = family.old_obstruction()
    assert d["old_R_at_zero"] > sp.Rational(1, 2) > d["old_R_at_quarter"]
    assert d["old_R_at_quarter"] < sp.Rational(1, 2) < d["old_R_at_half"]


def test_selected_regular_domain_margin_and_canonical_scale():
    d = calibration.point()
    assert d["kappa"] == 10**800
    assert d["tensor_factor_lower"] == sp.Rational(4097, 8192)
    assert d["additional_quotient_factor_lower"] == sp.Rational(1, 4096)
    assert calibration.point(1026, 0)["kappa"] == sp.Rational(1026, 1024) * 10**800


@pytest.mark.parametrize(
    "X", (-sp.Rational(1, 8192), 0, sp.Rational(1, 4), 1, sp.Rational(11, 10))
)
def test_exact_connected_domain_calibrations(X):
    assert calibration.point(1024, X)["physical_X"] == X


def test_general_null_hessian_not_zero_gradient_hessian():
    a = covariance.hessian((1, 1, 0, 0), sp.Rational(1, 2), 0, 512)
    b = covariance.hessian((0, 0, 0, 0), sp.Rational(1, 2), 0, 512)
    assert a != b
    assert a.shape == b.shape == (64, 64)


def test_actual_affine_removable_values_at_vacuum():
    d = local.data()
    u = family.u
    assert d["source_c_over_X_at_vacuum"].subs(u, 0) == -512
    assert d["source_F4_at_vacuum"].subs(u, 0) == 512
    assert d["source_J3_over_X_at_vacuum"].subs(u, 0) == 0
    assert d["regular_auxiliary_q_order"] == 4


def test_W_is_strict_away_from_R_one():
    r = sp.Symbol("positive_R_variable", positive=True)
    W = lower.data()["boundary_obstruction_nonnegative_function"]
    assert W.subs(r, sp.Rational(1, 16)) == sp.Rational(22, 3)
    assert W.subs(r, 16) == sp.Rational(17, 6)
    assert W.subs(r, 1) == 0


def test_exact_switch_and_full_lower_norms_are_finite_rationals():
    rows = bounds.tube()["weighted_four_jet_error_bounds"]
    assert set(rows) == {"F", "F2", "A3", "A4", "A5"}
    assert all(
        isinstance(v, sp.Rational) and 0 < v < sp.Rational(1, 10**400)
        for v in rows.values()
    )


def test_native_counts_and_original_scope():
    assert len(audit.residuals()) == 70
    assert (
        sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in audit.residuals().values()
        )
        == 4930
    )
    assert len(audit.gates()) == 51
    assert len(audit.controls()) == 13
    assert audit.rejected_inputs() == 54
