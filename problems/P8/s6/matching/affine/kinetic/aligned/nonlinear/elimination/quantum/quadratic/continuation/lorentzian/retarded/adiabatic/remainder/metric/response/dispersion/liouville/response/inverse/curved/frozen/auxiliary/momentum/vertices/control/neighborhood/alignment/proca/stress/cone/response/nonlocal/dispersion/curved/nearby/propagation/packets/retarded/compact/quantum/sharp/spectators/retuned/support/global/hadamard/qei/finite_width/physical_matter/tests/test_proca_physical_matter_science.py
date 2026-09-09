"""Actual lapse, geometric, fourteen-channel and independent canonical checks."""

import pytest
import sympy as sp
from p8_proca_physical_matter import (
    audit,
    bridge,
    center,
    coefficients,
    controls,
    observable,
    spatial,
    verify,
)

ROWS = audit.residuals()


@pytest.mark.parametrize("name", list(ROWS))
def test_native_exact_identity(name):
    val = ROWS[name]
    assert all(
        v == 0 for v in (list(val) if isinstance(val, sp.MatrixBase) else [val])
    ), name


@pytest.mark.parametrize("name", list(audit.gates()))
def test_written_and_exact_proof_gate(name):
    assert bool(audit.gates()[name]), name


@pytest.mark.parametrize(
    "name,call,args", controls.bad_cases(), ids=[row[0] for row in controls.bad_cases()]
)
def test_invalid_actual_observable_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_actual_manifest_and_exact_counts():
    assert len(ROWS) == 51
    assert (
        sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in ROWS.values()
        )
        == 1097
    )
    assert len(audit.gates()) == 19
    assert controls.rejected_inputs() == 29
    assert len(verify.source_files()) == 20


@pytest.mark.parametrize("key", ("rho", "pressure", "missing_if_lapse_truncated"))
def test_all_independent_fourteen_channel_Hessian_entries(key):
    assert spatial.matrices()[key] == center.data()["matrices"][key]
    assert spatial.matrices()[key].shape == (14, 14)


@pytest.mark.parametrize("point", (0, sp.Rational(-1, 100), sp.Rational(1, 100)))
def test_nonzero_output_mixed_control_does_not_erase_second_lapse(point):
    data = spatial.nonzero_output_control(point)
    assert all(data["checks"].values())
    assert set(data["nonzero_output_coefficients"]) == {3, 5, 6}
    assert all(
        row["missing_if_lapse_truncated"] != 0
        for row in data["nonzero_output_coefficients"].values()
    )


@pytest.mark.parametrize("bad", (0.1, sp.Float("0.1"), sp.oo, sp.nan))
def test_inexact_or_nonfinite_report_values_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        verify.serialize({"bad": bad})


def test_principal_density_not_confused_with_positive_evolution_Hamiltonian():
    c = center.data()
    assert c["negative_principal_configuration_test"] == -sp.Rational(134, 135)
    assert c["matrices"]["rho"][5, 5] == -sp.Rational(2, 135)
    assert c["matrices"]["rho"][12, 12] == -sp.Rational(200000, 81)
    assert c["matrices"]["missing_if_lapse_truncated"] != sp.zeros(14)


def test_full_positive_classical_density_not_refuted_by_quadratic_correction():
    d = controls.data()
    assert d["pure_tensor_implicit_branch_lapse_Hessian"] == -sp.Rational(243, 80)
    expression = d["actual_full_classical_center_density"]
    N = next(s for s in expression.free_symbols if str(s) == "positive_physical_lapse")
    dc = next(
        s for s in expression.free_symbols if str(s) == "matter_density_deviation"
    )
    G = next(
        s
        for s in expression.free_symbols
        if str(s) == "nonnegative_spatial_gradient_squared"
    )
    for point in (sp.Rational(1, 2), 1, 2):
        assert expression.subs({N: point, dc: 0, G: 0}) > 0


def test_both_boundary_shifts_and_scaled_symplectic_factor_retained():
    d = bridge.data()
    assert d["missing_matter_boundary_shift_matrix"] != sp.zeros(4)
    assert d["missing_metric_boundary_shift_matrix"] != sp.zeros(4)
    assert (
        d["actual_parent_packet_to_fixed_spatial_phase"]
        == center.data()["center_packet_map"]
    )


def test_formal_Ward_order_accounting_is_not_independent_conservation():
    d = controls.data()
    assert len(sp.Add.make_args(d["second_order_divergence_terms"])) == 3
    assert controls.summary()["no_automatic_positive_square_QEI_or_full_Ward_transfer"]


def test_actual_action_lapse_jets_and_physical_normalization_not_reset():
    d = coefficients.data()
    assert d["added_background_row"][0] == d["added_background_row"][1] == 0
    assert d["added_background_row"][2] != 0
    assert observable.data()["physical_stress_multiplier"] != 1
