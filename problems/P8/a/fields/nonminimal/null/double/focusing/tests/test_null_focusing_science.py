"""Actual index, full costs, complete controls and scalar past replay."""

from fractions import Fraction

import pytest
import sympy as sp
from p8a_null_focusing import audit, control, costs, thermal, verify

ROWS = verify.residuals()


@pytest.mark.parametrize("name", list(ROWS))
def test_native_exact_identity(name):
    assert sp.simplify(ROWS[name]) == 0


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_invalid_null_budget_fails(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_independent_fraction_full_costs_and_continuous_control():
    assert verify.checked_constants()["independent_Fraction_reconstruction_agrees"]
    assert len(ROWS) == 35
    assert len(audit.proof_checks()) == 28
    assert all(map(bool, audit.proof_checks().values()))
    assert audit.rejected_inputs() == 24


@pytest.mark.parametrize("d", [0, Fraction(1, 10**12)])
@pytest.mark.parametrize("z", [0, Fraction(1, 10**6)])
@pytest.mark.parametrize("s", [0, Fraction(1, 10)])
def test_actual_whole_rectangular_gate_margin(d, z, s):
    point = audit.gate(d, z, s)
    assert point["strict_outgoing_index_margin"] > sp.Rational(1, 2)
    assert point["future_null_affine_endpoint_factor"] == 2
    assert point["future_caps_are_conditional_not_inferred_from_past"]
    assert not point["state_homogeneity_or_single_null_QEI_assumed"]


def test_complete_geometry_control_really_passes_all_continuous_caps():
    data = control.data()
    assert data["future_A_lower"] > sp.Rational(1, 2)
    assert data["past_absolute_A_minus_one_bound"] < sp.Rational(1, 10)
    assert min(data["past_derivative_margins"].values()) > 0
    assert min(data["future_derivative_margins"].values()) > 0


def test_finite_spatial_width_and_past_Ricci_cost_not_discarded():
    data = costs.data()
    assert data["past_Ricci_upper_cost"] == sp.Rational(104, 7875) > 0
    for name in ("past", "future"):
        row = data[name]
        assert row["quantum_cost"] > sp.Rational(5, 9) * row["P_norm_majorant"]
        assert row["state_cost"] > 0


def test_actual_scalar_thermal_past_not_future_witness():
    data = thermal.data()
    assert data["strict_null_field_budget_squared_margin"] > 0
    assert min(data["strict_affine_derivative_margins"]) > 0
    assert data["actual_parent_scalar_history"]["actual_scalar_state_and_SEE_past"]
    assert not data["actual_parent_scalar_history"][
        "future_state_cap_proved_on_every_shorter_segment"
    ]


@pytest.mark.parametrize("bad", [0.1, sp.Float("0.1"), sp.oo, sp.nan])
def test_inexact_or_nonfinite_serialization_fails(bad):
    with pytest.raises((TypeError, ValueError)):
        verify.serialize({"bad": bad})
