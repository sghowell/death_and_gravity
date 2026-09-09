"""Native scalar identities, independent constants and genuine bad-input controls."""

from fractions import Fraction

import pytest
import sympy as sp
from p8a_nonminimal import audit, cosmology, field, independent, thermal, verify

ROWS = [
    (group, name, value)
    for group, rows in audit.exact_groups().items()
    for name, value in rows.items()
]


@pytest.mark.parametrize("group,name,value", ROWS, ids=[name for _, name, _ in ROWS])
def test_native_exact_identity(group, name, value):
    entries = list(value) if isinstance(value, sp.MatrixBase) else [value]
    assert all(sp.simplify(v) == 0 for v in entries), (group, name)


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_bad_domain_is_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_independent_Fraction_whole_cost_and_actual_state_replay():
    result = verify.checked_constants()
    assert result["independent_Fraction_reconstruction_agrees"]
    assert len(audit.proof_checks()) == 23
    assert all(map(bool, audit.proof_checks().values()))
    assert audit.rejected_inputs() == 55
    assert sum(len(v) if isinstance(v, sp.MatrixBase) else 1 for _, _, v in ROWS) == 72


@pytest.mark.parametrize("xi", [0, Fraction(1, 12), Fraction(1, 6), Fraction(1, 4)])
def test_full_allowed_flat_interval_has_independent_spectral_normalization(xi):
    actual = audit.flat_parameters(xi)
    assert actual["quantum_coefficient_times_pi_squared_over_hbar"] == sp.Rational(
        independent.spectral_moment(xi)
    )
    assert bool(actual["curved_cosmological_transport_proved_at_this_coupling"]) == (
        xi == Fraction(1, 6)
    )


@pytest.mark.parametrize(
    "delta", [Fraction(1, 10**8), Fraction(1, 10**10), Fraction(1, 10**12)]
)
def test_actual_thermal_state_history_continuous_proof_gates_at_exact_inputs(delta):
    result = thermal.history(delta)
    assert min(result["strict_C3_margins"]) > 0
    assert result["strict_squared_field_gate_margin"] > 0
    assert result["lambda"] == sp.Rational(delta) / 360
    assert not result["future_state_cap_proved_on_every_shorter_segment"]
    point = thermal.branch_point(2, delta)
    assert point["a_fourth"] == 1
    assert point["zeta_squared"] == result["initial_and_past_zeta_squared_upper"]
    assert point["kappa_EED_times_tau_squared"] > 0


def test_state_penalty_and_finite_beta_are_not_deleted():
    worst = cosmology.gate(cosmology.DELTA_MAX, cosmology.ZETA_MAX, cosmology.SIGMA_MAX)
    no_field = cosmology.gate(cosmology.DELTA_MAX, 0, cosmology.SIGMA_MAX)
    assert (
        no_field["strict_focusing_margin"] - worst["strict_focusing_margin"]
        == cosmology.costs()["Wick_square_cost"] / 5000
    )
    assert cosmology.costs()["scalar_Cbeta"] > 0
    assert worst["strict_focusing_margin"] == sp.Rational(63325013, 350000000)


def test_coherent_control_trace_and_local_energy_are_actual_stress():
    data = field.data()
    assert data["checks"]["coherent_conformal_energy_is_EED"] == 0
    assert data["checks"]["coherent_local_example_solves_massless_equation"] == 0
    assert "amplitude**2" in str(data["coherent_line_energy"])


@pytest.mark.parametrize("value", [0.1, sp.Float("0.1"), sp.oo, sp.nan])
def test_inexact_nonfinite_certificate_entries_fail(value):
    with pytest.raises((ValueError, TypeError)):
        verify.serialize({"bad": value})
