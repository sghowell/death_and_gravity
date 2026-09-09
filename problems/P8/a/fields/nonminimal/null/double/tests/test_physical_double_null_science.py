"""Native physical identities, exact constants and actual domain controls."""

from fractions import Fraction

import pytest
import sympy as sp
from p8a_double_null import audit, curved, verify

ROWS = verify.residuals()


@pytest.mark.parametrize("name", list(ROWS))
def test_native_exact_identity(name):
    value = ROWS[name]
    entries = list(value) if isinstance(value, sp.MatrixBase) else [value]
    assert all(sp.simplify(v) == 0 for v in entries)


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_invalid_sampling_data_fails(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_independent_fraction_polynomial_and_full_constant_reconstruction():
    assert verify.checked_constants()["independent_Fraction_reconstruction_agrees"]
    assert len(ROWS) == 39
    assert verify.exact_checks()[1] == 51
    assert len(audit.proof_checks()) == 18
    assert all(map(bool, audit.proof_checks().values()))
    assert audit.rejected_inputs() == 39


@pytest.mark.parametrize("xi", [0, Fraction(1, 6), Fraction(1, 2)])
@pytest.mark.parametrize(
    "widths", [(1, 1), (Fraction(1, 10), 2), (2, Fraction(1, 100))]
)
def test_valid_exact_positive_widths_and_full_coupling_interval(xi, widths):
    point = audit.product_budget(*widths, xi, 3)
    assert point["quantum_cost_coefficient_of_hbar"] > 0
    assert point["state_cost"] >= 0
    assert point["boost_squared"] > 0
    assert not point["single_null_limit_taken"]
    assert not point["curved_or_singularity_bound_assigned"]


def test_both_widths_and_state_cost_have_exact_scaling():
    point = audit.product_budget(2, 3, Fraction(1, 6), 5)
    assert point["quantum_cost_coefficient_of_hbar"] == sp.Rational(14117, 1536) / (
        sp.pi**2 * 2**3 * 3
    )
    assert point["state_cost"] == 5
    other = audit.product_budget(4, 6, Fraction(1, 6), 5)
    assert (
        point["quantum_cost_coefficient_of_hbar"]
        == 16 * other["quantum_cost_coefficient_of_hbar"]
    )
    assert point["state_cost"] == 4 * other["state_cost"]


def test_single_null_shrinking_increases_quantum_cost_without_new_cutoff():
    first = audit.product_budget(1, 1, Fraction(1, 6), 0)
    next_point = audit.product_budget(1, Fraction(1, 100), Fraction(1, 6), 0)
    assert (
        next_point["quantum_cost_coefficient_of_hbar"]
        == 100 * first["quantum_cost_coefficient_of_hbar"]
    )
    assert next_point["state_cost"] == 0


def test_actual_curved_affine_normalization_and_reference_not_omitted():
    data = curved.data()
    assert str(data["physical_reference_null_numerator"]).count("beta_S") > 0
    assert data["physical_plane_volume"] == "a*dt*dz=a^2*deta*dz"


@pytest.mark.parametrize("bad", [0.1, sp.Float("0.1"), sp.oo, sp.nan])
def test_inexact_or_nonfinite_serialization_fails(bad):
    with pytest.raises((ValueError, TypeError)):
        verify.serialize({"bad": bad})
