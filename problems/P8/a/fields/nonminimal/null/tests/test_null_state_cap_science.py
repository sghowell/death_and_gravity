"""Actual-state identities, complete constants and genuine parameter controls."""

from fractions import Fraction

import pytest
import sympy as sp
from p8a_null_cap import audit, construction, normalization, verify

ROWS = {**construction.checks(), **normalization.data()["checks"]}


@pytest.mark.parametrize("name", list(ROWS))
def test_native_exact_identity(name):
    assert sp.simplify(ROWS[name]) == 0


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_invalid_actual_family_input_fails(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_independent_fraction_continuous_and_resource_bounds():
    assert verify.checked_constants()["independent_Fraction_reconstruction_agrees"]
    assert len(ROWS) == 29
    assert len(audit.proof_checks()) == 17
    assert all(map(bool, audit.proof_checks().values()))
    assert audit.rejected_inputs() == 21


@pytest.mark.parametrize("xi", [0, Fraction(1, 6), Fraction(1, 4)])
@pytest.mark.parametrize("R", [1, 10, 10000])
def test_actual_family_keeps_negative_Wick_and_full_nonminimal_stress(R, xi):
    point = audit.family_point(R, xi)
    assert point["Wick_upper_in_hbar_gamma_squared_units"] < -sp.Rational(3, 5) * R * R
    assert (
        point["null_stress_upper_in_hbar_gamma_squared_units"]
        < -sp.Rational(3, 10) * R * R
    )
    assert point["one_sided_nonnegative_Wick_cap_satisfied_on_entire_sample_segment"]
    assert not point["uniform_energy_or_two_sided_Wick_cap_claimed"]


def test_full_improvement_is_not_minimal_tensor_at_conformal_coupling():
    minimal = audit.family_point(1, 0)
    conformal = audit.family_point(1, Fraction(1, 6))
    assert (
        minimal["null_stress_upper_in_hbar_gamma_squared_units"]
        != conformal["null_stress_upper_in_hbar_gamma_squared_units"]
    )
    assert (
        normalization.data()["checks"][
            "actual_null_contraction_matches_parent_covariant_stress"
        ]
        == 0
    )


def test_positive_finite_energy_grows_with_actual_transverse_radius():
    resources = normalization.resource_bounds()
    assert resources["energy_lower_coefficient_of_hbar_R_squared"] > 0
    assert resources["positive_ANEC_support_sum_gap"] > 0
    assert not resources[
        "uniform_total_energy_or_transverse_momentum_cutoff_for_the_family"
    ]
    assert not resources["global_two_sided_Wick_square_cap_for_the_family"]


@pytest.mark.parametrize("bad", [0.1, sp.Float("0.1"), sp.oo, sp.nan])
def test_inexact_nonfinite_serialization_fails(bad):
    with pytest.raises((ValueError, TypeError)):
        verify.serialize({"bad": bad})
