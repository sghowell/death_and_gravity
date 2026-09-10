"""Independent actual Wick, generated-interaction, source and local-bound checks."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_offshell_vacuum import jets
from p8_vacuum_quantum_map import audit, calibration, interactions, ward, wick


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_derivative_map_and_one_loop_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_explicit_perturbative_or_continuous_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_reject_unsupported_map_input(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_control_and_original_scope(name, value):
    assert bool(value), name


@pytest.mark.parametrize("invariant", (0, 1, 2, Fraction(1, 2), sp.Rational(3, 2)))
def test_actual_real_slice_inside_complex_disc_bound(invariant):
    d = calibration.point(invariant)
    assert d["mass_squared_invariant"] == sp.Rational(invariant)
    # Use the rigorous pi>3 rational upper, not a decimal evaluation.
    finite = d["finite_linear_composite_mixing"]
    rational_numerator = sp.factor(16 * sp.pi**2 * finite)
    assert abs(rational_numerator) / 144 < d["whole_complex_disc_mixing_upper"]


def test_input_validation_not_bypassed_by_equal_cached_float():
    calibration.point(1)
    for value in (True, 1.0, sp.Float(1)):
        with pytest.raises((TypeError, ValueError)):
            calibration.point(value)


@pytest.mark.parametrize("axis", range(4))
def test_literal_two_derivative_covariance_signs(axis):
    zero = (0, 0, 0, 0)
    one = tuple(int(i == axis) for i in range(4))
    two = tuple(2 * int(i == axis) for i in range(4))
    assert wick.covariance(one, one) == sp.Rational(jets.SIGNS[axis], 4)
    assert wick.covariance(zero, two) == -sp.Rational(jets.SIGNS[axis], 4)
    assert wick.covariance(zero, one) == 0


def test_rank_four_metric_pairing_multiplicity():
    zero = (0, 0, 0, 0)
    assert wick.covariance(zero, (4, 0, 0, 0)) == sp.Rational(1, 8)
    assert wick.covariance(zero, (2, 2, 0, 0)) == -sp.Rational(1, 24)
    assert wick.covariance(zero, (0, 4, 0, 0)) == sp.Rational(1, 8)


def test_literally_contract_actual_R_not_surrogate():
    d = wick.data()
    assert d["actual_literal_jet_contraction"] == wick.contract(
        jets.data()["cubic_field_redefinition"]
    )
    assert d["literal_jet_count"] == 210


def test_evanescent_finite_term_is_nonzero():
    d = wick.data()
    B = sp.Symbol("free_Box_eigenvalue", real=True)
    gamma = sp.Symbol("quartic_gamma", real=True)
    naive = -d["four_dimensional_coefficient"].subs(B, -1) / (16 * sp.pi**2)
    difference = sp.expand(d["on_shell_one_loop_overlap_correction"] - naive)
    assert difference == 3 * gamma / (64 * sp.pi**2)
    assert difference != 0


def test_one_loop_sextic_topology_not_discarded():
    assert set(interactions.topologies(4, 1)) == {(2, 0, 0, 0, 0), (0, 1, 0, 0, 0)}
    assert all(row[2:] == (0, 0, 0) for row in interactions.topologies(4, 1))


@pytest.mark.parametrize("external", (2, 4, 6))
@pytest.mark.parametrize("loops", (0, 1, 2))
def test_topology_half_edge_euler_identity(external, loops):
    for row in interactions.topologies(external, loops):
        assert (
            sum((n - 2) * v for n, v in zip((4, 6, 8, 10, 12), row))
            == 2 * loops + external - 2
        )


def test_all_higher_generated_interactions_nonzero():
    rows = interactions.data()["all_generated_field_degrees"]
    assert set(rows) == {2, 4, 6, 8, 10, 12}
    assert all(value != 0 for value in rows.values())


def test_actual_positive_nontrivial_residue_and_disc_margin():
    d = calibration.data()
    assert 0 < d["actual_on_shell_overlap_rational_upper"] < sp.Rational(1, 10**404)
    assert 0 < d["mapped_residue_rational_lower"] < 1
    assert (
        0 < d["mapped_unit_disc_inverse_factor_error_upper"] < sp.Rational(3, 10**404)
    )


def test_ordinary_amputation_not_physical_amplitude():
    d = ward.data()
    assert (
        sp.expand(d["ordinary_amputated_four_point_vertex"] - d["proper_LSZ_amplitude"])
        != 0
    )


def test_gaussian_source_Jacobian_and_sextic_controls_are_not_zero():
    d = ward.data()
    for key in (
        "omitting_Jacobian_changes_partition_function",
        "omitting_transformed_source_changes_two_point_function",
        "omitting_generated_sextic_changes_partition_function",
    ):
        assert d[key] != 0


def test_off_shell_effective_action_counterexample():
    assert ward.data()["off_shell_effective_action_difference_quartic"] != 0


def test_exact_counts_and_unclosed_original_scope():
    assert len(audit.residuals()) == 55
    assert len(audit.gates()) == 32
    assert len(audit.controls()) == 11
    assert audit.rejected_inputs() == 38
