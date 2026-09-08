"""Scientific prepared full-metric finite response and fixed-profile controls."""
import sympy as sp
from p8_vector_metric_response import (
    envelopes,
    evolution,
    proofs,
    source,
    tadpole,
    tail,
    tangent,
    verify,
)


def zero(values):
    for value in values.values():
        if isinstance(value, sp.MatrixBase):
            assert all(entry == 0 for entry in value)
        else:
            assert value == 0


def test_full_dimensional_low_order_variation_is_reproduced():
    zero(tangent.low_checks())


def test_varied_Riccati_residual_starts_at_inverse_frequency_eighth():
    for sector in ("T", "L"):
        zero(tangent.reference(sector)["varied_residual_low_coefficients"])


def test_four_physical_readout_tails_cancel_through_order_four():
    for sector in ("T", "L"):
        for component in ("energy", "pressure"):
            zero(tail.reference_tail(sector, component)["low_tail_numerator_coefficients"])


def test_reference_variation_is_linear_in_both_prepared_sources():
    for sector in ("T", "L"):
        for j in range(1, 5):
            value = tangent.coefficient(sector, j)
            assert source.clean(sum(field*sp.diff(value, field) for field in source.n+source.v)-value) == 0


def test_isotropic_zero_momentum_reference_limit():
    for j in range(1, 5):
        assert source.clean((tangent.varied_P("T", j)-tangent.varied_P("L", j)).subs(source.z, 0)) == 0


def test_tenth_derivatives_are_retained_for_both_metric_sources():
    for sector in ("T", "L"):
        value = tangent.reference(sector)["tenth_source_derivative_fixtures"]
        assert value == {"N": -sp.Rational(7, 20736), "Z": -sp.Rational(1, 512)}


def test_sixth_and_eighth_reference_source_jet_orders():
    for sector in ("T", "L"):
        for j in (3, 4):
            assert source.linear_bound(tangent.coefficient(sector, j))["highest_source_derivative"] == 2*j


def test_continuous_source_and_readout_majorants():
    for sector in ("T", "L"):
        assert all(value == 0 for value in envelopes.reference(sector)["box_reconstructions"])
        for component in ("energy", "pressure"):
            assert tail.reference_tail(sector, component)["all_majorants_nonnegative"]


def test_closed_physical_weights_support_the_reference_bound():
    zero(evolution.weight_checks())


def test_exact_mixing_Wronskian_phase_and_radial_algebra():
    zero(evolution.algebra_checks())
    zero(tail.algebra_checks())


def test_existing_tadpole_has_fixed_profile_physical_vertices():
    data = tadpole.physical_vertices()
    zero(data["checks"])
    assert data["physical_stress_jacobian"] == sp.ImmutableMatrix(
        [[data["rho"]+data["pressure"], 0], [data["rho"]+data["pressure"], 0]])


def test_finite_full_vector_metric_operator_bound():
    data = evolution.bound(10**24, 1000)
    assert data["joint_nonlocal_C10_to_C0_upper_bound"] < sp.Rational(1, 10**42)
    assert data["joint_complete_metric_C10_to_C0_upper_bound"] < sp.Rational(1, 10**38)


def test_fixed_tadpole_is_included_in_background_cancelled_bound():
    data = tadpole.bound(10**24, 1000)
    assert data["joint_background_cancelled_metric_C10_to_C0_upper_bound"] < sp.Rational(1, 10**37)
    assert data["joint_background_cancelled_metric_C10_to_C0_upper_bound"] > evolution.bound(
        10**24, 1000)["joint_complete_metric_C10_to_C0_upper_bound"]


def test_exact_planck_scaling_for_vector_and_fixed_profile_response():
    key = "joint_complete_metric_C10_to_C0_upper_bound"
    assert evolution.bound(2*10**24, 1000)[key] == evolution.bound(10**24, 1000)[key]/4
    key = "joint_background_cancelled_metric_C10_to_C0_upper_bound"
    assert tadpole.bound(2*10**24, 1000)[key] == tadpole.bound(10**24, 1000)[key]/4


def test_strict_source_and_interface_controls_after_cache_population():
    assert verify.controls()["rejected_inputs"] == 204


def test_all_exact_identities_and_continuous_audit_gates():
    identities = proofs.residuals()
    zero(identities)
    assert len(identities) == 59
    assert sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values()) == 69
    checks = proofs.checks()
    assert len(checks) == 38
    assert all(value is True for value in checks.values())
