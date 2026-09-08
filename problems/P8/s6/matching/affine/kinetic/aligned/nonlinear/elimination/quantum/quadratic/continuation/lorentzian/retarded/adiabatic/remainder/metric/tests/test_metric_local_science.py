"""Scientific full-metric local identities and explicit failure controls."""
import sympy as sp
from p8_vector_metric_local import (
    bounds,
    canonical,
    consistency,
    controls,
    counterterms,
    jets,
    proofs,
    radial,
    variation,
    verify,
)


def zero(values):
    for value in values.values():
        if isinstance(value, sp.MatrixBase):
            assert all(entry == 0 for entry in value)
        else:
            assert value == 0


def test_independent_physical_Hamiltonian_vertices():
    zero(canonical.checks())
    assert len(canonical.checks()) == 22


def test_full_metric_varied_Riccati_and_source_linearity():
    zero(variation.checks())
    assert len(variation.checks()) == 20


def test_actual_second_mass_vertex_is_not_optional():
    assert controls.second_mass_vertex()["bounce_unit_source_fixture"] == -sp.Rational(826, 243)


def test_proper_time_geometry_before_metric_variation():
    zero(counterterms.controls())


def test_six_full_curved_poles_match_fixed_counterterms():
    zero(radial.pole_checks())
    assert len(radial.pole_checks()) == 6


def test_full_potential_lapse_and_logscale_hessian():
    zero(consistency.potential_checks())


def test_independent_ordinary_finite_heat_action_limit():
    zero(consistency.ordinary_checks())


def test_full_weighted_self_adjoint_matrix():
    zero(consistency.self_adjoint_checks())
    assert len(consistency.self_adjoint_checks()) == 51


def test_constant_spatial_rescaling_of_density_currents():
    zero(consistency.scale_checks())


def test_physical_stress_normalization_is_varied():
    zero(bounds.checks())


def test_fourth_derivative_symbol_is_local_not_resummed_spectrum():
    zero(controls.checks())
    data = controls.top_derivative()
    assert data["determinant"].subs(jets.u, 0) == -sp.Rational(50992, 59049)


def test_omitted_evanescent_counterterms_fail_mixed_adjoint_test():
    data = controls.omitted_counterterm()
    assert data[1]["bounce_unit_mass_fixtures"][0] == -sp.Rational(1600, 81)
    assert data[2]["bounce_unit_mass_fixtures"][0] == -sp.Rational(64, 3)
    assert data[2]["bounce_unit_mass_fixtures"][2] == -sp.Rational(512, 27)


def test_continuous_rational_majorants_not_point_samples():
    assert all(value >= 0 and value.is_Rational
               for table in bounds.envelopes().values() for sources in table.values()
               for row in sources.values() for value in row.values())


def test_full_local_bound_and_exact_planck_scaling():
    value = bounds.operator_bound(10**24, 1000)["joint_physical_C4_to_C0_bound"]
    assert 0 < value < sp.Rational(1, 10**38)
    assert bounds.operator_bound(2*10**24, 1000)["joint_physical_C4_to_C0_bound"] == value/4


def test_strict_input_controls_after_valid_caches_are_populated():
    assert verify.controls()["rejected_inputs"] == 126


def test_all_exact_identities_and_audit_gates():
    residuals = proofs.residuals()
    zero(residuals)
    assert len(residuals) == 132
    assert sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in residuals.values()) == 135
    gates = proofs.checks()
    assert len(gates) == 22
    assert all(value is True for value in gates.values())
