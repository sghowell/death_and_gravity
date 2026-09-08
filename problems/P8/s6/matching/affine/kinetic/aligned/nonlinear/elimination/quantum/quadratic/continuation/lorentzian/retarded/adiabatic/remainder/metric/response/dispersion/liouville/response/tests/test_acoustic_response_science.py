"""Scientific moving-clock, selected-state and physical-readout checks."""
import pytest
import sympy as sp
from p8_vector_acoustic_response import (
    dimensional,
    proofs,
    readout,
    source,
    transport,
    verify,
)
from p8_vector_liouville import reduction, variation


def zero(values):
    for value in values.values():
        if isinstance(value, sp.MatrixBase):
            assert all(entry == 0 for entry in value)
        else:
            assert value == 0


def test_covariance_and_tangent_transport_are_exact():
    zero(transport.checks())


def test_retarded_green_and_selected_variance_kernel():
    zero(transport.green_checks())


def test_dimensional_pumps_and_source_jets_are_retained():
    zero(dimensional.checks())


def test_both_sources_replay_original_physical_mode_variation():
    zero(source.checks())


def test_independent_action_vertices_and_generic_readout_conversion():
    zero(readout.checks())


def test_all_four_physical_readout_and_contact_rows():
    zero(readout.actual_checks())


def test_complete_exact_and_scope_counts():
    checks = proofs.residuals()
    zero(checks)
    assert len(checks) == 46
    assert sum(v.rows*v.cols if isinstance(v, sp.MatrixBase) else 1 for v in checks.values()) == 82
    assert len(proofs.checks()) == 13
    assert all(value is True for value in proofs.checks().values())


def test_two_prepared_clock_shifts_are_not_identified():
    assert source.data("T")["prepared_shift"] != source.data("L")["prepared_shift"]
    difference = source.clean(source.data("L")["log_rate_variation"]-source.data("T")["log_rate_variation"])
    assert difference == -4*variation.n[0]/(81*(1+source.u**2)**3)


def test_covariance_sources_have_only_first_source_jets_before_stress_subtraction():
    for sector in ("T", "L"):
        matrix = source.data(sector)["delta_matrix_fixed_sigma"]
        assert not matrix.has(variation.n[2], variation.zeta[2], variation.n[3], variation.zeta[3])


def test_zero_varied_initial_source_does_not_set_unvaried_covariance_to_zero():
    for sector in ("T", "L"):
        data = source.data(sector)
        initial = {field: 0 for field in source.fields}
        assert data["delta_matrix_fixed_sigma"].subs(initial) == sp.zeros(3)
        assert data["matrix"] != sp.zeros(3)
        assert data["D"] > 0


def test_physical_readout_rows_have_zero_cross_column_before_history():
    for sector in ("T", "L"):
        for component in ("energy", "pressure"):
            item = readout.actual_readouts()[sector][component]
            assert item["response_row"][0, 1] == 0
            assert item["fixed_output_contact_row"][0, 1] == 0
            assert item["prepared_history_row"] != sp.zeros(1, 3)


def test_independent_squeezing_changes_variance_response_fixture():
    f, g = (1-sp.I)/2, 1/sp.sqrt(2)
    a, b = sp.Rational(5, 4), sp.Rational(3, 4)
    def kernel(left, right):
        return 2*sp.im((left*sp.conjugate(right))**2)
    squeezed_f, squeezed_g = a*f+b*sp.conjugate(f), a*g+b*sp.conjugate(g)
    assert sp.simplify(kernel(squeezed_f, squeezed_g)-kernel(f, g)) == -sp.Rational(3, 2)


def test_same_squeezing_preserves_the_causal_green_function():
    f, g = (1-sp.I)/2, 1/sp.sqrt(2)
    a, b = sp.Rational(5, 4), sp.Rational(3, 4)
    def green(left, right):
        return sp.I*(left*sp.conjugate(right)-sp.conjugate(left)*right)
    assert a*a-b*b == 1
    assert sp.simplify(green(a*f+b*sp.conjugate(f), a*g+b*sp.conjugate(g))-green(f, g)) == 0


def test_transverse_dimensional_jet_survives_physical_pump_simplification():
    value = dimensional.pumps()["T"]
    assert value.subs(dimensional.D, 3) == 1
    assert sp.diff(value, dimensional.D).subs({dimensional.D: 3, reduction.a: 2}) == sp.log(2)
    assert sp.log(2) != 0


def test_generic_rate_contact_is_not_optional():
    value = readout.algebra()["contact_only"]
    substitutions = {symbol: 0 for symbol in value.free_symbols}
    names = {"acoustic_rate": 1, "varied_log_acoustic_rate": 1, "readout_A": 1,
             "readout_B": 1, "U": 1, "k": 1, "S": 1, "P": 1}
    substitutions.update({symbol: names[str(symbol)] for symbol in value.free_symbols if str(symbol) in names})
    assert sp.simplify(value.subs(substitutions, simultaneous=True)) == sp.Rational(3, 2)


def test_source_cleaning_rejects_foreign_or_nonlinear_terms():
    assert source.clean(0) == 0
    assert source.clean(source.xi_L) == source.xi_L
    for value in (True, sp.Float(1), variation.n[0]*source.xi_L, sp.Symbol("foreign")*variation.n[0]):
        with pytest.raises((TypeError, ValueError)):
            source.clean(value)


def test_negative_half_rate_is_not_a_failed_positive_clock():
    data = source.data("L")
    fixture = {source.u: -sp.Rational(1, 4), source.mass2: 10**6, reduction.k: 1000}
    assert data["bk"].subs(fixture) < 0
    assert data["rate"].subs(fixture) > 0
    assert data["D"].subs(fixture) > 0


def test_strict_interface_controls_after_cache_population():
    source.data("L")
    source.derivative(variation.n[0], "T")
    assert verify.controls()["rejected_inputs"] == 40
