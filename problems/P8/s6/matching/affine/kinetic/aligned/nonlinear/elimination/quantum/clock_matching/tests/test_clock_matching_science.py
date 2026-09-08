"""Proposed first-variation matching: no premature full quantum verdict."""
import sympy as sp
from p8_vector_clock_matching import bounds, continuation, counterterms, variation


def test_proposed_pole_continuation_and_finite_first_variations():
    assert all(value == 0 for value in counterterms.checks().values())


def test_continuous_local_coefficient_reconstruction():
    assert all(value == 0 for value in bounds.checks().values())
    assert sp.Rational(5, 2)+sp.Rational(22, 27) < 4


def test_candidate_one_loop_scale_bound_is_not_full_quantum_closure():
    d = bounds.physical_bounds(10**12, 1000)
    assert d["energy_candidate_matched_total_over_reference_density"] < sp.Rational(1, 10**14)
    assert d["pressure_candidate_matched_total_over_reference_density"] < sp.Rational(1, 10**14)
    assert d["unique_scheme_or_full_quantum_bounce_claim"] is False


def test_lapse_variation_precedes_clock_restriction_and_pressure_mass_term_vanishes():
    assert all(value == 0 for value in variation.checks().values())


def test_exact_complex_dimension_residual_and_analytic_mode_pair():
    for value in continuation.checks().values():
        assert all(entry == 0 for entry in value) if isinstance(value, sp.MatrixBase) else value == 0


def test_uniform_complex_dimension_bounds_make_the_subtracted_integral_convergent():
    assert all(value is True for value in continuation.proof_checks().values())


def test_actual_local_coefficients_at_mu_mass():
    u = continuation.u
    actual = bounds.actual_coefficients()
    expected = {
        1: {"energy": -32*(1215*u**8+3645*u**6+3645*u**4+1136*u**2+7)/(243*(1+u**2)**5),
            "pressure": 80*(5*u**2+1)/(3*(1+u**2)**2)},
        2: {"energy": -32*(1701*u**10+3969*u**8+1782*u**6-1009*u**4-869*u**2+22)/(81*(1+u**2)**7),
            "pressure": 64*(7*u**4-1)/(1+u**2)**4}}
    for n, values in expected.items():
        for name, value in values.items():
            assert sp.factor(actual[n][name]-value) == 0


def test_complex_dimension_box_counts_modulus_not_real_sign():
    u, z, D = continuation.u, continuation.z, continuation.local.dimension
    expression = (3*u*z*(D-3)**2-2*(D-3)+1)/(1+u**2)**2
    result = continuation.box_bound(expression)
    assert result["reconstruction"] == 0
    assert result["absolute_upper"] == sp.Rational(51, 32)
    at_complex_boundary = expression.subs({u: sp.Rational(1, 2), z: 1, D: 3+sp.I/4})
    assert sp.simplify(sp.Abs(at_complex_boundary)**2) < result["absolute_upper"]**2


def test_finite_bound_decreases_with_planck_scale_at_fixed_mass_time():
    first = bounds.physical_bounds(10**12, 1000)
    second = bounds.physical_bounds(2*10**12, 1000)
    for name in ("energy", "pressure"):
        key = name+"_candidate_matched_total_over_reference_density"
        assert 4*second[key] == first[key]
