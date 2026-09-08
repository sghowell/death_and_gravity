"""Scientific acoustic mode, source-domain and momentum-bound checks."""
from fractions import Fraction

import pytest
import sympy as sp
from p8_vector_liouville import clock, proofs, reduction, variation, verify


def zero(values):
    for value in values.values():
        if isinstance(value, sp.MatrixBase):
            assert all(entry == 0 for entry in value)
        else:
            assert value == 0


def test_exact_two_sector_Liouville_and_canonical_identities():
    zero(reduction.checks())


def test_original_clock_potentials_replay_physical_modes():
    zero(clock.checks())


def test_prepared_acoustic_variation_and_history_identities():
    zero(variation.checks())


def test_all_exact_entries_and_continuous_gates():
    checks = proofs.residuals()
    zero(checks)
    assert len(checks) == 30
    assert sum(x.rows*x.cols if isinstance(x, sp.MatrixBase) else 1 for x in checks.values()) == 36
    assert len(proofs.checks()) == 12
    assert all(value is True for value in proofs.checks().values())


def test_independent_fraction_majorant_sum_includes_history_twice():
    varied = [Fraction(3343645513, 5308416), Fraction(40047147671, 28311552), Fraction(61181640625, 50331648)]
    history = [Fraction(18203125, 16384), Fraction(830078125, 262144), Fraction(91552734375, 33554432)]
    actual = variation.estimates()
    assert tuple(map(sp.Rational, varied)) == actual["varied_numerator_C2_bounds"]
    assert tuple(map(sp.Rational, history)) == actual["history_numerator_bounds"]
    assert sum(varied)+2*sum(history) == Fraction(11735921060393, 679477248)
    assert actual["correction_C2_coefficient"] == sp.Rational(sum(varied)+2*sum(history))
    assert sum(varied)+2*sum(history) < 18000


def test_exact_momentum_decay_and_nonzero_zero_momentum_limit():
    value0 = variation.bound(0, 1000)["prepared_remainder_C2_to_C0_upper"]
    assert variation.bound(1000, 1000)["prepared_remainder_C2_to_C0_upper"] == value0/2
    assert variation.bound(10000, 1000)["prepared_remainder_C2_to_C0_upper"] == value0/101
    assert value0 > 0


def test_both_third_source_derivatives_are_outside_norm():
    for value in (variation.n[3], variation.zeta[3]):
        with pytest.raises(ValueError, match="Third"):
            variation.coefficient_bound(value)


def test_fixed_acoustic_time_history_is_not_a_zero_term():
    data = variation.data()
    difference = data["delta_remainder_fixed_acoustic_time"]-data["delta_remainder_fixed_physical_time"]
    fixture = difference.subs({variation.xi: 1, clock.u: sp.Rational(1, 4), clock.mass2: 10**6, reduction.k: 1000})
    assert sp.factor(fixture) != 0
    assert difference.subs(variation.xi, 0) == 0


def test_source_chart_inverse_is_regular_and_small_at_rational_points():
    T = variation.source()["physical_to_acoustic_sources"]
    for u in (0, sp.Rational(1, 4), -sp.Rational(1, 2)):
        fixed = T.subs(clock.u, u)
        inverse = fixed.inv()
        assert fixed.det() < -1
        assert sum(abs(x) for x in inverse.row(0)) < 2
        assert sum(abs(x) for x in inverse.row(1)) == 1
        assert fixed*inverse == sp.eye(2)


def test_canonical_zero_momentum_is_not_invertible_but_potential_limit_is_regular():
    pump = reduction.data("L")["acoustic_pump_squared"]
    assert pump.subs(reduction.k, 0) == 0
    assert sp.factor(clock.data()["U"]) == clock.mass2*clock.a**2


def test_omitting_acoustic_rate_variation_has_nonzero_bounce_fixture():
    src, bg = variation.source(), clock.data()
    missing = bg["b"]*variation.acoustic(src["acoustic_log_rate"])+2*src["acoustic_log_rate"]*bg["pump_curvature"]
    value = sp.expand(missing.subs(clock.u, 0))
    assert value.coeff(variation.n[0]) == sp.Rational(616, 81)
    assert value.coeff(variation.zeta[0]) == -8


def test_physical_readout_prefactor_is_not_discardable():
    data = reduction.data("T")
    substitutions = {reduction.N: 2, reduction.a: 1, reduction.k: 1, reduction.m: 1, reduction.bm: 1}
    rate = data["acoustic_rate"].subs(substitutions)
    omega2 = data["physical_frequency_squared"].subs(substitutions)
    original = (1+omega2)/2
    intrinsic = (1/rate+omega2/rate)/2
    assert original == rate*intrinsic
    assert original-intrinsic == sp.Rational(9, 4)


def test_zero_polynomials_are_accepted_but_constants_and_mixed_sources_are_not():
    x = variation.n[0]
    assert variation.clean((x+1)**2-x*x-2*x-1) == 0
    assert variation.coefficient_bound(0) == 0
    for value in (1, variation.n[0]*variation.zeta[0], True):
        with pytest.raises((TypeError, ValueError)):
            variation.coefficient_bound(value)


def test_selected_mass_scalar_principal_is_strictly_positive():
    assert clock.bounds(0, 1000)["clock_scalar_principal_potential_lower"] == sp.Rational(15999725, 16)
    assert clock.bounds(0, 1000)["clock_scalar_principal_potential_lower"] > 0


def test_background_even_polynomial_endpoint_bounds():
    data = clock.data()
    assert (data["remainder_first_numerator"]/clock.mass2).subs(clock.u, sp.Rational(1, 2)) == sp.Rational(359375, 4096)
    assert data["pump_curvature"].subs(clock.u, sp.Rational(1, 2)) == sp.Rational(275, 16)
    assert clock.bounds(1000, 1000)["coarse_clock_remainder_upper"] == 44


def test_prepared_shift_and_potential_scope_are_explicit():
    assert variation.estimates()["initial_zero_time_shift_C0_upper_per_source_C0"] == 2
    assert "not the renormalized stress" in variation.estimates()["scope"]


def test_interface_controls_after_valid_cache_population():
    reduction.data("L")
    clock.bounds(0, 1000)
    variation.bound(1000, 1000)
    assert verify.controls()["rejected_inputs"] == 93
