"""Exact phase pairs, actual endpoint map, and separated scalar readout."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_own_forced import loading, phase
from p8_own_scattering import connection


def test_exact_primary_phase_identities():
    assert phase.identities()
    assert all(value == 0 for value in phase.identities().values())


def test_all_primary_phase_margins():
    assert all(phase.checks().values())


def test_independent_recomputed_half_width_budget():
    f = Fraction
    a, delta = f(1, 200), f(1, 10**6)
    central = 44 * (a**2 + 2 * a / 2000)
    tails = 3 * delta / (32 * a**2)
    budget = f(4, 5) * (central + tails)
    assert (central, tails, budget) == (f(33, 25000), f(3, 800), f(507, 125000))
    assert budget < f(1, 200)
    assert 4 * (1 / (1 - f(1, 400)) - 1) == f(4, 399)
    assert phase.calibration()["perturbation_budget_upper"] == budget


def test_direct_complex_frame_gives_the_real_frequency_dictionary():
    ar, ai, beta = sp.symbols("ar ai beta", real=True)
    aa = ar + sp.I * ai
    v0 = sp.Matrix([[1, 1], [sp.I, -sp.I]]) / sp.sqrt(2)
    transfer = sp.Matrix([[sp.conjugate(aa), sp.I * beta], [-sp.I * beta, aa]])
    real = v0 * transfer * v0.H
    assert sp.simplify(real - sp.Matrix([[ar, -ai - beta], [ai - beta, ar]])) == sp.zeros(2)


def test_phase_separation_retains_constant_B_and_flips_A():
    data = phase.generic_reference()
    theta = data["theta"]
    plus = data["reference"].subs(theta, 0)
    minus = data["reference"].subs(theta, sp.pi / 2)
    assert sp.simplify((plus + minus) / 2 - data["constant"]) == sp.zeros(2)
    assert sp.simplify((plus - minus) / 2 - data["oscillatory"].subs(theta, 0)) == sp.zeros(2)


def test_alignment_can_target_original_Q_not_just_a_chosen_vector_norm():
    # A=(3+4i)/5 and w=(2,1); this uses a non-axis-aligned loaded state.
    aa = (3 + 4 * sp.I) / 5
    ww = 2 + sp.I
    z = sp.expand(aa * ww)
    norm = sp.sqrt(sp.simplify(z * sp.conjugate(z)))
    phasor = sp.conjugate(z) / norm
    assert sp.simplify(phasor * sp.conjugate(phasor)) == 1
    assert sp.simplify(aa * phasor * ww - norm) == 0
    assert sp.simplify(-aa * phasor * ww + norm) == 0


def test_mixing_B_is_not_the_reason_for_the_phase_nonconvergence():
    data = phase.generic_reference()
    zero_mix = data["reference"].subs({data["ar"]: 1, data["ai"]: 0, data["beta"]: 0})
    assert zero_mix.subs(data["theta"], 0) == sp.eye(2)
    assert zero_mix.subs(data["theta"], sp.pi / 2) == -sp.eye(2)


@pytest.mark.parametrize("n", [2, 3, 7])
def test_exact_phase_sequence_formulas(n):
    values = phase.phase_sequences(n, sp.pi / 7)
    assert sp.simplify(phase.phase_value(values["delta_plus"]) - values["theta_plus"]) == 0
    assert sp.simplify(phase.phase_value(values["delta_minus"]) - values["theta_minus"]) == 0
    assert values["theta_minus"] - values["theta_plus"] == sp.pi / 2


def test_sequence_domain_has_a_rational_analytic_proof():
    f = Fraction
    sinh_lower = f(9, 2) + f(9, 2)**3 / 6
    assert sinh_lower == f(315, 16)
    upper = 8 * f(1, 200)**2 / sinh_lower**2
    assert upper < f(1, 10**6)
    assert phase.calibration()["sequence_delta_upper"] == upper


@pytest.mark.parametrize("bad", [True, 2.0, sp.Float(2), 1, 0, -1])
def test_phase_index_guard(bad):
    with pytest.raises((TypeError, ValueError)):
        phase.phase_sequences(bad)


@pytest.mark.parametrize("bad", [True, 0.0, sp.Float(".1"), sp.oo, sp.nan, -1])
def test_phase_offset_guard(bad):
    with pytest.raises((TypeError, ValueError)):
        phase.phase_sequences(2, bad)


@pytest.mark.parametrize("side", [-1, 1])
def test_actual_endpoint_limit_matches_frozen_full_map(side):
    delta, k = sp.symbols("delta k", positive=True)
    ku = sp.Symbol("k_u", real=True)
    full = connection.endpoint_map(side * sp.Rational(loading.A), delta, k, ku)
    limit = full.subs(delta, 0)
    ours = phase.endpoint_limit(k, ku, side=side)
    assert sp.simplify(limit - ours) == sp.zeros(2)
    assert sp.simplify(ours.det() - k / phase.RHO) == 0


def test_Q_first_row_inversion_does_not_omit_velocity_from_full_state():
    k = sp.Symbol("k", positive=True)
    ku = sp.Symbol("k_u", real=True)
    matrix = phase.endpoint_limit(k, ku)
    assert sp.simplify(matrix.inv()[0, 0] - sp.sqrt(sp.Rational(loading.A) / k)) == 0
    assert matrix.inv()[0, 1] == 0
    assert sp.diff(matrix[1, 0], ku) != 0


def test_independent_scalar_gap_arithmetic():
    f = Fraction
    paired_error = 2 * f(4, 399)
    separation = 2 - paired_error
    loading_floor = f(9, 1280)
    reciprocal_root_kinetic_floor = f(2, 5)
    gap = separation * loading_floor * reciprocal_root_kinetic_floor
    assert separation == f(790, 399)
    assert gap == f(237, 42560) > f(1, 200)
    assert phase.calibration()["scalar_Q_gap_lower_over_eta"] == gap


def test_zero_loaded_data_control_has_no_oscillatory_response():
    data = phase.generic_reference()
    assert data["reference"] * sp.zeros(2, 1) == sp.zeros(2, 1)
    assert data["alignment_squared_norm"].subs({data["w1"]: 0, data["w2"]: 0}) == 0


def test_no_exact_subsequence_or_corner_convergence_is_assumed():
    data = phase.calibration()
    assert data["subsequential_reference_limits_asserted"] is False
    assert data["central_remainder_convergence_assumed"] is False
    assert data["matter_source_or_EFT_verdict"] is False
