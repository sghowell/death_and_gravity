"""Independent physical-loading, phase, and bounded-error audit for S6.36."""

import hashlib
from fractions import Fraction
from pathlib import Path

import pytest
import sympy as sp
from p8_own_forced import audit, loading, phase
from p8_own_scattering import connection

PARENT_SHA = "c77d6d3e41d81c5d2bda33aeba201fe438620d872f3e9f2d49fde5469194641e"


def _zero_matrix(matrix):
    return all(sp.simplify(sp.expand_complex(value)) == 0 for value in matrix)


def test_frozen_connection_report_is_the_named_parent():
    parent = Path(__file__).resolve().parents[2]
    report = parent / "certificates/own-fixed-window-transfer.json"
    assert hashlib.sha256(report.read_bytes()).hexdigest() == PARENT_SHA


def test_independent_fraction_source_and_positive_Green_bounds():
    f = Fraction
    length, a, delta = f(1, 100), f(1, 200), f(1, 10**6)
    coframe_low, coframe_high = f(19, 10), f(21, 10)
    k_low, k_high = coframe_low**3 / coframe_high, coframe_high**3 / coframe_low
    pole_low = 16 / (delta + 8 * length**2)
    ceiling = 16 / (8 * a**2) + 44
    derivative = 1 - ceiling * (length / 4)**2 / 2
    ratio = 1 - ceiling * (length / 4)**2 / 6
    assert k_low > 1 and k_high < 5 < f(25, 4)
    assert pole_low - 44 > 0
    assert pole_low - 44 - 22 > 1 / length**2
    assert ceiling == 80044
    assert derivative == f(59989, 80000) > f(2, 3)
    assert ratio == f(219989, 240000) > f(9, 10)
    primary, separate = loading.calibration(), audit.constants()
    assert primary["kinetic_lower"] == separate["strict_kinetic_lower"] == k_low
    assert primary["kinetic_upper"] == separate["strict_kinetic_upper"] == k_high
    assert primary["green_derivative_lower"] == separate["Green_derivative_floor"] == derivative
    assert primary["green_ratio_lower"] == separate["Green_ratio_floor"] == ratio
    assert primary["mass_lower"] == pole_low - 66
    assert primary["canonical_source_lower"] == 1 / length**2


def test_canonical_source_and_pump_derive_from_the_own_metric_operator():
    u = sp.symbols("u", real=True)
    k = sp.Function("k", positive=True)(u)
    mass = sp.Function("m")(u)
    q, yy = sp.Function("q")(u), sp.Function("Y")(u)
    qq = yy / sp.sqrt(k)
    literal = (sp.diff(k * sp.diff(qq, u), u) + k * mass * (qq - q)) / sp.sqrt(k)
    potential = mass - sp.diff(sp.sqrt(k), u, 2) / sp.sqrt(k)
    canonical = sp.diff(yy, u, 2) + potential * yy - sp.sqrt(k) * mass * q
    assert sp.simplify(literal - canonical) == 0
    # Neither the actual source sqrt(k)*mass*q nor the pump may be deleted.
    omitted_source = sp.diff(yy, u, 2) + potential * yy - q
    assert sp.simplify(literal - omitted_source - (1 - sp.sqrt(k) * mass) * q) == 0
    omitted_pump = sp.diff(yy, u, 2) + mass * yy - sp.sqrt(k) * mass * q
    assert sp.simplify(literal - omitted_pump + sp.diff(sp.sqrt(k), u, 2) * yy / sp.sqrt(k)) == 0


def test_Green_integrals_and_the_support_dependent_lag_are_independent():
    h, t, vmax = sp.symbols("h t Vmax", nonnegative=True)
    position = h - sp.integrate((h - t) * vmax * t, (t, 0, h))
    velocity = 1 - sp.integrate(vmax * t, (t, 0, h))
    assert sp.expand(position - h + vmax * h**3 / 6) == 0
    assert sp.expand(velocity - 1 + vmax * h**2 / 2) == 0
    assert sp.diff(position, h) == velocity
    length = Fraction(1, 100)
    ceiling = 8 / length**2 + 44
    assert loading.volterra_bounds(ceiling, length / 4)["positive_derivative_bootstrap"]
    # The source starts later than the zero-data slice. Extending this crude
    # bootstrap to the full L/2 initial-to-loading interval would fail.
    assert not loading.volterra_bounds(ceiling, length / 2)["positive_derivative_bootstrap"]


def test_the_explicit_fixed_pulse_has_a_proved_area_not_a_sampled_integral():
    f = Fraction
    length = f(1, 100)
    left, right, ramp = -95 * length / 128, -81 * length / 128, length / 64
    assert -3 * length / 4 < left < right < -5 * length / 8
    assert right - left - 2 * ramp == 5 * length / 64 > length / 16
    assert length / 8 < -length / 2 - right
    assert -length / 2 - left < length / 4
    actual = loading.pulse()
    assert (actual["support_left"], actual["support_right"], actual["ramp"]) == (left, right, ramp)
    u, eta = actual["u"], actual["eta"]
    assert sp.simplify(actual["q"].subs(u, sp.Rational((left + right) / 2)) - eta) == 0
    assert actual["q"].subs(u, sp.Rational(-length)) == 0
    assert actual["q"].subs(u, sp.Rational(-length / 2)) == 0
    assert actual["q"].free_symbols == {u, eta}


def test_stronger_uniform_loading_preserves_strictness_after_the_limit():
    f = Fraction
    length = f(1, 100)
    sharp_y = f(219989, 240000) * (length / 8) * (length / 16) / length**2
    sharp_velocity = f(59989, 80000) * (length / 16) / length**2
    assert sharp_y - f(9, 1280) == f(3989, 30720000) > 0
    assert sharp_velocity > 1 / (24 * length)
    assert loading.calibration()["Y_load_sharper_lower_over_eta"] == sharp_y
    assert audit.constants()["sharp_loaded_Y_floor_per_eta"] == sharp_y
    # These fixed gaps, not lim(strict)>strict, justify the loaded-limit claim.
    assert sharp_y > loading.calibration()["Y_load_lower_over_eta"]


def test_loading_coefficients_have_a_fixed_punctured_denominator():
    v, delta = sp.symbols("v delta", nonnegative=True)
    denominator = 2 + delta - 2 / (1 + v)**4
    assert sp.diff(denominator, delta) == 1
    assert sp.diff(denominator, v) == 8 / (1 + v)**5
    a = Fraction(1, 200)
    floor = 2 - 2 / (1 + a**2)**4
    assert floor > 0
    assert loading.source_bounds()["punctured_D_lower"] == floor
    # The continuous limit used here is not a literal delta=0 action at u=0.
    assert denominator.subs({v: 0, delta: 0}) == 0


def test_two_endpoint_errors_and_half_width_budget_are_recomputed():
    f = Fraction
    a, delta, radius = f(1, 200), f(1, 10**6), f(1, 2000)
    assert delta / 8 < radius**2
    central = 44 * (a**2 + 2 * a * radius)
    tails = 3 * delta / (32 * a**2)
    budget = (central + tails) * f(4, 5)
    assert central == f(33, 25000)
    assert tails == f(3, 800)
    assert budget == f(507, 125000) < f(1, 200)
    error = 4 * (1 / (1 - f(1, 400)) - 1)
    assert error == f(4, 399)
    assert 2 - 2 * error == f(790, 399)
    primary, separate = phase.calibration(), audit.constants()
    assert primary["perturbation_budget_upper"] == separate["budget"] == budget
    assert primary["transfer_error_upper"] == separate["transfer_error_upper"] == error
    assert primary["normalized_pair_separation_factor"] == separate["balanced_pair_separation_factor"] == 2 - 2 * error


def test_frozen_time_transfer_and_both_endpoint_frames_give_the_real_map():
    ar, ai, beta, theta = sp.symbols("ar ai beta theta", real=True)
    aa = ar + sp.I * ai
    frozen = connection.coefficients()
    tt = connection.time_transfer().xreplace({
        frozen["A"]: aa, sp.conjugate(frozen["A"]): sp.conjugate(aa),
    }).subs(sp.sinh(sp.pi * connection.RHO), 1 / beta)
    vplus = connection.plane_wave_frame(theta / connection.RHO)
    vminus = connection.plane_wave_frame(-theta / connection.RHO)
    actual = vplus * tt * vminus.H
    cc, ss = sp.cos(2 * theta), sp.sin(2 * theta)
    rr, ii = ar * cc + ai * ss, ar * ss - ai * cc
    expected = sp.Matrix([[rr, ii - beta], [-ii - beta, rr]])
    assert _zero_matrix(actual - expected)
    generic = phase.generic_reference()
    assert _zero_matrix(generic["reference"].subs({
        generic["ar"]: ar, generic["ai"]: ai,
        generic["beta"]: beta, generic["theta"]: theta,
    }) - expected)
    independent = audit.phase_algebra()
    symbols = {str(symbol): symbol for symbol in independent["real_map"].free_symbols}
    assert _zero_matrix(independent["real_map"].subs({
        symbols["ar"]: ar, symbols["ai"]: ai, symbols["beta"]: beta,
        symbols["cos_phase"]: cc, symbols["sin_phase"]: ss,
    }) - expected)


def test_reference_cauchy_API_uses_time_not_radial_transfer():
    theta = sp.symbols("theta", real=True)
    direct = (connection.plane_wave_frame(theta / connection.RHO)
              * connection.time_transfer()
              * connection.plane_wave_frame(-theta / connection.RHO).H)
    assert all(sp.simplify(value) == 0 for value in phase.reference_cauchy(theta) - direct)
    ar, ai, beta = sp.symbols("ar ai beta", real=True)
    time = sp.Matrix([[ar - sp.I * ai, sp.I * beta], [-sp.I * beta, ar + sp.I * ai]])
    radial = sp.Matrix([[sp.I * beta, ar - sp.I * ai], [ar + sp.I * ai, -sp.I * beta]])
    assert sp.expand(time.det() + radial.det()) == 0
    assert sp.expand(time.det() - ar**2 - ai**2 + beta**2) == 0


def test_exact_Q_phase_offset_and_constant_mixing_cancellation():
    ar, ai, w1, w2, beta = sp.symbols("ar ai w1 w2 beta", real=True)
    cc, ss = sp.symbols("c s", real=True)
    c0, s0 = ar * w1 - ai * w2, ai * w1 + ar * w2
    amplitude_squared = (ar**2 + ai**2) * (w1**2 + w2**2)
    assert sp.expand(c0**2 + s0**2 - amplitude_squared) == 0
    first = cc * c0 + ss * s0 - beta * w2
    amplitude = sp.symbols("amplitude", positive=True)
    aligned = first.subs({cc: c0 / amplitude, ss: s0 / amplitude})
    assert sp.expand(amplitude * (aligned + beta * w2) - amplitude_squared) == 0
    flipped = first.subs({cc: -cc, ss: -ss}, simultaneous=True)
    assert sp.expand(first - flipped - 2 * (cc * c0 + ss * s0)) == 0
    # The constant B term cancels. Nonzero B is not necessary for this phase
    # phenomenon; it must not be identified with the oscillatory A term.
    assert sp.diff(first - flipped, beta) == 0


def test_a_missing_half_in_the_offset_fails_an_exact_readout_control():
    # A=(3+4i)/5, w=(1,0), beta=0. Correct 2*theta=arg A gives +1.
    ar, ai = sp.Rational(3, 5), sp.Rational(4, 5)
    correct = ar**2 + ai**2
    doubled_phase = ar * (ar**2 - ai**2) + ai * (2 * ar * ai)
    assert correct == 1
    assert doubled_phase == sp.Rational(3, 5) != correct


@pytest.mark.parametrize("side", [-1, 1])
def test_endpoint_limit_and_original_Q_readout_against_frozen_map(side):
    delta, k = sp.symbols("delta k", positive=True)
    ku = sp.symbols("ku", real=True)
    a = sp.Rational(1, 200)
    frozen = connection.endpoint_map(side * a, delta, k, ku)
    limit = frozen.applyfunc(lambda value: sp.limit(value, delta, 0, dir="+"))
    actual = phase.endpoint_limit(k, ku, side=side)
    assert _zero_matrix(actual - limit)
    assert sp.simplify(actual.det() - k / connection.RHO) == 0
    assert _zero_matrix(actual.inv()[0, :] - sp.Matrix([[sp.sqrt(a / k), 0]]))


def test_physical_Q_is_not_a_delta_dependent_chosen_projection():
    q, qu, ku = sp.symbols("q qu ku", real=True)
    k, delta, tau = sp.symbols("k delta tau", positive=True)
    a = sp.Rational(1, 200)
    actual = connection.endpoint_map(a, delta, k, ku)
    state = actual * sp.Matrix([q, qu])
    radius = sp.sqrt(a**2 + delta / 8)
    assert sp.simplify(sp.sqrt(radius / k) * state[0] - q) == 0
    physical_time_input = actual * sp.diag(1, tau)
    assert sp.simplify(physical_time_input.det() - tau * k / connection.RHO) == 0
    assert sp.simplify(physical_time_input[0, 0] - actual[0, 0]) == 0
    assert physical_time_input[0, 1] == 0
    # The k_u jet is required in the second normalized coordinate even though
    # the exact scalar Q readout uses only the first one.
    assert sp.simplify(sp.diff(actual[1, 0], ku)) != 0


def test_scalar_separation_retains_two_errors_and_the_kinetic_prefactor():
    f = Fraction
    error = f(4, 399)
    pair_gap = (2 - 2 * error) * f(9, 1280) * f(2, 5)
    assert pair_gap == f(237, 42560)
    assert pair_gap - f(1, 200) == f(121, 212800) > 0
    assert phase.calibration()["scalar_Q_gap_lower_over_eta"] == pair_gap
    assert audit.constants()["Q_pair_separation_per_eta"] == pair_gap
    a, k, yy = sp.symbols("a k Y", positive=True)
    loaded_w1 = yy / sp.sqrt(a)
    assert sp.simplify(sp.sqrt(a / k) * loaded_w1 - yy / sp.sqrt(k)) == 0


def test_bounded_errors_cannot_be_counted_once_or_assumed_to_vanish():
    # Norm-admissible comparison errors, not a new physical ODE fixture.
    # At aligned +/-I reference phases and w=e1, errors -e I and +e I
    # realize the full two-error loss; a claimed single-error bound is false.
    error = sp.Rational(4, 399)
    plus = sp.eye(2) - error * sp.eye(2)
    minus = -sp.eye(2) + error * sp.eye(2)
    difference = ((plus - minus) * sp.Matrix([1, 0]))[0]
    assert difference == 2 - 2 * error < 2 - error
    # Alternating a bounded error with zero satisfies the same uniform bound
    # and has no limit. The proof needs only the fixed comparison envelopes.
    envelopes = [sp.Rational(0), error]
    assert max(envelopes) == error and min(envelopes) == 0
    assert phase.calibration()["subsequential_reference_limits_asserted"] is False
    assert phase.calibration()["central_remainder_convergence_assumed"] is False


@pytest.mark.parametrize("offset", [0, sp.pi / 7])
def test_exact_parameter_sequence_inverts_the_asinh_clock(offset):
    sequence = phase.phase_sequences(2, offset)
    for sign in ("plus", "minus"):
        assert sp.simplify(phase.phase_value(sequence[f"delta_{sign}"])
                           - sequence[f"theta_{sign}"]) == 0
    assert sp.simplify(sequence["theta_minus"] - sequence["theta_plus"] - sp.pi / 2) == 0
    assert sp.simplify(sp.exp(2 * sp.I * sequence["theta_minus"])
                       + sp.exp(2 * sp.I * sequence["theta_plus"])) == 0


def test_sequence_monotonicity_limit_and_all_index_domain_bound():
    theta, a, rho = sp.symbols("theta a rho", positive=True)
    delta = 8 * a**2 / sp.sinh(theta / rho)**2
    positive_factor = 16 * a**2 * sp.cosh(theta / rho) / (rho * sp.sinh(theta / rho)**3)
    assert sp.simplify(sp.diff(delta, theta) + positive_factor) == 0
    assert sp.limit(delta, theta, sp.oo) == 0
    f = Fraction
    # For n>=2, offset>=0, rho<4/3 and pi>3: theta/rho>9/2.
    sinh_lower = f(9, 2) + f(9, 2)**3 / 6
    assert sinh_lower == f(315, 16)
    upper = 8 * f(1, 200)**2 / sinh_lower**2
    assert 0 < upper < f(1, 10**6)
    assert phase.calibration()["sequence_delta_upper"] == upper
    assert f(7, 4) < f(4, 3)**2


@pytest.mark.parametrize("bad", [True, False, 1, 0, -1, sp.Rational(5, 2), 2.0])
def test_sequence_index_cannot_bypass_the_stated_common_domain(bad):
    with pytest.raises((TypeError, ValueError)):
        phase.phase_sequences(bad)


@pytest.mark.parametrize("bad", [True, -1, 0, 1.0, sp.oo, sp.nan])
def test_literal_phase_delta_is_positive_and_exact(bad):
    with pytest.raises((TypeError, ValueError)):
        phase.phase_value(bad)


def test_declared_source_state_and_physical_scope_are_retained():
    data = loading.calibration()
    assert data["initial_data"].startswith("Q(-L)=Q_u(-L)=0")
    assert "tau is fixed" in data["physical_clock"]
    assert "punctured interval only" in data["loading_limit"]
    assert data["matter_source_response"] is False
    assert phase.calibration()["matter_source_or_EFT_verdict"] is False
    assert all(audit.checks().values())
    assert all(value == 0 for value in audit.phase_algebra()["identities"].values())
