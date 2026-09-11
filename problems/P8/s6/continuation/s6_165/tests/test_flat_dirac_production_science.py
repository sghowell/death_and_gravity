"""Independent mode ODE, full-profile and complete energy-integral diagnostics."""

from math import cos, hypot, sin

import mpmath as mp
import pytest
import sympy as s
from p8_exceptional_vacuum import analytic
from p8_vacuum_clock_transparent_map import audit as previous
from p8_vacuum_flat_dirac_production import (
    audit,
    calibration,
    energy,
    mode,
    profile,
    transition,
)


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_named_exact_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_each_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_rejected_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("p,m", [(0, 2), (s.Rational(1, 3), 2), (3, 5)])
def test_independent_eigenprojectors_and_both_helicities(p, m):
    sigma1 = s.Matrix([[0, 1], [1, 0]])
    sigma3 = s.diag(1, -1)
    omega = s.sqrt(p * p + m * m)
    for helicity in (-1, 1):
        H = helicity * p * sigma1 + m * sigma3
        occupied = (s.eye(2) - H / omega) / 2
        assert s.simplify(occupied * occupied - occupied) == s.zeros(2)
        assert s.simplify(H * occupied + omega * occupied) == s.zeros(2)
        assert s.trace(occupied) == 1
        assert occupied.H == occupied


def test_independent_frequency_derivatives_before_absolute_bounds():
    t, p = s.symbols("t p", real=True)
    M = s.Function("M")(t)
    omega = s.sqrt(p * p + M * M)
    assert s.simplify(s.diff(omega, t) - M * s.diff(M, t) / omega) == 0
    exact_second = M * s.diff(M, t, 2) / omega + p * p * s.diff(M, t) ** 2 / omega**3
    assert s.simplify(s.diff(omega, t, 2) - exact_second) == 0


def test_full_profile_variation_by_independent_quadrature():
    with mp.workdps(60):
        xstar = mp.power(mp.mpf(7) / 10, mp.mpf(1) / 8)
        first = lambda x: mp.power(1 + x**8, -mp.mpf(9) / 8)
        second = lambda x: -9 * x**7 / mp.power(1 + x**8, mp.mpf(17) / 8)
        third = lambda x: (
            -9 * x**6 * (7 - 10 * x**8) / mp.power(1 + x**8, mp.mpf(25) / 8)
        )
        pieces = [-mp.inf, -xstar, 0, xstar, mp.inf]
        assert abs(mp.quad(first, pieces) - 2) < mp.mpf("1e-55")
        assert abs(mp.quad(lambda x: abs(first(x) * second(x)), pieces) - 1) < mp.mpf(
            "1e-55"
        )
        observed = mp.quad(lambda x: abs(third(x)), pieces)
        assert abs(observed - 4 * abs(second(xstar))) < mp.mpf("1e-55")
        assert observed < 36
        assert 0 < mp.quad(lambda x: first(x) ** 3, pieces) < 2


@pytest.mark.parametrize("T", [1, 2, 3, 10, 100])
def test_mass_tail_and_finite_interval_mixing_tail(T):
    with mp.workdps(70):
        T = mp.mpf(T)
        f = lambda x: x / mp.power(1 + x**8, mp.mpf(1) / 8)
        derivative = lambda x: mp.power(1 + x**8, -mp.mpf(9) / 8)
        tail = 1 - f(T)
        assert 0 < tail <= 1 / (8 * T**8)
        assert abs(mp.quad(derivative, [T, mp.inf]) - tail) < mp.mpf("1e-60")
        # Two time tails: their integrated mixing is <= eta * (1-s(T)).
        eta = mp.mpf("0.003")
        assert eta * tail <= eta / (8 * T**8)
        assert abs(f(-T) + f(T)) < mp.mpf("1e-65")


@pytest.mark.parametrize("floor", [s.Rational(1, 2), 1, 3, 10])
def test_radial_integrals_independently_over_full_momentum_range(floor):
    with mp.workdps(60):
        m = mp.mpf(str(s.N(floor, 65)))
        first = mp.quad(
            lambda p: p**4 / (p * p + m * m) ** mp.mpf("3.5"), [0, m, mp.inf]
        )
        higher = mp.quad(
            lambda p: p**8 / (p * p + m * m) ** mp.mpf("5.5"), [0, m, mp.inf]
        )
        assert abs(first - 1 / (5 * m * m)) < mp.mpf("1e-55")
        assert abs(higher - 1 / (9 * m * m)) < mp.mpf("1e-55")


def test_energy_composition_retains_cross_term_allowance_and_dirac_counts():
    d, m, tau, N = s.symbols("d m tau N", positive=True)
    # Spherical density of modes * 2 helicities * particle/antiparticle.
    measure = 4 * s.pi / (2 * s.pi) ** 3 * 2 * 2 * N
    # omega_out <= 2E, and (A+B)^2 <= 2 A^2 + 2 B^2.
    integrated = (
        measure * 2 * 2 * (25 * d * d / (tau**4 * 5 * m * m) + d**6 / (25 * 9 * m * m))
    )
    expected = 8 * N / s.pi**2 * (5 * d * d / (tau**4 * m * m) + d**6 / (225 * m * m))
    assert s.simplify(integrated - expected) == 0
    assert 2 * 2 * 6 == 24
    assert 12 * 3 == 36
    assert (2 + 12) * 3 == 42
    assert (s.Integer(1) + 1) ** 2 > s.Integer(1) ** 2 + 1


@pytest.mark.parametrize("p", [0, 1, 3, 10**6])
def test_constant_mass_has_no_transition_or_energy(p):
    result = transition.amplitude_terms(p, 10, 0, 1)
    assert all(value == 0 for value in result.values())
    e = energy.enclosure(10, 0, 1)
    assert e["complete_free_out_particle_energy_density_upper"] == 0
    assert e["uniform_exact_transition_amplitude_upper"] == 0


def test_zero_momentum_exact_mixing_and_ultraviolet_power():
    assert (
        transition.amplitude_terms(0, 10, s.Rational(1, 100), 1)[
            "exact_transition_amplitude_upper"
        ]
        == 0
    )
    p = s.Symbol("p", positive=True)
    A = 5 * p * s.Rational(1, 100) / (p * p + 100) ** 2
    B = p**3 * s.Rational(1, 100) ** 3 / (5 * (p * p + 100) ** 3)
    assert s.limit(p**3 * A, p, s.oo) == s.Rational(1, 20)
    assert s.limit(p**3 * B, p, s.oo) == s.Rational(1, 5000000)


@pytest.mark.parametrize(
    "eta", [s.Rational(1, 1000), s.Rational(1, 500), s.Rational(1, 100)]
)
def test_complete_nonzero_tail_for_exact_constant_phase_rotation(eta):
    with mp.workdps(70):
        x = mp.mpf(int(eta.p)) / int(eta.q)
        # For constant phase, the exact off-diagonal evolution is sin(eta),
        # not its first term eta. Every higher odd order is included.
        difference = abs(mp.sin(x) - x)
        assert 0 < difference < x**3 / 5
        assert mp.sinh(x) - x < x**3 / 5
        assert abs(mp.cos(x) ** 2 + mp.sin(x) ** 2 - 1) < mp.mpf("1e-65")
        assert transition.mixing_tail(eta) == eta**3 / 5


def rk4_mode(p, tau, amplitude=0.009, steps=8000, endpoint=20.0):
    """Finite-interval diagnostic, not an enclosure or actual-scale integrator."""

    def rhs(x, state):
        ar, ai, br, bi, phase, _, _ = state
        shape = x / (1 + x**8) ** 0.125
        derivative = (1 + x**8) ** -1.125
        M = 1 + amplitude * shape
        omega = hypot(p, M)
        rate = p * amplitude * derivative / (2 * omega * omega)
        c, sn = cos(2 * phase), sin(2 * phase)
        return (
            -rate * (c * br - sn * bi),
            -rate * (c * bi + sn * br),
            rate * (c * ar + sn * ai),
            rate * (c * ai - sn * ar),
            tau * omega,
            rate * c,
            -rate * sn,
        )

    h = 2 * endpoint / steps
    state = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for j in range(steps):
        x = -endpoint + j * h
        a = rhs(x, state)
        b = rhs(x + h / 2, [v + h * k / 2 for v, k in zip(state, a)])
        c = rhs(x + h / 2, [v + h * k / 2 for v, k in zip(state, b)])
        d = rhs(x + h, [v + h * k for v, k in zip(state, c)])
        state = [
            v + h * (ka + 2 * kb + 2 * kc + kd) / 6
            for v, ka, kb, kc, kd in zip(state, a, b, c, d)
        ]
    alpha = complex(state[0], state[1])
    beta = complex(state[2], state[3])
    first = complex(state[5], state[6])
    return alpha, beta, first


@pytest.mark.parametrize("p", [0.25, 1, 3])
@pytest.mark.parametrize("tau", [0.5, 1, 2])
def test_independent_mode_solution_refines_and_retains_higher_transitions(p, tau):
    alpha, beta, first = rk4_mode(p, tau, steps=16000)
    _, coarse, _ = rk4_mode(p, tau, steps=8000)
    assert abs(beta - coarse) < 1e-10
    assert abs(abs(alpha) ** 2 + abs(beta) ** 2 - 1) < 1e-10
    eta = p * 0.009 / (p * p + 0.99**2)
    assert 0 < abs(beta - first) <= eta**3 / 5 + 2e-12
    E = hypot(p, 0.99)
    A, B = 5 * p * 0.009 / (tau * tau * E**4), eta**3 / 5
    interval_tail = eta / (8 * 20**8)
    assert abs(first) <= A + interval_tail + 2e-12
    assert abs(beta) <= A + B + interval_tail + 2e-12


@pytest.mark.parametrize("p,tau", [(0.25, 0.5), (1, 1), (3, 2)])
def test_opposite_mass_switch_has_same_probability(p, tau):
    _, plus, _ = rk4_mode(p, tau, amplitude=0.009)
    _, minus, _ = rk4_mode(p, tau, amplitude=-0.009)
    assert abs(abs(plus) ** 2 - abs(minus) ** 2) < 1e-12


def test_zero_momentum_numerical_control():
    alpha, beta, first = rk4_mode(0, 1, steps=1000)
    assert alpha == 1
    assert beta == first == 0


def test_actual_parameters_and_all_energy_terms_without_extreme_oscillation_sampling():
    d = calibration.data()
    enclosure = d["complete_quadratic_transition_and_out_energy_enclosure"]
    assert d["quadratic_subsystem_mean_mass"] == 10**200
    assert d["same_profile_time_scale"] == s.Rational(1, 10**100)
    assert d["profile_amplitude_rational_upper"] == 3 * 10**197
    first = enclosure["first_transition_part_of_energy_upper"]
    tail = enclosure["higher_transition_part_of_energy_upper"]
    rho = enclosure["complete_free_out_particle_energy_density_upper"]
    assert 0 < first < 10**397
    assert 0 < tail < 10**784
    assert rho == first + tail
    assert rho < 10**784
    assert rho / analytic.KAPPA < s.Rational(1, 10**16)
    assert enclosure["uniform_exact_transition_amplitude_upper"] < s.Rational(1, 10**9)
    assert "NOT a relative error" in d["reference_scale_scope"]
    assert "bounce" in d["reference_scale_scope"]


def test_free_state_scope_and_unchanged_parent_obligations():
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert len(audit.matching()) == len(previous.matching()) + 1
    assert "NOT_INTERACTING_OR_CURVED" in audit.matching()[-1]["status"]
    assert "global infinite-volume Fock unitary" in mode.data()["state"]
    assert "Hadamard" in mode.data()["scope"]
    assert "transient renormalized local stress" in energy.data()["scope"]
    assert profile.data()["first_derivative_L1"] == 2
    assert all(
        value is True
        for key, value in audit.controls().items()
        if key != "rejected_inputs"
    )
