"""Independent derivative, convention and norm-limit diagnostics."""

from math import cos, expm1, hypot, log1p, sin, sqrt

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_flat_dirac_hadamard import audit, scattering, symbols
from p8_vacuum_flat_dirac_production import audit as previous


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_exact_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_rejected_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("n", [0, 1, 2, 3, 4, 7, 9, 12])
@pytest.mark.parametrize("x", [-10, -3, -2, 2, 3, 10])
def test_independent_high_order_real_derivatives_and_tail_caps(n, x):
    with mp.workdps(90):
        x = mp.mpf(x)
        f = lambda z: z / mp.power(1 + z**8, mp.mpf(1) / 8)
        observed = mp.diff(f, x, n)
        P = s.Poly(symbols.polynomial(n), symbols.X)
        numerator = sum(mp.mpf(int(c)) * x ** powers[0] for powers, c in P.terms())
        expected = numerator / mp.power(1 + x**8, n + mp.mpf(1) / 8)
        assert abs(observed - expected) < mp.mpf("1e-75") * max(1, abs(expected))
        tail = observed - (mp.sign(x) if n == 0 else 0)
        upper = 5 * mp.factorial(n) * 4**n / abs(x) ** (8 + n)
        assert abs(tail) <= upper
        reflected = mp.diff(f, -x, n)
        assert abs(reflected - (-1) ** (n + 1) * observed) < mp.mpf("1e-75") * max(
            1, abs(observed)
        )


@pytest.mark.parametrize("n", [0, 1, 2, 4, 8, 12])
@pytest.mark.parametrize("x", [s.Rational(-3, 2), 0, s.Rational(1, 2), 2])
def test_compact_derivative_polynomial_cap(n, x):
    P = symbols.polynomial(n)
    actual = abs(P.subs(symbols.X, x)) / (1 + x**8) ** (s.Rational(1, 8) + n)
    assert bool(actual <= symbols.compact_coefficient(n))


@pytest.mark.parametrize("x", [2, 3, 10])
def test_complex_disk_branch_bound_and_cauchy_coefficients(x):
    with mp.workdps(75):
        x = mp.mpf(x)
        radius = x / 4
        sample_count = 192
        samples = []
        for j in range(sample_count):
            angle = 2 * mp.pi * j / sample_count
            z = x + radius * mp.exp(1j * angle)
            w = z**-8
            observed = mp.power(1 + w, -mp.mpf(1) / 8) - 1
            assert abs(w) < mp.mpf(1) / 2
            assert abs(observed) < 5 * x**-8
            samples.append(observed)
        for n in (0, 1, 3, 8):
            # Fourier coefficient on the complex circle gives f^(n)/n!.
            coefficient = (
                sum(
                    value * mp.exp(-2j * mp.pi * n * j / sample_count)
                    for j, value in enumerate(samples)
                )
                / sample_count
                / radius**n
            )
            direct = mp.diff(
                lambda z: mp.power(1 + z**-8, -mp.mpf(1) / 8) - 1, x, n
            ) / mp.factorial(n)
            assert abs(coefficient - direct) < mp.mpf("1e-65")


def test_wrong_eighth_root_is_not_the_real_profile_branch():
    with mp.workdps(60):
        x = mp.mpf(3)
        correct = mp.power(1 + x**-8, -mp.mpf(1) / 8)
        wrong = correct * mp.exp(2j * mp.pi / 8)
        assert abs(wrong - 1) > 5 * x**-8
        assert abs(correct - x / mp.power(1 + x**8, mp.mpf(1) / 8)) < mp.mpf("1e-55")


@pytest.mark.parametrize("n", [0, 1, 2, 4, 9])
def test_physical_time_rescaling_of_the_tail(n):
    d, tau, T = s.Rational(3, 1000), s.Rational(1, 10), 2
    dimensionless = symbols.tail_coefficient(n) / (s.Rational(T) / tau) ** (8 + n)
    assert symbols.mass_derivative_tail(n, d, tau, T) == d * tau ** (-n) * dimensionless


def test_source_clifford_and_Dirac_operator_dictionary_independently():
    d = scattering.data()
    physical, source = d["physical_gamma"], d["source_frame_gamma"]
    gp, gs = s.diag(1, -1, -1, -1), s.diag(-1, 1, 1, 1)
    for a in range(4):
        for b in range(4):
            assert physical[a] * physical[b] + physical[b] * physical[a] == 2 * gp[
                a, b
            ] * s.eye(4)
            assert source[a] * source[b] + source[b] * source[a] == 2 * gs[
                a, b
            ] * s.eye(4)
    derivatives = s.symbols("d0:4")
    mass = s.Symbol("M")
    physical_operator = sum(
        (s.I * physical[a] * derivatives[a] for a in range(4)), s.zeros(4)
    ) - mass * s.eye(4)
    source_operator = sum(
        (gs[a, a] * source[a] * derivatives[a] for a in range(4)), s.zeros(4)
    ) - mass * s.eye(4)
    assert source_operator == physical_operator
    # The unconverted source mass sign would give a different operator.
    assert source_operator + 2 * mass * s.eye(4) != physical_operator


def test_spectral_labels_exchange_but_occupied_projector_is_preserved():
    d = scattering.data()
    H = d["physical_Hamiltonian"]
    variables = {v: {"M": 3, "p1": 1, "p2": 2, "p3": 4}[str(v)] for v in H.free_symbols}
    H = H.subs(variables)
    source = d["source_Hamiltonian"].subs(variables)
    omega = s.sqrt(30)
    occupied = (s.eye(4) - H / omega) / 2
    source_positive = (s.eye(4) + source / omega) / 2
    assert s.simplify(occupied - source_positive) == s.zeros(4)
    assert s.simplify(occupied * occupied - occupied) == s.zeros(4)
    assert s.trace(occupied) == 2
    assert H.det() != 0


@pytest.mark.parametrize("T", [1, 2, 4, 8])
def test_tail_integral_and_projector_factor_without_mode_sampling(T):
    with mp.workdps(70):
        T = mp.mpf(T)
        d = mp.mpf("0.009")
        observed = mp.quad(
            lambda x: -d * mp.expm1(-mp.log1p(x**-8) / 8), [T, 2 * T, mp.inf]
        )
        bound = d / (56 * T**7)
        assert 0 < observed <= bound
        exact = scattering.norm_caps(s.Rational(9, 1000), 1, int(T))
        assert abs(
            mp.mpf(str(exact["one_particle_Moller_operator_norm_error"].evalf(75)))
            - bound
        ) < mp.mpf("1e-65")
        assert (
            exact["transported_asymptotic_projector_norm_error"]
            == 2 * exact["one_particle_Moller_operator_norm_error"]
        )


def multiply(a, b):
    return (
        a[0] * b[0] + a[1] * b[2],
        a[0] * b[1] + a[1] * b[3],
        a[2] * b[0] + a[3] * b[2],
        a[2] * b[1] + a[3] * b[3],
    )


def adjoint(a):
    return a[0].conjugate(), a[2].conjugate(), a[1].conjugate(), a[3].conjugate()


def difference(a, b):
    return tuple(x - y for x, y in zip(a, b))


def operator_norm(a):
    trace = sum(abs(v) ** 2 for v in a)
    determinant = abs(a[0] * a[3] - a[1] * a[2]) ** 2
    return sqrt((trace + sqrt(max(0.0, trace * trace - 4 * determinant))) / 2)


IDENTITY = (1 + 0j, 0j, 0j, 1 + 0j)


def tail_unitary(p, T, steps):
    """Exact-asymptotic interaction picture; finite-interval RK4 diagnostic."""
    delta = 0.009
    mass = 1 + delta
    omega = hypot(p, mass)
    beta = (1, 0, 0, -1)

    def rhs(t, V):
        # Integrate V-I so a tiny diagonal change is not repeatedly added
        # to one. This avoids cumulative cancellation at machine precision.
        V = tuple(v + i for v, i in zip(V, IDENTITY))
        theta = omega * (t - T)
        c, q = cos(theta), 1j * sin(theta) / omega
        free = (c + q * mass, q * p, q * p, c - q * mass)
        rotated = multiply(multiply(free, beta), adjoint(free))
        dm = delta * expm1(-log1p(t**-8) / 8)
        return tuple(-1j * dm * v for v in multiply(rotated, V))

    h = T / steps
    V = (0j, 0j, 0j, 0j)
    for j in range(steps):
        t = T + j * h
        a = rhs(t, V)
        b = rhs(t + h / 2, tuple(v + h * k / 2 for v, k in zip(V, a)))
        c = rhs(t + h / 2, tuple(v + h * k / 2 for v, k in zip(V, b)))
        d = rhs(t + h, tuple(v + h * k for v, k in zip(V, c)))
        V = tuple(
            v + h * (ka + 2 * kb + 2 * kc + kd) / 6
            for v, ka, kb, kc, kd in zip(V, a, b, c, d)
        )
    V = tuple(v + i for v, i in zip(V, IDENTITY))
    projector = (
        (1 - mass / omega) / 2,
        -p / (2 * omega),
        -p / (2 * omega),
        (1 + mass / omega) / 2,
    )
    moved = multiply(multiply(adjoint(V), projector), V)
    return V, projector, moved


@pytest.mark.parametrize("p", [0, 1, 3])
@pytest.mark.parametrize("T", [2, 4, 8, 16])
def test_finite_endpoint_operator_and_projector_errors_with_refinement(p, T):
    V, P, moved = tail_unitary(p, T, 800)
    coarse, _, _ = tail_unitary(p, T, 400)
    cap = 0.009 / (56 * T**7) * (1 - 2**-7)
    assert operator_norm(difference(V, coarse)) < max(2e-15, cap * 1e-5)
    assert operator_norm(difference(multiply(adjoint(V), V), IDENTITY)) < 2e-14
    assert operator_norm(difference(V, IDENTITY)) <= cap + 2e-15
    assert operator_norm(difference(moved, P)) <= 2 * cap + 4e-15
    if p == 0:
        assert operator_norm(difference(moved, P)) < 2e-14


def test_exact_commuting_tail_rotation_and_nonzero_wave_error():
    with mp.workdps(60):
        T = mp.mpf(4)
        integral = mp.quad(
            lambda t: mp.mpf("0.009") * mp.expm1(-mp.log1p(t**-8) / 8), [T, 2 * T]
        )
        expected = complex(mp.exp(-1j * integral))
        V, _, _ = tail_unitary(0, 4, 800)
        assert abs(V[0] - expected) < 2e-15
        assert 0 < abs(expected - 1) < 0.009 / (56 * 4**7)


def test_actual_operator_errors_and_no_stress_inference():
    d = scattering.data()
    caps = d["actual_unit_time_norm_enclosures"]
    assert caps["one_particle_Moller_operator_norm_error"] == s.Rational(
        3, 56 * 10**603
    )
    assert caps["transported_asymptotic_projector_norm_error"] == s.Rational(
        3, 28 * 10**603
    )
    assert "not a second-quantized global Fock" in d["norm_scope"]
    assert "numerical stress" in d["norm_scope"]
    assert "no numerical stress" in audit.application()["quantitative_boundary"]


def test_parent_scope_and_exact_application_are_preserved():
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert len(audit.matching()) == len(previous.matching()) + 1
    assert audit.validate_application(audit.application())
    assert "HADAMARD_FOR_SPECIFIED_FREE" in audit.matching()[-1]["status"]
    assert "external theorem" in audit.application()["quantitative_boundary"]
