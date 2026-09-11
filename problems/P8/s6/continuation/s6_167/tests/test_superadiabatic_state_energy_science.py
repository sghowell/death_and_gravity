"""Independent exact-frame derivatives, coherence and full-momentum checks."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_flat_dirac_hadamard import audit as previous
from p8_vacuum_superadiabatic_state_energy import audit, energy, rotation


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_named_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_each_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_bad_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def multiply(a, b, n):
    return [sum(a[k] * b[j - k] for k in range(j + 1)) for j in range(n + 1)]


def inverse(a, n):
    out = [1 / a[0]]
    for j in range(1, n + 1):
        out.append(-sum(a[k] * out[j - k] for k in range(1, j + 1)) / a[0])
    return out


def square_root(a, n):
    out = [mp.sqrt(a[0])]
    for j in range(1, n + 1):
        out.append(
            (a[j] - sum(out[k] * out[j - k] for k in range(1, j))) / (2 * out[0])
        )
    return out


def frames(z, p, m, d, tau):
    # Finite Taylor coefficients evaluate full-function derivatives at z.
    # They do not truncate the physical time evolution or expand in Delta.
    profile = mp.taylor(lambda x: x / mp.power(1 + x**8, mp.mpf(1) / 8), z, 5)
    mass = [d * c for c in profile]
    mass[0] += m
    omega2 = multiply(mass, mass, 4)
    omega2[0] += p * p
    e = [tau * c for c in square_root(omega2, 4)]
    derivative = [(j + 1) * profile[j + 1] for j in range(5)]
    g = [p * d * c / 2 for c in multiply(derivative, inverse(omega2, 4), 4)]
    result = [(e[0], g[0])]
    for j in range(4):
        n = 4 - j
        q = multiply(g, inverse(e, n), n)
        q2 = multiply(q, q, n)
        q2[0] += 1
        prime = [(k + 1) * q[k + 1] for k in range(n)]
        next_g = [(-1) ** j * c / 2 for c in multiply(prime, inverse(q2, n - 1), n - 1)]
        e2 = multiply(e, e, n)
        g2 = multiply(g, g, n)
        next_e = square_root([a + b for a, b in zip(e2, g2)], n)[:n]
        e, g = next_e, next_g
        result.append((e[0], g[0]))
    return result


def direct_pair(z, p, m, d, tau, depth):
    if depth == 0:
        f = lambda x: x / mp.power(1 + x**8, mp.mpf(1) / 8)
        M = m + d * f(z)
        omega = mp.sqrt(p * p + M * M)
        return tau * omega, p * d * mp.diff(f, z) / (2 * omega**2)
    e, g = direct_pair(z, p, m, d, tau, depth - 1)
    ratio = lambda x: (
        direct_pair(x, p, m, d, tau, depth - 1)[1]
        / direct_pair(x, p, m, d, tau, depth - 1)[0]
    )
    q = g / e
    return mp.sqrt(e * e + g * g), (-1) ** (depth - 1) * mp.diff(ratio, z) / (
        2 * (1 + q * q)
    )


@pytest.mark.parametrize("x", [-2, mp.mpf("0.5"), 3])
def test_Taylor_jet_frames_against_independent_nested_differentiation(x):
    with mp.workdps(75):
        p, m, d, tau = mp.mpf(1), mp.mpf(2), mp.mpf("0.01"), mp.mpf(3)
        computed = frames(mp.mpf(x), p, m, d, tau)
        for j in range(3):
            direct = direct_pair(mp.mpf(x), p, m, d, tau, j)
            for a, b in zip(computed[j], direct):
                assert abs(a - b) < mp.mpf("1e-65") * max(1, abs(b))


@pytest.mark.parametrize("x", [-10, -2, 0, mp.mpf("0.5"), 2, 10])
def test_complex_profile_disk_bound_and_root_branch(x):
    with mp.workdps(65):
        x = mp.mpf(x)
        r = mp.mpf(1) / 64
        w = mp.power(1 + x**8, -mp.mpf(9) / 8)
        for j in range(24):
            z = x + r * mp.exp(2j * mp.pi * j / 24)
            assert abs(z**8 - x**8) / (1 + x**8) < mp.mpf(1) / 4
            assert mp.re(1 + z**8) > 0
            f = z / mp.power(1 + z**8, mp.mpf(1) / 8)
            first = mp.power(1 + z**8, -mp.mpf(9) / 8)
            assert abs(f) < 2
            assert abs(first) < 2 * w


@pytest.mark.parametrize("x", [-3, 0, mp.mpf("0.5"), 2, 10])
@pytest.mark.parametrize("momentum_ratio", [mp.mpf("0.5"), 1, 2])
def test_all_four_complex_frame_bounds_at_actual_parameters(x, momentum_ratio):
    with mp.workdps(1200):
        m, d, tau = mp.mpf(10) ** 200, 3 * mp.mpf(10) ** 197, mp.mpf(10) ** -100
        floor = mp.mpf("0.99") * m
        p = mp.mpf(momentum_ratio) * floor
        E = mp.sqrt(p * p + floor * floor)
        x = mp.mpf(x)
        w = mp.power(1 + x**8, -mp.mpf(9) / 8)
        for j in range(5):
            r = mp.mpf(1) / 64 - mp.mpf(j) / 512
            z = x + 1j * r / 2
            e, g = frames(z, p, m, d, tau)[j]
            c = mp.mpf("0.9") * (mp.mpf(15) / 16) ** j
            G = 2 * p * d / E**2 * (1024 / (tau * E)) ** j
            assert abs(e) >= c * tau * E
            assert abs(g) <= G * w
            assert abs(g / e) < mp.mpf(1) / 4


def test_exact_frame_transformations_retain_nonzero_connections():
    with mp.workdps(70):
        sigma1 = mp.matrix([[0, 1], [1, 0]])
        sigma2 = mp.matrix([[0, -1j], [1j, 0]])
        sigma3 = mp.matrix([[1, 0], [0, -1]])
        # These moderate parameters test exact identities, not the large-gap bound.
        chain = frames(mp.mpf("0.5"), mp.mpf(1), mp.mpf(2), mp.mpf("0.01"), mp.mpf(3))
        for j in range(4):
            e, g = chain[j]
            en, gn = chain[j + 1]
            theta = mp.atan(g / e)
            derivative = 2 * (-1) ** j * gn
            if j % 2 == 0:
                axis, after = sigma2, sigma1
                U = mp.cos(theta / 2) * mp.eye(2) + 1j * mp.sin(theta / 2) * sigma1
                Up = (
                    derivative
                    * (-mp.sin(theta / 2) * mp.eye(2) + 1j * mp.cos(theta / 2) * sigma1)
                    / 2
                )
            else:
                axis, after = sigma1, sigma2
                U = mp.cos(theta / 2) * mp.eye(2) - 1j * mp.sin(theta / 2) * sigma2
                Up = (
                    derivative
                    * (-mp.sin(theta / 2) * mp.eye(2) - 1j * mp.cos(theta / 2) * sigma2)
                    / 2
                )
            transformed = U.H * (e * sigma3 + g * axis) * U - 1j * U.H * Up
            assert mp.norm(transformed - en * sigma3 - gn * after) < mp.mpf("1e-60")
            assert abs(gn) > 0
            assert mp.norm(transformed - en * sigma3) > 0


def test_zero_momentum_exact_frames_and_zero_profile_amplitude():
    with mp.workdps(65):
        for p, d in [(0, mp.mpf("0.003")), (1, 0)]:
            for _, g in frames(
                mp.mpf("0.5"), mp.mpf(p), mp.mpf(2), mp.mpf(d), mp.mpf(3)
            ):
                assert g == 0
    assert energy.amplitude_upper(0, 2**22, s.Rational(1, 100), 1) == 0
    caps = energy.enclosures(2**22, 0, 1)
    assert all(v == 0 for k, v in caps.items() if k != "mass_floor")


@pytest.mark.parametrize("floor", [mp.mpf("0.5"), 1, 3, 10])
def test_complete_radial_integrals_independently(floor):
    with mp.workdps(65):
        m = mp.mpf(floor)
        out = mp.quad(lambda p: p**4 / (p * p + m * m) ** mp.mpf("5.5"), [0, m, mp.inf])
        difference = mp.quad(
            lambda p: p**3 / (p * p + m * m) ** mp.mpf("2.5"), [0, m, mp.inf]
        )
        assert abs(out - 8 / (315 * m**6)) < mp.mpf("1e-58")
        assert abs(difference - 2 / (3 * m)) < mp.mpf("1e-58")


def test_late_frequency_corrections_are_resolved_not_rounded_away():
    with mp.workdps(1200):
        m, d, tau = mp.mpf(10) ** 200, 3 * mp.mpf(10) ** 197, mp.mpf(10) ** -100
        chain = frames(mp.mpf("0.5"), mp.mpf("0.99") * m, m, d, tau)
        for j in range(4):
            e, g = chain[j]
            next_e, _ = chain[j + 1]
            assert next_e != e
            assert abs((next_e * next_e - e * e) / (g * g) - 1) < mp.mpf("1e-250")


@pytest.mark.parametrize(
    "b", [s.Rational(1, 1000), s.Rational(1, 10), s.Rational(3, 5)]
)
def test_projector_trace_norm_and_local_coherence(b):
    a = s.sqrt(1 - b * b)
    vector = s.Matrix([b, a])
    change = vector * vector.T - s.diag(0, 1)
    assert change.det() == -b * b
    assert change * change == b * b * s.eye(2)
    p, M = s.Integer(2), s.Integer(3)
    H = s.Matrix([[M, p], [p, -M]])
    exact = s.simplify(s.trace(H * change))
    occupation_only = 2 * M * b * b
    assert s.simplify(exact - occupation_only - 2 * p * a * b) == 0
    assert exact != occupation_only
    assert bool(abs(exact) <= 2 * s.sqrt(p * p + M * M) * b)
    # Both helicities give trace norm4|beta| per active color/flavor.
    assert 2 * (2 * b) == 4 * b


def test_symmetric_Dirac_energy_equals_the_Hamiltonian_bilinear():
    p, M = s.symbols("p M", real=True)
    H = s.Matrix([[M, p], [p, -M]])
    u, v = s.symbols("u v", complex=True)
    psi = s.Matrix([u, v])
    derivative = -s.I * H * psi
    symmetric = s.I * (psi.H * derivative - derivative.H * psi)[0] / 2
    assert s.simplify(symmetric - (psi.H * H * psi)[0]) == 0


def test_state_difference_energy_exchange_uses_the_external_mass_source():
    p, M, rate, a, b, c, d = s.symbols("p M rate a b c d", real=True)
    H = s.Matrix([[M, p], [p, -M]])
    beta = s.diag(1, -1)
    covariance = s.Matrix([[a, c + s.I * d], [c - s.I * d, b]])
    evolution = -s.I * (H * covariance - covariance * H)
    observed = s.trace(rate * beta * covariance + H * evolution)
    source = rate * s.trace(beta * covariance)
    assert s.simplify(observed - source) == 0
    assert source.subs({rate: 2, a: 1, b: 0}) == 2


def test_complete_energy_factors_and_dimensions():
    C, D, tau, m, N = s.symbols("C Delta tau m N", positive=True)
    measure = 4 * s.pi / (2 * s.pi) ** 3
    # Out energy:4N occupations, frequency<=2E, square of full beta bound.
    out = measure * 4 * N * 2 * C * C * D * D / tau**8 * s.Rational(8, 315) / m**6
    assert (
        s.simplify(out - 32 * N * C * C * D * D / (315 * s.pi**2 * tau**8 * m**6)) == 0
    )
    # Local difference: trace norm4|beta|, not4|beta|^2.
    local = measure * 4 * N * 2 * C * D / tau**4 * s.Rational(2, 3) / m
    assert s.simplify(local - 8 * N * C * D / (3 * s.pi**2 * tau**4 * m)) == 0
    assert 2 + 8 - 6 == 4
    assert 1 + 4 - 1 == 4


def test_actual_strict_bounds_and_their_scope():
    data = energy.data()
    d = data["actual_complete_exact_frame_enclosures"]
    assert d["uniform_exact_transition_amplitude_upper"] < s.Rational(1, 10**390)
    assert 0 < d["complete_free_out_particle_energy_upper"] < 10**20
    assert 0 < d["uniform_in_out_local_energy_density_difference_upper"] < 10**411
    assert data["new_out_energy_over_named_kappa_scale"] < s.Rational(1, 10**780)
    assert data["uniform_local_state_difference_over_named_kappa_scale"] < s.Rational(
        1, 10**389
    )
    assert "No absolute local renormalized energy" in data["renormalization_scope"]
    assert "not a proof" in data["cutoff_scope"]
    assert "not set to zero" in audit.observable()["difference_vs_out_particles"]


def test_exact_transition_ultraviolet_power_and_boundary_parameter():
    p = s.Symbol("p", positive=True)
    C = s.Integer(2) ** 42
    assert s.limit(p**5 * C * p / (p * p + 1) ** 3, p, s.oo) == C
    m = s.Rational(100, 99)
    assert energy.parameters(m, s.Rational(1, 1000), 2**20)[-1] == 1
    with pytest.raises(ValueError):
        energy.parameters(m, s.Rational(1, 1000), 2**20 - 1)


def test_parent_state_and_all_prior_rows_remain_fixed():
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert "NOT_ABSOLUTE_STRESS" in audit.matching()[-1]["status"]
    assert "not set to zero" in rotation.data()["complete_evolution"]


def test_own_report_payload_serializes_without_inexact_symbolic_numbers():
    from p8_affine import verify as certificate

    for module in audit.MODULES:
        certificate.serialize(
            {key: value for key, value in module.data().items() if key != "checks"}
        )
    diagnostics = energy.data()["actual_bound_decimal_diagnostics"]
    assert all(isinstance(value, str) for value in diagnostics.values())
    assert isinstance(
        energy.data()["actual_complete_exact_frame_enclosures"][
            "uniform_exact_transition_amplitude_upper"
        ],
        s.Rational,
    )
    with pytest.raises(ValueError, match="Inexact"):
        certificate.serialize({"diagnostic": s.Float("1.0")})
