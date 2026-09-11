"""Independent complete-source, constraint, covariance and profile fixtures."""

import mpmath as mp
import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_vacuum_causal_noise import audit, energy, noise, source, vacuum
from scipy.integrate import solve_ivp


def mpq(value):
    value = s.Rational(value)
    return mp.mpf(int(value.p)) / int(value.q)


def regular_A(u, x, n=1024):
    den = x**n + (1 - x) ** n
    T = x**n / den
    return (
        (x - 1)
        / (1 + u * u) ** 3
        * (x ** (n - 2) / den + (1 - T) * n * mp.exp(-n * x * x))
    )


def regular_A_X(u, x, n=1024):
    den = x**n + (1 - x) ** n
    denx = n * (x ** (n - 1) - (1 - x) ** (n - 1))
    T = x**n / den
    Tx = n * x ** (n - 1) * (1 - x) ** (n - 1) / den**2
    exp = mp.exp(-n * x * x)
    f = x ** (n - 2) / den + (1 - T) * n * exp
    fx = (
        (n - 2) * x ** (n - 3) / den
        - x ** (n - 2) * denx / den**2
        - Tx * n * exp
        - (1 - T) * 2 * n * n * x * exp
    )
    return (f + (x - 1) * fx) / (1 + u * u) ** 3


def literal_source(jet, k0, kappa, n=1024):
    phi = jet[0]
    g = jet[1:5]
    h = [[mp.mpf(0) for _ in range(4)] for _ in range(4)]
    at = 5
    for i in range(4):
        for j in range(i, 4):
            h[i][j] = h[j][i] = jet[at]
            at += 1
    sign = [1, -1, -1, -1]
    Y = sum(sign[i] * g[i] ** 2 for i in range(4))
    box = sum(sign[i] * h[i][i] for i in range(4))
    Z = sum(
        sign[i] * sign[j] * g[i] * g[j] * h[i][j] for i in range(4) for j in range(4)
    )
    u = phi / mp.sqrt(k0)
    x = Y / k0
    a = regular_A(u, x, n) / k0
    ap = -6 * u * regular_A(u, x, n) / (1 + u * u) / k0 ** mp.mpf("1.5")
    ay = regular_A_X(u, x, n) / k0**2
    R = 1 + Y * Y * a / kappa
    H = 4 * u / (1 + u * u)
    common = a * (Y * box - Z) - 3 * Y * Y * a * H / mp.sqrt(kappa)
    common += 3 * Y * Y * a * (Y * Y * ap / 4 + (2 * a + Y * ay) * Z / 2) / (kappa * R)
    return [gi * common / kappa for gi in g]


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_exact_residual(name):
    v = audit.residuals()[name]
    assert all(x == 0 for x in (list(v) if isinstance(v, s.MatrixBase) else [v])), name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_unsupported_scope_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("name", list(audit.packets()))
def test_all_continuous_proof_packets(name):
    assert all(audit.packets()[name]["gates"].values())


@pytest.mark.parametrize("n", (4, 8, 1024))
def test_full_complex_switch_and_removable_coefficient(n):
    with mp.workdps(90):
        rho = mp.mpf(1) / (16 * n)
        for i in range(12):
            x = rho * mp.exp(2j * mp.pi * i / 12)
            u = mp.mpf("0.25") * mp.exp(2j * mp.pi * (i + 3) / 12)
            den = x**n + (1 - x) ** n
            T = x**n / den
            assert abs(den) > mp.mpf(7) / 8
            assert abs(T) <= abs(x) ** 2
            assert abs(regular_A(u, x, n)) < 20 * n
            inner = x / 2
            actual_X = mp.diff(lambda y, u=u: regular_A(u / 2, y, n), inner)
            assert abs(actual_X) < 640 * n * n
            assert abs(actual_X - regular_A_X(u / 2, inner, n)) < mp.mpf("1e-60")
            assert (
                abs(mp.diff(lambda t, inner=inner: regular_A(t, inner, n), u / 2))
                < 160 * n
            )
        assert regular_A(mp.mpf(0), mp.mpf(0), n) == -n


@pytest.mark.parametrize("ratio", (1, 2, 100))
@pytest.mark.parametrize("hierarchy", ("minimum", "actual"))
def test_literal_full_source_on_complex_jet_domain(ratio, hierarchy):
    # The minimal hierarchy is an inequality fixture, not a new physical parent.
    with mp.workdps(100):
        k0 = mp.mpf(512 * 1024) if hierarchy == "minimum" else mp.mpf("1e800")
        k = ratio * k0
        upper = mp.mpf(20000 * 1024) / (k * k0)
        for fixture in range(4):
            jet = [
                mp.mpf("1.8") * mp.exp(2j * mp.pi * (i + fixture) / 17)
                for i in range(15)
            ]
            assert max(abs(x) for x in literal_source(jet, k0, k)) < upper


@pytest.mark.parametrize("order", (1, 2))
def test_full_multijet_Frechet_derivative_fixture(order):
    with mp.workdps(90):
        k0 = mp.mpf(512 * 1024)
        k = 2 * k0
        jet = [mp.mpf(i - 7) / 10 for i in range(15)]
        direction = [mp.mpf((-1) ** i) / 2 for i in range(15)]
        upper = mp.mpf(20000 * 1024) / (k * k0)
        for component in range(4):
            f = lambda e, component=component: literal_source(
                [x + e * d for x, d in zip(jet, direction)], k0, k
            )[component]
            actual = abs(mp.diff(f, mp.mpf(0), order))
            assert actual < mp.factorial(order) * upper


def test_nonzero_null_gradient_full_source_not_divided_away():
    with mp.workdps(90):
        jet = [mp.mpf(0)] * 15
        jet[1] = jet[2] = mp.mpf("0.5")
        jet[5] = mp.mpf(1)
        k0 = mp.mpf(512 * 1024)
        actual = literal_source(jet, k0, k0)
        expected = mp.mpf(1024) / (8 * k0 * k0)
        assert abs(actual[0] - expected) < expected * mp.mpf("1e-70")
        assert actual[0] == actual[1] and actual[2] == actual[3] == 0


@pytest.mark.parametrize("q", (0, 1, 1000))
def test_complete_longitudinal_covariance_and_temporal_constraint(q):
    d = noise.data()
    symbols = {str(x): x for x in d["complete_polarization_covariance"].free_symbols}
    P = d["complete_polarization_covariance"].subs(
        {symbols["m"]: 1000, symbols["q"]: q}
    )
    assert P == P.T
    assert set(P.eigenvals()).issubset(
        {s.S.Zero, s.S.One, 1 + s.Rational(2 * q * q, 1000**2)}
    )
    if q > 0:
        assert P[0, 0] > 0
    assert P.rank() == 3


@pytest.mark.parametrize("mass", (2, 1000))
@pytest.mark.parametrize("component", ("time", "spatial"))
def test_complete_spectral_covariance_below_time_space_norm(mass, component):
    # A rectangular time pulse is in the norm completion of smooth compact
    # tests; Gaussian spatial Fourier data remain Schwartz. No cutoff.
    with mp.workdps(45):
        m = mp.mpf(mass)
        radial = lambda q: q * q * mp.exp(-q * q) / (2 * mp.pi**2)

        def integrand(q):
            w = mp.sqrt(m * m + q * q)
            pol = q * q / (m * m) if component == "time" else 1 + q * q / (3 * m * m)
            return radial(q) * pol / (2 * w) * (2 * mp.sin(w / 2) / w) ** 2

        value = mp.quad(integrand, [0, 1, 4, mp.inf])
        spatial = mp.quad(radial, [0, 1, 4, mp.inf])
        gradient = mp.quad(lambda q: q * q * radial(q), [0, 1, 4, mp.inf])
        upper = spatial / (2 * m) + gradient / m**3
        assert 0 < value < upper


@pytest.mark.parametrize("momentum", (0.0, 10.0, 1000.0))
def test_forced_longitudinal_energy_and_actual_temporal_readout(momentum):
    mass = 1000.0
    q = momentum

    def source_values(t):
        x = t / 0.01
        b = np.exp(4 - 1 / (x * (1 - x))) if 0 < x < 1 else 0.0
        return b, 0.4 * b

    def rhs(t, y):
        J0, J = source_values(t)
        A, p = y[:2]
        F = np.sqrt((1 + q * q / mass**2) * J * J + q * q * J0 * J0 / mass**2)
        return [(1 + q * q / mass**2) * p - q * J0 / mass**2, -(mass**2) * A + J, F]

    sol = solve_ivp(
        rhs,
        (0, 0.01),
        [0.0, 0.0, 0.0],
        method="DOP853",
        rtol=1e-10,
        atol=1e-14,
        max_step=0.00005,
        dense_output=True,
    )
    assert sol.success
    for t in np.linspace(0, 0.01, 101):
        A, p, integral = sol.sol(t)
        J0, _ = source_values(t)
        E = ((1 + q * q / mass**2) * p * p + mass**2 * A * A) / 2
        A0 = (q * p - J0) / mass**2
        free = (p * p + mass**2 * A * A + mass**2 * A0 * A0) / 2
        assert E <= integral**2 / 2 + 1e-13
        assert free <= 2 * E + J0 * J0 / mass**2 + 1e-15
    assert max(abs(sol.sol(0.01)[:2])) > 0


@pytest.mark.parametrize("y", (-4, -1, 1, 4))
def test_full_fixed_vacuum_profile_value_and_derivative(y):
    with mp.workdps(120):
        k0 = mp.mpf("1e800")
        eps = mp.mpf("1e-770")
        phi = mp.mpf("0.25")

        # Bounded reference-profile fixtures, not substitute quantum states.
        def coefficient(ph, Y):
            u = ph / mp.sqrt(k0)
            x = Y / k0
            r = eps * (mp.mpf("0.25") + u)
            p = eps * (mp.mpf("0.5") - u)
            pv = eps * mp.mpf("1e-19")
            T = x**1024 / (x**1024 + (1 - x) ** 1024)
            return k0 * T * (-p + pv - (r + p) * (x - 1) / 2)

        actual = abs(coefficient(phi, mp.mpf(y)))
        derivative = abs(mp.diff(lambda Y: coefficient(phi, Y), mp.mpf(y)))
        x = abs(mp.mpf(y) / k0)
        value_upper = 8 * k0 * eps * x**1024
        derivative_upper = 8 * 1024 * eps * x**1023 + 2 * eps * x**1024
        assert 0 < actual < value_upper
        assert 0 < derivative < derivative_upper
        assert mp.log10(actual) < -818500
        assert mp.log10(derivative) < -818500


@pytest.mark.parametrize("ratio", (1, 2, 10, 1000))
def test_actual_family_decay_and_fixed_scalar_boundary(ratio):
    k, k0 = source.K, source.K0
    e = energy.data()
    n = noise.data()
    for key in (
        "reference_energy_coefficient",
        "full_scalar_causal_mean_force_coefficient",
    ):
        expr = e[key]
        assert s.factor(expr.subs(k, ratio * k0) * ratio - expr.subs(k, k0)) == 0
    expr = n["variance_coefficient"]
    assert s.factor(expr.subs(k, ratio * k0) * ratio - expr.subs(k, k0)) == 0
    assert "persists" in audit.observable()["fixed_profile"]
    assert "not a regulator-removed interacting" in vacuum.data()["quantum_scope"]


def test_exact_exponent_certificate_without_rounded_logarithms():
    d = vacuum.data()
    assert 4**1024 < 10**617
    assert d["anchor_action_decimal_exponent_upper"] == -818552
    assert d["anchor_force_decimal_exponent_upper"] == -818547
    assert d["anchor_force_decimal_exponent_upper"] < -818500
