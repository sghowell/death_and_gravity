"""Independent complete-source, constrained-energy and causal-force controls."""

from functools import cache
from itertools import product

import mpmath as mp
import numpy as np
import pytest
import sympy as s
from p8_affine_vacuum_domain import family
from p8_vacuum_affine_retarded_energy import (
    audit,
    clock,
    coefficients,
    estimates,
    hamiltonian,
)
from scipy.integrate import solve_ivp

PACKETS = {
    "hamiltonian": hamiltonian.data,
    "modes": hamiltonian.modes,
    "coefficients": coefficients.data,
    "clock": clock.data,
    "clock_bounds": clock.bounds,
    "estimates": estimates.data,
}


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_exact_residual(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_unsupported_scope_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("packet", list(PACKETS))
def test_all_quantitative_gates(packet):
    assert all(PACKETS[packet]().get("gates", {}).values())


@pytest.mark.parametrize(
    "x",
    [
        s.Rational(9, 10),
        s.Rational(19, 20),
        s.S.One,
        s.Rational(21, 20),
        s.Rational(11, 10),
    ],
)
def test_literal_full_switch_derivatives_on_both_sides(x):
    X = family.X
    n = family.N
    d = family.data()
    q = (1 - x) / x
    t1 = n * q ** (n - 1) / (x * x * (1 + q**n) ** 2)
    t2 = (
        -n * (n - 1) * q ** (n - 2) / (x**4 * (1 + q**n) ** 2)
        - 2 * n * q ** (n - 1) / (x**3 * (1 + q**n) ** 2)
        + 2 * n * n * q ** (2 * n - 2) / (x**4 * (1 + q**n) ** 3)
    )
    assert s.cancel(s.diff(d["T"], X).subs(X, x) - t1) == 0
    assert s.cancel(s.diff(d["T"], X, 2).subs(X, x) - t2) == 0
    assert abs(t1) < 1 and abs(t2) < 1
    w1 = s.diff(d["bump"], X).subs(X, x)
    w2 = s.diff(d["bump"], X, 2).subs(X, x)
    assert abs(s.N(w1, 90)) < 1
    assert abs(s.N(w2, 90)) < 1
    bounds = coefficients.data()
    assert abs(t1) <= bounds["rational_step_derivative_bounds"][0]
    assert abs(t2) <= bounds["rational_step_derivative_bounds"][1]


ALPHAS = [a for a in product(range(4), repeat=4) if sum(a) <= 3]


def index(*directions):
    a = [0] * 4
    for d in directions:
        a[d] += 1
    return tuple(a)


@cache
def exact_R_functions():
    d = family.data()
    return s.lambdify(
        (family.u, family.X), (d["R"], d["RX"], s.diff(d["R"], family.u)), "mpmath"
    )


def literal_source(t, jets):
    get = lambda *ds: jets.get(index(*ds), mp.mpf("0"))
    a = (1 + t * t) ** 2
    H = 4 * t / (1 + t * t)
    u = t + get()
    p = mp.matrix([1 + get(0), get(1), get(2), get(3)])
    inv = mp.diag([1, -1 / (a * a), -1 / (a * a), -1 / (a * a)])
    Hess = mp.matrix(4)
    Hess[0, 0] = get(0, 0)
    for i in range(1, 4):
        Hess[0, i] = Hess[i, 0] = get(0, i) - H * get(i)
        for j in range(1, 4):
            Hess[i, j] = get(i, j) - (a * a * H * (1 + get(0)) if i == j else 0)
    X = (p.T * inv * p)[0]
    Q = sum((inv * Hess)[i, i] for i in range(4))
    Z = (p.T * inv * Hess * inv * p)[0]
    R, RX, Ru = exact_R_functions()(u, X)
    Href = 4 * u / (1 + u * u)
    scalar = (R - 1) * (
        -3 * Href + 3 * Ru / (4 * R) + Q / X + (-1 / (X * X) + 3 * RX / (2 * R * X)) * Z
    )
    return p * scalar


def directional_source(t, jets, direction):
    def shifted(h):
        return {a: jets[a] + h * direction.get(a, 0) for a in ALPHAS}

    return mp.matrix(
        [mp.diff(lambda h, i=i: literal_source(t, shifted(h))[i], 0) for i in range(4)]
    )


@pytest.mark.parametrize("t", ["-0.4", "0", "0.3"])
@pytest.mark.parametrize("seed", [1, 2, 3])
def test_actual_full_source_spatial_and_arbitrary_Frechet_tame_bounds(t, seed):
    with mp.workdps(80):
        t = mp.mpf(t)
        delta = mp.mpf("0.01")
        jets = {
            a: delta * mp.mpf(((i * seed + 2) % 5) - 2) / 2
            for i, a in enumerate(ALPHAS)
        }
        eta = {
            a: mp.mpf(((i * (seed + 1) + 1) % 7) - 3) / 3 for i, a in enumerate(ALPHAS)
        }
        G = mp.sqrt(sum(v * v for v in jets.values()))
        Geta = mp.sqrt(sum(v * v for v in eta.values()))
        S = literal_source(t, jets)
        assert mp.norm(S) <= 512 * delta * G
        DS = directional_source(t, jets, eta)
        assert mp.norm(DS) <= 2048 * delta * Geta
        for i in (1, 2, 3):
            direction = {}
            for alpha in ALPHAS:
                beta = list(alpha)
                beta[i] += 1
                direction[alpha] = jets.get(tuple(beta), 0)
            dS = directional_source(t, jets, direction)
            assert mp.norm(dS) <= 2048 * delta * G


@pytest.mark.parametrize("x", ["0", "0.2", "-0.3"])
def test_actual_nonclosed_quadratic_source_coefficient(x):
    with mp.workdps(85):
        x = mp.mpf(x)
        c = 1 + x

        def profile(eps):
            js = {a: mp.mpf("0") for a in ALPHAS}
            js[index(0)] = eps * c
            js[index(0, 0)] = eps * c
            js[index(0, 1)] = eps
            js[index(0, 0, 1)] = eps
            return literal_source(mp.mpf("0"), js)

        actual = mp.matrix(
            [mp.diff(lambda e, i=i: profile(e)[i], 0, 2) / 2 for i in range(4)]
        )
        assert abs(actual[0] - 3 * c * c) < mp.mpf("1e-65")
        assert mp.norm(mp.matrix([actual[i] for i in (1, 2, 3)])) < mp.mpf("1e-65")


@pytest.mark.parametrize("kind", ["T", "L"])
@pytest.mark.parametrize("momentum", [0.0, 2.0, 1000.0])
def test_forced_actual_CD_modes_obey_finite_time_energy_bound(kind, momentum):
    m = 1000.0
    k = momentum
    B = float(s.Rational(5, 4) ** 6)

    def sources(t):
        q = 1 - 4 * t * t
        bump = np.exp(-1 / q) if q > 0 else 0.0
        return 2 * bump, (1 + t) * bump

    def rhs(t, y):
        A, pi, _ = y
        a = (1 + t * t) ** 2
        j0, j = sources(t)
        if kind == "L":
            dA = pi / a + k * k * pi / (a**3 * m * m) - k * j0 / (m * m)
            dpi = -a * m * m * A + a * j
            force = np.sqrt(
                a * j * j + k * k * j * j / (a * m * m) + a * k * k * j0 * j0 / (m * m)
            )
        else:
            dA = pi / a
            dpi = -(a * m * m + k * k / a) * A + a * j
            force = np.sqrt(a) * abs(j)
        return [dA, dpi, force]

    times = np.linspace(-0.5, 0.5, 81)
    sol = solve_ivp(
        rhs,
        (-0.5, 0.5),
        [0.0, 0.0, 0.0],
        method="DOP853",
        t_eval=times,
        rtol=1e-10,
        atol=[1e-14, 1e-12, 1e-12],
        max_step=0.002,
    )
    assert sol.success and sol.status == 0 and sol.t[-1] == 0.5
    for t, A, pi, integ in zip(sol.t, *sol.y):
        a = (1 + t * t) ** 2
        E = pi * pi / (2 * a) + a * m * m * A * A / 2
        if kind == "L":
            E += k * k * pi * pi / (2 * a**3 * m * m)
        else:
            E += k * k * A * A / (2 * a)
        assert np.sqrt(2 * E) <= B * integ + 1e-9
    assert sol.y[2, -1] > 0


@pytest.mark.parametrize(
    "a,m,k,j0,j",
    [
        (1, 1000, 0, 2, 3),
        (s.Rational(25, 16), 1000, 1, 2, -3),
        (s.Rational(121, 100), 1000, 1000, -2, 1),
    ],
)
def test_literal_longitudinal_temporal_constraint_and_contact(a, m, k, j0, j):
    A0, A, Ad, pi = s.symbols("A0 A A_dot pi", real=True)
    # Spatial A,pi,J in the normalized sin mode; J0,A0 in cos mode.
    L = (
        a * (Ad + k * A0) ** 2 / 2
        + a**3 * m * m * A0 * A0 / 2
        - a * m * m * A * A / 2
        - a**3 * A0 * j0
        + a * A * j
    )
    momentum = s.diff(L, Ad)
    velocity = s.solve(momentum - pi, Ad)[0]
    Ham = s.expand((pi * Ad - L).subs(Ad, velocity))
    solved = s.solve(s.diff(Ham, A0), A0)[0]
    expected = (j0 - k * pi / a**3) / (m * m)
    assert s.factor(solved - expected) == 0
    reduced = s.factor(Ham.subs(A0, solved))
    E = pi * pi / (2 * a) + a * m * m * A * A / 2 + k * k * pi * pi / (2 * a**3 * m * m)
    assert (
        s.factor(
            reduced
            - E
            + j0 * k * pi / (m * m)
            - a**3 * j0 * j0 / (2 * m * m)
            + a * A * j
        )
        == 0
    )
    assert (
        s.factor(
            s.diff(reduced, pi)
            - (pi / a + k * k * pi / (a**3 * m * m) - k * j0 / (m * m))
        )
        == 0
    )


@pytest.mark.parametrize("S0,variation", [(1, 2), (-2, 3), (s.Rational(2, 3), -1)])
def test_homogeneous_source_full_force_cancels_only_with_contact(S0, variation):
    K = s.Integer(10) ** 800
    z = s.Rational(1, 10**6)
    J0 = s.sqrt(K / z) * S0
    A0 = z * J0
    W0 = A0 / s.sqrt(K * z)
    full = K * (S0 - W0) * variation
    assert full == 0
    assert -K * W0 * variation != 0
    assert K * S0 * variation != 0


def test_causal_force_is_not_single_branch_retarded_action_variation():
    x, y = s.symbols("source_early source_late", real=True)
    source = s.Matrix([x, y])
    ret = s.Matrix([[1, 0], [2, 1]])
    causal = (s.eye(2) - ret) * source
    action = (source.T * (s.eye(2) - ret) * source)[0] / 2
    gradient = s.Matrix([s.diff(action, v) for v in (x, y)])
    assert s.diff(causal[0], y) == 0
    assert s.diff(gradient[0], y) != 0
    assert gradient != causal


@pytest.mark.parametrize("omega,time", [(1, 1), (1000, s.Rational(1, 100)), (1000, 1)])
def test_finite_time_resonance_energy_does_not_diverge(omega, time):
    t, w = s.symbols("t w", positive=True)
    q = (s.sin(w * t) - w * t * s.cos(w * t)) / (2 * w * w)
    p = s.diff(q, t)
    E = (p * p + w * w * q * q).subs({w: omega, t: time}) / 2
    assert s.N(s.sqrt(2 * E), 80) <= s.N(time, 80)
    assert s.simplify(q.subs(t, 0)) == 0


def test_energy_and_causal_force_constants_do_not_use_support_volume():
    d = estimates.data()
    assert d["normalized_reference_energy_coefficient"] < 3 * 10**12
    assert d["normalized_full_causal_light_force_constant"] < 4 * 10**10
    assert d["normalized_free_mean_energy_instantaneous_coefficient"] == 10**6
    assert "pointwise" in d["energy_result"]
    assert "volume" in clock.bounds()["norm_definition"]
