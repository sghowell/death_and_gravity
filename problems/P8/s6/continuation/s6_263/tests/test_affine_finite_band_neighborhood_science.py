"""Independent whole-Hamiltonian, source, interval and comparison-flow tests."""

import math
from functools import cache

import numpy as np
import pytest
import sympy as s
from p8_affine import verify as certificate
from p8_vacuum_affine_coupled_principal_obstruction import coupled as c
from p8_vacuum_affine_finite_band_neighborhood import (
    audit,
    background,
    source,
    window,
)
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate
from scipy import integrate, linalg


@pytest.mark.parametrize("name", tuple(audit.packets()))
def test_all_complete_exact_packets(name):
    data = audit.packets()[name]
    assert certificate.certify_residuals(data["checks"])
    assert all(bool(value) for value in data["gates"].values())


@pytest.mark.parametrize(
    "name,call,args",
    audit.bad_cases(),
    ids=lambda value: value if isinstance(value, str) else None,
)
def test_every_unsupported_input_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_full_unchanged_frontier_and_corrected_physical_scope():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 119
    assert audit.matching()[:-1] == audit.previous.matching()
    assert all(audit.gates().values())
    assert "OPEN" in audit.observable()["original_problem"]
    assert "refuted" in audit.observable()["fixed_reference"]
    assert audit.require_band(source.LOW, source.HIGH) == (source.LOW, source.HIGH)
    assert audit.require_times(-source.TIME, source.TIME) == (-source.TIME, source.TIME)


@pytest.mark.parametrize("value", [-source.EPSILON, 0, source.EPSILON])
def test_closed_signed_initial_family(value):
    assert audit.require_epsilon(value) == value


@pytest.mark.parametrize("order", range(1, 6))
def test_complete_lapse_chain_bound_against_independent_Bell_polynomials(order):
    x = [s.Integer(value) for value in (3, 7, 25, 121, 1000)]
    independent = sum(s.bell(order, k, x[: order - k + 1]) for k in range(1, order + 1))
    assert independent == source.data()["whole_lapse_Bell_chain_factors"][order]


def comparison_values(P, seed):
    # Finite algebra fixtures, not a different physical state or a sample
    # proving the actual huge-P interval bound.
    values = {
        c.P: P,
        c.a: 1.1 + 0.03 * seed,
        c.D: 1.05,
        c.Z: 0.95,
        c.J: 1.52,
        c.Y: 1.04,
        c.Yv: 0.87,
        c.zeta: 0.13,
        c.C: 0.48,
        c.L2: -0.96,
        c.r: -0.013,
        c.rN: -2.02,
        c.th: 0.08,
        c.dH: 0.003,
        c.H: 0.04,
        c.cs[0]: 0.11,
        c.cs[1]: 0.006,
        c.ws[0]: 0.052,
        c.ws[1]: -0.004,
        c.ds[1]: 0.008,
        c.L0: 0.017,
        c.Vvv: -0.014,
        c.vs[0]: 0.006,
        c.vs[1]: -0.003,
    }
    named = {
        "fixed_heavy_mass2": 91 + 17 * seed,
        "positive_mass_ratio": 0.98,
        "K_over_zeta": 0.93,
        "heavy_energy_weight": math.sqrt(P * P + 91 + 17 * seed),
        "fixed_reference_Sbar00": 0.002,
        "fixed_reference_Sbar01": 0.003,
        "fixed_reference_Sbar11": -0.004,
    }
    packet = window.matrices()
    A = packet["whole_sixteen_phase_generator_without_common_reference_connection"]
    values.update(
        {
            symbol: named[symbol.name]
            for symbol in A.free_symbols
            if symbol not in values
        }
    )
    return values, named


def independent_coupled_hamiltonian(x, values, named):
    Qb, sm, sh, QL, Pb, pm, ph, PL = x
    get = values.__getitem__
    a, P, D, Z, J, zeta = [get(key) for key in (c.a, c.P, c.D, c.Z, c.J, c.zeta)]
    r, rN, theta, dH = [get(key) for key in (c.r, c.rN, c.th, c.dH)]
    charge = get(c.cs[0]) * sm + get(c.cs[1]) * sh
    linear = -2 * a * theta * P * P * Qb / D + 3 * a**3 * theta * charge / D
    linear -= (get(c.L2) / 2 + get(c.L0) * a * a / (2 * P * P)) * Pb
    linear -= (get(c.ws[0]) * pm + get(c.ws[1]) * ph) / Z + a**3 * get(c.ds[1]) * sh
    linear -= 3 * math.sqrt(zeta) * (r * theta / D + rN * dH) * PL
    result = linear**2 / (4 * J * a**3) + (pm * pm + ph * ph) / (2 * Z * a**3)
    result += (
        charge * a * P * P * Qb / D
        + 3 * charge * r * math.sqrt(zeta) * PL / (2 * D)
        - 3 * a**3 * charge**2 / (4 * D)
    )
    result += a * get(c.Y) * P * P * (sm * sm + sh * sh) / 2
    result -= get(c.C) * Pb * Pb / (2 * a * P * P) + a * get(c.Vvv) * Pb * Pb / (
        4 * P**4
    )
    result -= a * a * (get(c.vs[0]) * sm + get(c.vs[1]) * sh) * Pb / (2 * P * P)
    result += (
        a**3 * named["positive_mass_ratio"] * named["fixed_heavy_mass2"] * sh * sh / 2
    )
    result -= r * math.sqrt(zeta) * P * P * Qb * PL / (D * a * a)
    result += zeta * (1 - 3 * Z * r * r / (2 * D)) * PL * PL / (2 * Z * a**3)
    result += PL * PL / (2 * named["K_over_zeta"] * P * P * a**3)
    result += a**3 * get(c.Yv) * P * P * QL * QL / (2 * zeta) - get(c.H) * Qb * Pb
    return result


@pytest.mark.parametrize("P", [7.0, 19.0, 41.0])
@pytest.mark.parametrize("seed", range(3))
def test_all_sixteen_generator_entries_against_independent_quadratic_polarization(
    P, seed
):
    values, named = comparison_values(P, seed)
    basis = np.eye(8)
    single = np.array(
        [independent_coupled_hamiltonian(e, values, named) for e in basis]
    )
    H = np.empty((8, 8))
    for i in range(8):
        H[i, i] = 2 * single[i]
        for j in range(i):
            H[i, j] = H[j, i] = (
                independent_coupled_hamiltonian(basis[i] + basis[j], values, named)
                - single[i]
                - single[j]
            )
    zeta = values[c.zeta]
    R = np.eye(8)
    R[4:6, :2] = P * P * np.array([[0.002, 0.003], [0.003, -0.004]])
    W = np.diag(
        [
            P,
            P,
            named["heavy_energy_weight"],
            P / math.sqrt(zeta),
            1,
            1,
            1,
            math.sqrt(zeta),
        ]
    )
    omega = np.array(c.OMEGA, float)
    coupled = (
        W
        @ linalg.solve(R, omega @ H @ R, check_finite=True)
        @ linalg.inv(W, check_finite=True)
    )
    a, D, Z, Y, Yv = [values[key] for key in (c.a, c.D, c.Z, c.Y, c.Yv)]
    kt = named["K_over_zeta"]
    TT = np.array([[0, P / (a**3 * D)], [-2 * a * values[c.C] * P, 0]])
    PT = np.array(
        [[0, P / (a**3 * kt)], [-a * kt * Y * P / Z - a**3 * Yv / (zeta * P), 0]]
    )
    wanted = linalg.block_diag(coupled, TT, TT, PT, PT)
    actual = np.array(
        window.matrices()[
            "whole_sixteen_phase_generator_without_common_reference_connection"
        ].subs(values),
        float,
    )
    np.testing.assert_allclose(actual, wanted, rtol=4e-9, atol=3e-7)
    assert np.linalg.norm(actual[:8, :8] - np.diag(np.diag(actual[:8, :8]))) > 1
    # Every heavy contact is nonzero in the independent fixture.
    assert all(values[key] != 0 for key in (c.cs[1], c.ws[1], c.ds[1]))


@pytest.mark.parametrize("parameter", [c.D, c.J, c.r, c.th, c.L0, c.ws[1]])
def test_complete_generator_derivatives_against_independent_differences(parameter):
    values, _ = comparison_values(13.0, 1)
    A = window.matrices()[
        "whole_sixteen_phase_generator_without_common_reference_connection"
    ]
    step = 1e-5
    plus, minus = dict(values), dict(values)
    plus[parameter] += step
    minus[parameter] -= step
    difference = (np.array(A.subs(plus), float) - np.array(A.subs(minus), float)) / (
        2 * step
    )
    expected = np.array(A.diff(parameter).subs(values), float)
    np.testing.assert_allclose(difference, expected, rtol=2e-7, atol=2e-7)


def test_reference_time_connection_is_retained_but_cancels_in_difference():
    data = window.matrices()
    connection = data["whole_common_fixed_reference_time_connection"]
    assert connection != s.zeros(16)
    assert all(
        not connection.has(variable) for variable in (c.D, c.J, c.r, c.th, c.a, c.H)
    )
    assert (
        len(data["whole_all_24_coefficient_generator_derivative_sum_enclosures"]) == 24
    )
    assert data["whole_generator_coefficient_Lipschitz_sum"] < 10**180


@pytest.mark.parametrize("sign", [-1, 1])
def test_centered_exact_root_resolves_actual_tiny_initial_lapse_radius(sign):
    epsilon = sign * source.EPSILON
    lapse = 1 + epsilon
    m2 = (1120 * lapse**4 - 2847 * lapse**2 + 1728) / (100 * lapse**2)
    difference = m2 - s.Rational(1, 100)
    assert difference != 0
    assert s.Rational(9, 100) ** 2 < m2 < s.Rational(11, 100) ** 2
    # Work with the rational numerator before the square root. Binary
    # floating-point subtraction would silently identify both lapses with1.
    ratio = difference / epsilon
    assert -13 < ratio < -12
    centered_bound = abs(difference) / s.Rational(19, 100)
    assert centered_bound < s.Rational(1, 10**225)
    assert float(lapse) == 1.0 and float(m2) == 0.01


def test_pointwise_zero_Ru_does_not_erase_the_initial_Ruu_constraint_contact():
    lapse = s.Rational(9999, 10000)
    n = s.Symbol("independent_initial_N", positive=True)
    R = n**-2
    Ruu = -6 * (n**-2 - 1)
    F = -(624 * n**-4 + 753 * n**-2 + 224) / 200
    U = R ** (-s.Rational(3, 4))
    Z = U / n
    contact = 3 * U * Ruu * s.diff(R, n) / (4 * R * n)
    G = s.diff(n * U * F, n) - contact
    m2 = (1120 * n**4 - 2847 * n**2 + 1728) / (100 * n**2)
    actual = s.factor(G + s.diff(Z, n) * m2 / 2)
    erased = s.factor(s.diff(n * U * F, n) + s.diff(Z, n) * m2 / 2)
    assert actual == 0
    assert Ruu.subs(n, lapse) != 0
    assert erased.subs(n, lapse) != 0
    assert s.factor(erased - contact) == 0


@pytest.mark.parametrize("order", range(1, 6))
def test_independent_point_samples_inside_complete_lapse_derivative_intervals(order):
    n = background.N
    expression = s.diff(n**-2, n, order)
    interval = evaluate(
        expression, {n: I(s.Rational(999, 1000), s.Rational(1001, 1000))}
    )
    for numerator in range(9990, 10011):
        point = s.Rational(numerator, 10000)
        value = expression.subs(n, point)
        assert interval.lo <= value <= interval.hi
    # Independent falling-power derivative formula.
    assert (
        s.factor(
            expression - (-1) ** order * s.factorial(order + 1) * n ** (-2 - order)
        )
        == 0
    )


@pytest.mark.parametrize("P", [1e64, 1.5e64, 2e64])
@pytest.mark.parametrize("a", [1.0, 1.3, 2.0])
def test_positive_full_reference_energy_with_noncommuting_scalar_matrices_and_actual_mass(
    P, a
):
    K = np.array([[1.8, 0.31], [0.31, 0.7]])
    G = np.array([[0.9, 0.25], [0.25, 1.6]])
    gyro = np.array([[0.0, 7e17], [-7e17, 0.0]])
    assert np.linalg.norm(K @ G - G @ K) > 0.1
    Kinv = linalg.inv(K, check_finite=True)
    velocity = np.column_stack((-Kinv @ gyro / P, Kinv / a**3))
    Qlight = velocity.T @ K @ velocity
    Qlight[:2, :2] += G / a**2
    mass = float(background.b.source.MASS2)
    zeta = 1e-6
    # Whole normalized oscillators: heavy, longitudinal Proca, two TT,
    # two transverse Proca pairs. No mass expansion or commuting K/G.
    weights = [
        [(a * P * P + a**3 * mass) / (P * P + mass), a**-3],
        [a, a**-3 + 1 / (a * zeta * P * P)],
        [a, a**-3],
        [a, a**-3],
        [a**-1 + a / (zeta * P * P), a**-1],
        [a**-1 + a / (zeta * P * P), a**-1],
    ]
    Q = linalg.block_diag(Qlight, *[np.diag(pair) for pair in weights])
    spectrum = linalg.eigvalsh(Q, check_finite=True)
    assert Q.shape == (16, 16) and spectrum[0] > 1e-12 and spectrum[-1] < 1e12
    x = np.linspace(-1.2, 1.1, 16)
    y = x[:2] / P
    Pi = x[2:4]
    ydot = Kinv @ (Pi / a**3 - gyro @ y)
    direct = (ydot @ K @ ydot + (P / a) ** 2 * (y @ G @ y)) / 2
    direct += sum(
        np.dot(pair, x[4 + 2 * i : 6 + 2 * i] ** 2) / 2
        for i, pair in enumerate(weights)
    )
    assert float(x @ Q @ x / 2) == pytest.approx(direct, rel=2e-14)
    assert mass > 1e197 and mass > P * P


@pytest.mark.parametrize("seed", range(3))
@pytest.mark.parametrize("direction", [-1, 1])
def test_independent_sixteen_phase_nonstationary_energy_flow(seed, direction):
    # A diagnostic comparison theorem fixture, NOT an integration of the
    # actual giant-mass parent or a replacement for exact interval proofs.
    rng = np.random.default_rng(3100 + seed)
    orthogonal, _ = linalg.qr(rng.normal(size=(16, 16)), check_finite=True)
    rates = np.linspace(-0.3, 0.3, 16)
    core = linalg.block_diag(
        *[np.array([[0.0, w], [-w, 0.0]]) for w in range(9, 81, 10)]
    )
    perturbation = 0.02 * np.eye(16) + 0.005 * np.diag(np.ones(15), 1)

    def moving(t):
        B = np.exp(rates * t)[:, None] * orthogonal
        return B, rates[:, None] * B

    def rhs(t, x):
        B, Bdot = moving(t)
        return (
            linalg.solve(B, core @ (B @ x) - Bdot @ x, check_finite=True)
            + perturbation @ x
        )

    start, stop = -direction * 0.1, direction * 0.1
    initial = rng.normal(size=16)
    solution = integrate.solve_ivp(
        rhs,
        (start, stop),
        initial,
        method="DOP853",
        rtol=2e-12,
        atol=2e-13,
        dense_output=True,
    )
    assert solution.success
    B0, _ = moving(start)
    root0 = linalg.norm(B0 @ initial)
    normdelta = linalg.norm(perturbation, 2)
    rootrate = 1.2 * normdelta
    for t in np.linspace(start, stop, 41):
        B, Bdot = moving(t)
        x = solution.sol(t)
        root = linalg.norm(B @ x)
        assert root <= root0 * math.exp(rootrate * abs(t - start)) * (1 + 2e-10)
        # Connection is really nonzero. Omitting it loses exact cancellation.
        assert linalg.norm(Bdot @ x) > 1e-3
        exact_derivative = (B @ x) @ (Bdot @ x + B @ rhs(t, x))
        wanted = (B @ x) @ (B @ perturbation @ x)
        assert exact_derivative == pytest.approx(wanted, rel=2e-10, abs=2e-10)


@pytest.mark.parametrize("P", [11.0, 31.0])
def test_complete_original_canonical_endpoint_maps_are_mutual_inverses(P):
    values, _ = comparison_values(P, 1)
    data = window.matrices()
    forward = np.array(
        data["whole_original_canonical_phase_to_energy_coordinates"].subs(values), float
    )
    inverse = np.array(
        data["whole_energy_coordinates_to_original_canonical_phase"].subs(values), float
    )
    np.testing.assert_allclose(forward @ inverse, np.eye(16), rtol=1e-10, atol=2e-9)
    np.testing.assert_allclose(inverse @ forward, np.eye(16), rtol=1e-10, atol=2e-9)


@cache
def bare_nonlinear_comparison_functions():
    # Independent finite-time algebra/ODE fixture only. Actual sources
    # remain present in every certified estimate and are not evaluated by
    # this bare-function diagnostic.
    b = background.b
    substitutions = {}
    for function, tree in [
        (b.R, background.Rtree),
        (b.F, background.Ftree),
        (b.j, s.Integer(0)),
    ]:
        for i in range(6):
            for j in range(6 - i):
                substitutions[s.diff(function, b.u, i, b.N, j)] = s.diff(
                    tree, b.u, i, b.N, j
                )
    substitutions.update({b.h: 0, b.mh: 0})
    J = s.factor(background.Jactual.xreplace(substitutions))
    expressions = [b.expressions[name] for name in ("Ndot", "Hdot", "M1_acceleration")]
    expressions.extend(
        [b.a * b.H, b.constraint, background.physical()["whole_actual_physical_Hubble"]]
    )
    expressions = [
        s.factor(expr.xreplace(substitutions).subs(b.J, J)) for expr in expressions
    ]
    assert not any(expr.has(s.Derivative) for expr in expressions)
    return s.lambdify((b.u, b.N, b.H, b.mc, b.a), expressions, "numpy", cse=True)


@pytest.mark.parametrize("epsilon", [-1e-6, 0.0, 1e-6])
@pytest.mark.parametrize("direction", [-1, 1])
def test_bare_comparison_nonlinear_flow_preserves_constraint_and_physical_turn(
    epsilon, direction
):
    # These finite-size numerical fixtures are outside the certified
    # signed epsilon family; they do not refute S256's full-system growth.
    rates = bare_nonlinear_comparison_functions()
    lapse = 1 + epsilon
    m2 = (1120 * lapse**4 - 2847 * lapse**2 + 1728) / (100 * lapse**2)
    initial = np.array([lapse, 0.0, math.sqrt(m2), 1.0])
    solution = integrate.solve_ivp(
        lambda t, x: np.asarray(rates(t, *x)[:4], float),
        (0.0, direction * 1e-4),
        initial,
        method="DOP853",
        rtol=2e-12,
        atol=2e-13,
        dense_output=True,
    )
    assert solution.success
    for t in np.linspace(0, direction * 1e-4, 17):
        state = solution.sol(t)
        values = rates(t, *state)
        assert state[0] > 0.99 and state[3] > 0.99
        assert abs(values[4]) < 2e-10
        if t != 0:
            assert direction * values[5] > 0
        if epsilon == 0:
            wanted = np.array(
                [1, 4 * t / (1 + t * t), 1 / (10 * (1 + t * t) ** 6), (1 + t * t) ** 2]
            )
            np.testing.assert_allclose(state, wanted, rtol=2e-10, atol=2e-11)
    if epsilon:
        with pytest.raises(ValueError):
            audit.require_epsilon(s.Rational(str(epsilon)))
