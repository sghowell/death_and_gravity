"""Independent full-function jets, curved mode physics and complete integrals."""

import json
from math import comb
from pathlib import Path

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_curved_dirac_state import (
    audit,
    geometry,
    hadamard,
    transition,
    tube,
)


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_named_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_unsupported_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def mul(a, b, n):
    return [sum(a[k] * b[j - k] for k in range(j + 1)) for j in range(n + 1)]


def power(a, alpha, n):
    # Solve a*b'=alpha*a'*b coefficient by coefficient for the full local germ.
    out = [mp.power(a[0], alpha)]
    for j in range(1, n + 1):
        out.append(
            sum(((alpha + 1) * k - j) * a[k] * out[j - k] for k in range(1, j + 1))
            / (j * a[0])
        )
    return out


def derivative(a):
    return [(j + 1) * a[j + 1] for j in range(len(a) - 1)]


def exact_frames(t, p, m, d, tau, count):
    n = count + 1
    den = [
        mp.mpf(comb(8, j)) * (t / tau) ** (8 - j) / tau**j if j <= 8 else mp.mpf(0)
        for j in range(n + 1)
    ]
    den[0] += 1
    line = [t / tau, 1 / tau] + [mp.mpf(0)] * (n - 1)
    profile = mul(line, power(den, -mp.mpf(1) / 8, n), n)
    mass = [d * v for v in profile]
    mass[0] += m
    geom = [1 + t * t, 2 * t, mp.mpf(1)] + [mp.mpf(0)] * (n - 2)
    momentum = [p * v for v in power(geom, -2, n)]
    omega2 = [a + b for a, b in zip(mul(mass, mass, n), mul(momentum, momentum, n))]
    omega = power(omega2, mp.mpf(1) / 2, count)
    first = mul(momentum, derivative(mass), count)
    second = mul(mass, derivative(momentum), count)
    coupling = [
        v / 2
        for v in mul(
            [a - b for a, b in zip(first, second)], power(omega2, -1, count), count
        )
    ]
    result = [(omega[0], coupling[0])]
    e, g = omega, coupling
    for j in range(count):
        k = count - j
        e2 = mul(e, e, k)
        g2 = mul(g, g, k)
        total = [a + b for a, b in zip(e2, g2)]
        numerator = [
            a - b
            for a, b in zip(mul(e, derivative(g), k - 1), mul(g, derivative(e), k - 1))
        ]
        next_g = [
            (-1) ** j * v / 2 for v in mul(numerator, power(total, -1, k - 1), k - 1)
        ]
        e, g = power(total, mp.mpf(1) / 2, k - 1), next_g
        result.append((e[0], g[0]))
    return result


def direct_pair(t, p, m, d, tau, depth):
    if depth == 0:
        mass = lambda z: m + d * (z / tau) / mp.power(1 + (z / tau) ** 8, mp.mpf(1) / 8)
        q = lambda z: p / (1 + z * z) ** 2
        w = mp.sqrt(q(t) ** 2 + mass(t) ** 2)
        g = (q(t) * mp.diff(mass, t) - mass(t) * mp.diff(q, t)) / (2 * w * w)
        return w, g
    e, g = direct_pair(t, p, m, d, tau, depth - 1)
    ratio = lambda z: (
        direct_pair(z, p, m, d, tau, depth - 1)[1]
        / direct_pair(z, p, m, d, tau, depth - 1)[0]
    )
    return mp.sqrt(e * e + g * g), (-1) ** (depth - 1) * mp.diff(ratio, t) / (
        2 * (1 + (g / e) ** 2)
    )


@pytest.mark.parametrize("t", ["-2", "0.13", "3"])
@pytest.mark.parametrize("d", ["0", "0.01", "-0.01"])
def test_full_function_jets_against_independent_nested_derivatives(t, d):
    with mp.workdps(80):
        args = (mp.mpf(t), mp.mpf("1.4"), mp.mpf(2), mp.mpf(d), mp.mpf("0.7"))
        jets = exact_frames(*args, 4)
        for n in range(3):
            direct = direct_pair(*args, n)
            for a, b in zip(jets[n], direct):
                assert abs(a - b) < mp.mpf("1e-65") * max(1, abs(b))


@pytest.mark.parametrize("t", ["-100", "-2", "-0.2", "0", "0.2", "2", "100"])
def test_initial_variable_complex_disk_bounds(t):
    with mp.workdps(75):
        t = mp.mpf(t)
        tau = mp.mpf("0.2")
        m = mp.mpf("1e8")
        d = mp.mpf("3e5")
        p = mp.mpf("7e8")
        L = mp.sqrt(t * t + tau * tau)
        D = 1 + t * t
        q = p / D**2
        m0 = mp.mpf(".99") * m
        E = mp.sqrt(q * q + m0 * m0)
        w = mp.power(1 + (t / tau) ** 8, -mp.mpf(9) / 8)
        G = 16 * q / E**2 * (d / tau * w + m * L / D)
        for j in range(24):
            z = t + L / 1024 * mp.exp(2j * mp.pi * j / 24)
            profile = (z / tau) / mp.power(1 + (z / tau) ** 8, mp.mpf(1) / 8)
            assert abs(profile) < 2
            assert abs(mp.power(1 + (z / tau) ** 8, -mp.mpf(9) / 8)) < 2 * w
            qz = p / (1 + z * z) ** 2
            M = m + d * profile
            wz = mp.sqrt(qz * qz + M * M)
            Md = d / tau * mp.power(1 + (z / tau) ** 8, -mp.mpf(9) / 8)
            g = qz * (Md + 4 * z / (1 + z * z) * M) / (2 * wz * wz)
            assert abs(wz) > mp.mpf(".9") * E
            assert abs(g) <= G
            assert abs(qz * qz - q * q) < q * q / 16


@pytest.mark.parametrize("count", [1, 4, 20, 100, 1000])
def test_arbitrary_order_gap_and_floor(count):
    assert tube.radius(count, count) == s.Rational(1, 2048)
    assert tube.floor(count, count) > s.Rational(1, 2)
    assert s.Rational(1, 256 * count) ** 2 < s.Rational(1, 4 * count)


@pytest.mark.parametrize("t", ["-3", "0.3", "4"])
@pytest.mark.parametrize("amplitude", ["0", "1000"])
def test_twenty_full_frames_inside_nested_disks(t, amplitude):
    with mp.workdps(500):
        t = mp.mpf(t)
        p = mp.mpf("7e8")
        m = mp.mpf("1e8")
        d = mp.mpf(amplitude)
        tau = mp.mpf(".2")
        L = mp.sqrt(t * t + tau * tau)
        q = p / (1 + t * t) ** 2
        m0 = mp.mpf(".99") * m
        E = mp.sqrt(q * q + m0 * m0)
        w = mp.power(1 + (t / tau) ** 8, -mp.mpf(9) / 8)
        G0 = 16 * q / E**2 * (d / tau * w + m * L / (1 + t * t))
        ratio = mp.mpf(163840) / (L * E)
        for z in (t, t + 1j * L / 2048):
            frames = exact_frames(z, p, m, d, tau, 20)
            for j, (e, g) in enumerate(frames):
                assert abs(e) > mp.mpf(".9") * (1 - mp.mpf(1) / 80) ** j * E
                assert abs(g) <= G0 * ratio**j
            assert frames[-1][0] != frames[-2][0]


def test_actual_scale_twenty_frames_precision_refinement():
    answers = []
    for digits in (4800, 5200):
        with mp.workdps(digits):
            m = mp.mpf("1e200")
            p = m
            d = mp.mpf("3e197")
            tau = mp.mpf("1e-100")
            t = tau * mp.mpf(".7")
            values = exact_frames(t, p, m, d, tau, 20)
            L = mp.sqrt(t * t + tau * tau)
            q = p / (1 + t * t) ** 2
            E = mp.sqrt(q * q + (mp.mpf(".99") * m) ** 2)
            w = mp.power(1 + (t / tau) ** 8, -mp.mpf(9) / 8)
            G0 = 16 * q / E**2 * (d / tau * w + m * L / (1 + t * t))
            ratio = mp.mpf(163840) / (L * E)
            for j, (e, g) in enumerate(values):
                assert abs(g) <= G0 * ratio**j
                assert e >= E
            assert values[-1][1] != 0
            assert values[-1][0] - values[-2][0] > 0
            answers.append((values[-1][1], values[-1][0] - values[-2][0]))
    with mp.workdps(150):
        for a, b in zip(*answers):
            assert abs(a / b - 1) < mp.mpf("1e-100")


@pytest.mark.parametrize("count", [1, 4, 20, 80])
def test_high_momentum_condition_for_every_time_not_fixed_mass_gap(count):
    with mp.workdps(70):
        tau = mp.mpf(".1")
        m0 = mp.mpf(1)
        K = mp.mpf(tube.factor(count))
        p = max(8 * K / tau, 4 * (2 * K) ** 4 / m0**3)
        scale = (p / m0) ** (mp.mpf(1) / 4)
        for t in (0, tau, 1, scale / 2, scale, scale * 2, scale * 100):
            Z = mp.sqrt(t * t + tau * tau) * mp.sqrt(m0 * m0 + p * p / (1 + t * t) ** 4)
            assert Z >= 2 * K
        assert tau * m0 < 2 * K


@pytest.mark.parametrize("count", [1, 4, 20])
@pytest.mark.parametrize("p", ["0.1", "2", "100"])
def test_full_tail_substitution_and_integral_bounds(count, p):
    with mp.workdps(70):
        p = mp.mpf(p)
        m0 = mp.mpf("1.7")
        n = mp.mpf(count)
        direct = mp.quad(
            lambda t: t ** (-n - 5) / (m0 * m0 + p * p / (16 * t**8)) ** ((n + 2) / 2),
            [1, 2, mp.inf],
        )
        yform = (
            (4 / p) ** (n / 4)
            / p
            * mp.quad(
                lambda y: y ** (n / 4) / (y * y + m0 * m0) ** ((n + 2) / 2), [0, p / 4]
            )
        )
        assert abs(direct / yform - 1) < mp.mpf("1e-55")
        assert direct <= 1 / ((n + 4) * m0 ** (n + 2))


@pytest.mark.parametrize("power_value", ["1", "2"])
def test_complete_radial_integrals(power_value):
    with mp.workdps(70):
        m0 = mp.mpf(power_value)
        tail = mp.quad(lambda y: y**5 / (y * y + m0 * m0) ** 11, [0, m0, mp.inf])
        compact = mp.quad(
            lambda y: y**3 / (y * y + m0 * m0) ** (mp.mpf(21) / 2), [0, m0, mp.inf]
        )
        assert abs(tail / (1 / (720 * m0**16)) - 1) < mp.mpf("1e-60")
        assert abs(compact / (2 / (323 * m0**17)) - 1) < mp.mpf("1e-60")


def test_physical_curved_Ward_identity_with_external_mass():
    q, M, qd, Md, H = s.symbols("q M qdot Mdot H", real=True)
    c0, c1, c2, c3 = s.symbols("c0 c1 c2 c3", real=True)
    sx = s.Matrix([[0, 1], [1, 0]])
    sy = s.Matrix([[0, -s.I], [s.I, 0]])
    sz = s.diag(1, -1)
    c = (c0 * s.eye(2) + c1 * sx + c2 * sy + c3 * sz) / 2
    h = q * sx + M * sz
    cd = -s.I * (h * c - c * h)
    rho = s.trace(h * c)
    pressure = s.trace(q * sx * c) / 3
    dot = -3 * H * rho + s.trace((qd * sx + Md * sz) * c + h * cd)
    assert (
        s.simplify(
            (dot + 3 * H * (rho + pressure) - Md * s.trace(sz * c)).subs(qd, -H * q)
        )
        == 0
    )


def test_full_four_component_Sylvester_uniqueness():
    entries = s.symbols("d:16")
    D = s.Matrix(4, 4, entries)
    P = s.diag(1, 1, 0, 0)
    h = s.diag(1, 1, -1, -1)
    equations = list(P * D + D * P - D) + list(h * D - D * h)
    assert s.linsolve(equations, entries) == s.FiniteSet((0,) * 16)


@pytest.mark.parametrize("t", ["-2", "0", "0.4"])
@pytest.mark.parametrize("amplitude", ["-0.01", "0", "0.01"])
def test_physical_mass_rotation_and_connection_sign(t, amplitude):
    with mp.workdps(75):
        t = mp.mpf(t)
        m = mp.mpf(2)
        p = mp.mpf("1.3")
        d = mp.mpf(amplitude)
        tau = mp.mpf(".7")
        sx = mp.matrix([[0, 1], [1, 0]])
        sy = mp.matrix([[0, -1j], [1j, 0]])
        sz = mp.matrix([[1, 0], [0, -1]])
        mass = lambda z: m + d * (z / tau) / mp.power(1 + (z / tau) ** 8, mp.mpf(1) / 8)
        q = lambda z: p / (1 + z * z) ** 2
        angle = lambda z: mp.atan(q(z) / mass(z))
        U = lambda z: mp.eye(2) * mp.cos(angle(z) / 2) - 1j * sy * mp.sin(angle(z) / 2)
        h = q(t) * sx + mass(t) * sz
        value = U(t).H * h * U(t) - 1j * U(t).H * mp.diff(U, t)
        omega = mp.sqrt(q(t) ** 2 + mass(t) ** 2)
        g = q(t) * (mp.diff(mass, t) + 4 * t / (1 + t * t) * mass(t)) / (2 * omega**2)
        assert mp.norm(value - omega * sz - g * sy) < mp.mpf("1e-60")
        assert mp.norm(U(t).H * h * U(t) - omega * sz) < mp.mpf("1e-60")
        if t != 0:
            assert abs(q(t) * 4 * t / (1 + t * t) * mass(t) / (2 * omega**2)) > 0


@pytest.mark.parametrize("T", ["1", "2", "10"])
def test_complete_fixed_mode_Moller_tail_integral(T):
    with mp.workdps(70):
        T = mp.mpf(T)
        p = mp.mpf("1.3")
        d = mp.mpf(".01")
        tau = mp.mpf(".4")
        geometric = mp.quad(lambda t: p / (1 + t * t) ** 2, [T, 2 * T, mp.inf])
        mass_tail = mp.quad(
            lambda t: d * (-mp.expm1(-mp.log1p((tau / t) ** 8) / 8)), [T, 2 * T, mp.inf]
        )
        assert geometric + mass_tail <= p / (3 * T**3) + d * tau**8 / (56 * T**7)


def test_CD_not_D_and_state_scope_not_parent_solution():
    t = s.Symbol("t", real=True)
    d = geometry.data()
    assert s.simplify(d["scale_factor"] - (1 + t * t)) != 0
    assert d["witness_background"]["tail_power"] == "4"
    assert "NOT the geometric bounce time" in d["units"]
    assert (
        "not asserted"
        in Path(__file__).parents[1].joinpath("FORMULATION.md").read_text()
    )


def test_prescribed_curved_state_scope_and_strict_serializer():
    from p8_vacuum_curved_dirac_state import verify

    for mod in audit.MODULES:
        payload = {k: v for k, v in mod.data().items() if k != "checks"}
        json.dumps(verify.serialize(payload), sort_keys=True)
    with pytest.raises(ValueError):
        verify.serialize({"not_an_exact_payload": s.Float("0.1")})
    observable = audit.observable()
    assert "42copies" in observable["quantity"]
    assert "absolute" in observable["finite_reference"]
    assert "Hadamard" in hadamard.data()["same_in_out_state_identified"]
    assert (
        "not uniform"
        in Path(__file__).parents[1].joinpath("notes/geometry.md").read_text()
    )
    actual = geometry.data()
    assert (
        transition.bounds(1, actual["mean_mass"], 0, actual["mass_transition_time"])[
            "total"
        ]
        > 0
    )
