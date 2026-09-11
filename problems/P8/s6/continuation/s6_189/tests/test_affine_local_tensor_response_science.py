"""Independent rotating TT curvature, actual action and full-spacetime norm fixtures."""

from functools import cache

import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_local_tensor_response import (
    audit,
    bounds,
    geometry,
    operator,
    prescription,
)
from scipy.integrate import quad


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_exact_residual(name):
    assert audit.residuals()[name] == 0, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_unsupported_scope_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("name", list(audit.packets()))
def test_all_continuous_proof_packets(name):
    assert all(audit.packets()[name].get("gates", {}).values())


def test_scope_and_counts():
    assert len(audit.residuals()) == 34 and audit.scalar_entry_count() == 34
    assert len(audit.gates()) == 28 and audit.rejected_inputs() == 111
    assert len(audit.controls()) == 9 and all(audit.gates().values())
    assert audit.validate_scope(audit.frontier(), audit.matching())


@pytest.mark.parametrize("fixture", range(10))
def test_independent_full_spatial_rotated_two_polarization_Weyl(fixture):
    rng = np.random.default_rng(16000 + fixture)
    k = rng.normal(size=3) * (1 + fixture / 3)
    n = k / np.linalg.norm(k)
    e = rng.normal(size=3)
    e -= n * np.dot(n, e)
    e /= np.linalg.norm(e)
    f = np.cross(n, e)
    plus = (np.outer(e, e) - np.outer(f, f)) / np.sqrt(2)
    cross = (np.outer(e, f) + np.outer(f, e)) / np.sqrt(2)
    assert np.allclose(np.sum(plus * plus), 1) and np.allclose(np.sum(cross * cross), 1)
    jets = rng.normal(
        size=(2, 3)
    )  # tt, t-theta, theta-theta for each physical unit polarization
    second = np.zeros((4, 4, 4, 4))
    for E, (tt, tx, xx) in zip((plus, cross), jets):
        second[0, 0, 1:, 1:] -= tt * E
        for i in range(3):
            second[0, i + 1, 1:, 1:] -= k[i] * tx * E
            second[i + 1, 0, 1:, 1:] -= k[i] * tx * E
            for j in range(3):
                second[i + 1, j + 1, 1:, 1:] -= k[i] * k[j] * xx * E
    R = np.zeros((4, 4, 4, 4))
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    R[a, b, c, d] = (
                        second[c, b, a, d]
                        + second[d, a, b, c]
                        - second[c, a, b, d]
                        - second[d, b, a, c]
                    ) / 2
    sign = np.array([1.0, -1.0, -1.0, -1.0])
    Ric = np.array(
        [
            [sum(sign[a] * R[a, b, a, d] for a in range(4)) for d in range(4)]
            for b in range(4)
        ]
    )
    scalar = np.dot(sign, np.diag(Ric))
    Ric2 = sum(sign[a] * sign[b] * Ric[a, b] ** 2 for a in range(4) for b in range(4))
    Riem2 = sum(
        sign[a] * sign[b] * sign[c] * sign[d] * R[a, b, c, d] ** 2
        for a in range(4)
        for b in range(4)
        for c in range(4)
        for d in range(4)
    )
    Weyl = Riem2 - 2 * Ric2 + scalar * scalar / 3
    k2 = np.dot(k, k)
    expected = sum(((tt + k2 * xx) ** 2 - 4 * k2 * tx**2) / 2 for tt, tx, xx in jets)
    box = sum((tt - k2 * xx) ** 2 / 2 for tt, _, xx in jets)
    boundary = sum(2 * k2 * (tt * xx - tx**2) for tt, tx, xx in jets)
    assert abs(scalar) < 1e-10
    assert np.isclose(Weyl, expected, rtol=1e-10, atol=1e-10)
    assert np.isclose(Weyl - box, boundary, rtol=1e-10, atol=1e-10)


@cache
def compact_functions(kind):
    t = operator.t
    radius = s.Rational(2, 5) if kind == 0 else s.Rational(7, 20)
    core = s.exp(-1 / (1 - (t / radius) ** 2))
    field = (
        core * s.cos(s.Rational(3, 2) * t)
        if kind == 0
        else core * s.sin(t + s.Rational(3, 10))
    )
    functions = [s.lambdify(t, s.diff(field, t, j), "numpy") for j in range(5)]
    return float(radius), functions


def values(time, kind=0):
    radius, functions = compact_functions(kind)
    return (
        np.array([float(f(time)) for f in functions])
        if abs(time) < radius
        else np.zeros(5)
    )


def coefficients(time):
    a = (1 + time * time) ** 2
    H = 4 * time / (1 + time * time)
    Hp = 4 * (1 - time * time) / (1 + time * time) ** 2
    Hpp = -8 * time * (3 - time * time) / (1 + time * time) ** 3
    R = 24 * (1 + 7 * time * time) / (1 + time * time) ** 2
    Rp = 48 * time * (5 - 7 * time * time) / (1 + time * time) ** 3
    A = 5 * 1000**2 / 12 - R / 36
    Ap = -Rp / 36
    return a, H, Hp, Hpp, A, Ap


def quadratic_density(time, k, h):
    a, H, _, _, A, _ = coefficients(time)
    q = k * k / a**2
    D = h[2] + H * h[1] + q * h[0]
    return a**3 * (A * (h[1] ** 2 - q * h[0] ** 2) - D * D / 60)


def local_readout(time, k, h):
    a, H, Hp, Hpp, A, Ap = coefficients(time)
    q = k * k / a**2
    fourth = (
        h[4]
        + 6 * H * h[3]
        + (4 * Hp + 11 * H * H + 2 * q) * h[2]
        + (Hpp + 7 * H * Hp + 6 * H**3 + 2 * H * q) * h[1]
        + q * q * h[0]
    )
    return A * (h[2] + 3 * H * h[1] + q * h[0]) + Ap * h[1] + fourth / 60


@pytest.mark.parametrize("k", (0.0, 0.1, 1.0, 10.0, 100.0, 1000.0))
def test_independent_compact_action_variation_and_weighted_adjoint(k):
    epsilon = 1e-5
    lhs = quad(
        lambda time: (
            (
                quadratic_density(time, k, values(time) + epsilon * values(time, 1))
                - quadratic_density(time, k, values(time) - epsilon * values(time, 1))
            )
            / (2 * epsilon)
        ),
        -0.35,
        0.35,
        epsabs=1e-3,
        epsrel=2e-8,
        limit=300,
    )[0]
    rhs = quad(
        lambda time: (
            -2
            * coefficients(time)[0] ** 3
            * values(time, 1)[0]
            * local_readout(time, k, values(time))
        ),
        -0.35,
        0.35,
        epsabs=1e-5,
        epsrel=2e-9,
        limit=300,
    )[0]
    assert np.isclose(lhs, rhs, rtol=2e-7, atol=0.02)
    action = quad(
        lambda time: quadratic_density(time, k, values(time)),
        -0.4,
        0.4,
        epsabs=1e-4,
        epsrel=1e-9,
        limit=300,
    )[0]
    paired = quad(
        lambda time: (
            coefficients(time)[0] ** 3
            * values(time)[0]
            * local_readout(time, k, values(time))
        ),
        -0.4,
        0.4,
        epsabs=1e-4,
        epsrel=1e-9,
        limit=300,
    )[0]
    assert np.isclose(action, -paired, rtol=1e-8, atol=0.01)


@pytest.mark.parametrize("k", (0.0, 1e-6, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0))
def test_full_spatial_joint_H4_mode_norm(k):
    # A modal fixture uses the z direction; the continuous proof integrates all Fourier directions on R3.
    grid = np.linspace(-0.4, 0.4, 3001)
    norm = []
    readout = []
    for time in grid:
        h = values(time)
        norm.append(
            sum(h[j] ** 2 * sum(k ** (2 * l) for l in range(5 - j)) for j in range(5))
        )
        readout.append((local_readout(time, k, h) / (8 * np.pi**2)) ** 2)
    H4 = np.trapezoid(norm, grid)
    output = np.trapezoid(readout, grid)
    assert output < (60000**2) * H4
    assert bounds.data()["local_response_upper_before_inverse_kappa"] < 60000


@pytest.mark.parametrize("time", (-0.5, -0.4, -0.2, 0.0, 0.2, 0.4, 0.5))
def test_actual_continuous_coefficient_enclosures(time):
    _, H, Hp, Hpp, A, Ap = coefficients(time)
    assert abs(H) <= 1.6 + 1e-14 and abs(Hp) <= 4 + 1e-14
    assert abs(Hpp) <= 12 and 0 < A < 500000 and abs(Ap) < 4


def test_wrong_unweighted_adjoint_has_nonzero_bounce_error():
    t = operator.t
    h = s.Integer(1)
    q = 1 / operator.a**2
    D = lambda f: s.diff(f, t, 2) + operator.H * s.diff(f, t) + q * f
    Dadj = lambda f: (
        s.diff(f, t, 2)
        + 5 * operator.H * s.diff(f, t)
        + (2 * s.diff(operator.H, t) + 6 * operator.H**2 + q) * f
    )
    assert s.simplify((Dadj(D(h)) - D(D(h))).subs(t, 0)) == 8


def test_fourth_derivative_cannot_be_discarded():
    # At the bounce this independent fourth jet isolates the Weyl contribution.
    h = np.array([0.0, 0.0, 0.0, 0.0, 1.0])
    assert np.isclose(local_readout(0.0, 0.0, h), 1 / 60)
    assert "not the full determinant response" in operator.data()["response_boundary"]


def test_no_new_prescription_or_complete_response_claim():
    p = prescription.data()
    assert "not the full state-dependent determinant" in p["same_prescription"]
    assert "not changed" in p["fixed_scalar_profile"]
    assert "both orthonormal TT polarizations" in geometry.data()["tensor_Weyl_rule"]
