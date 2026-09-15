"""Independent geometric, clock, operator and scope diagnostics."""

import itertools

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_finite_volume_turnaround import (
    audit,
    branch,
    geometry,
    moving,
    quantum,
    source,
)

pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning")
PACKETS = (source.data, moving.data, geometry.data, branch.data, quantum.data)


@pytest.mark.parametrize("packet", PACKETS)
def test_entire_exact_packet(packet):
    data = packet()
    for name, value in data["checks"].items():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.factor(entry) == 0 for entry in entries), name
    assert all(bool(value) for value in data["gates"].values())


@pytest.mark.parametrize("index", range(187))
def test_each_complete_translated_derivative_bound(index):
    name, expression = list(source.expressions().items())[index]
    value = source.magnitude(source.q.eliminate_N_primitive(expression))
    assert value == source.derivative_bounds()[name]
    assert value.is_Rational is True and 0 <= value < source.MAJORANT


@pytest.mark.parametrize("order", range(196))
def test_each_actual_longer_domain_phase_jet_ratio(order):
    ratio = quantum.jet(order + 1, quantum.VERROR) / quantum.jet(order, quantum.VERROR)
    assert ratio == 2056 * (order + 1) ** 3 / moving.CORE
    assert ratio < s.Rational(1, 10**9)


@pytest.mark.parametrize("z", source.q.COORDS)
def test_each_full_moving_center_gradient_by_independent_difference(z):
    # Formal full-source polynomial/analytic differentiation is tested before
    # clock restriction, not after replacing the implicit source by a germ.
    shift = s.Symbol("independent_invariant_shift", real=True)
    expression = source.q.CONSTRAINT.subs(z, z + shift)
    direct = s.diff(expression, shift).subs(shift, 0)
    assert s.expand(direct - s.diff(source.q.CONSTRAINT, z)) == 0
    assert source.clock_at(direct) == source.center()["gradients"][z]


@pytest.mark.parametrize(
    "value",
    (s.Rational(-1, 3), s.Rational(-1, 10), 0, s.Rational(1, 10), s.Rational(1, 3)),
)
def test_exact_clock_center_algebra_not_additional_quantitative_time_domain(value):
    rho, pressure = source.q.physical.RHO, source.q.physical.PRESSURE
    actual = source.center()["C"]
    expected = 3 * pressure / (2 * (1 + source.u**2) ** 3) - rho
    assert s.factor(actual - expected) == 0
    assert actual.subs({rho: 0, pressure: 0}).subs(source.u, value) == 0


@pytest.mark.parametrize("normal", ((1, 0, 0), (1, 2, 2), (2, -1, 2), (1, -2, -2)))
def test_full_TT_linear_scalar_curvature_vanishes(normal):
    k = s.Matrix(normal)
    P = s.eye(3) - k * k.T / k.dot(k)
    seed = s.Matrix([[1, 2, -1], [2, 3, 4], [-1, 4, -2]])
    tau = P * seed * P - P * s.trace(P * seed) / 2
    assert tau * k == s.zeros(3, 1)
    assert s.trace(tau) == 0
    assert -(k.T * tau * k)[0] + k.dot(k) * s.trace(tau) == 0


def numeric_full_curvature(normal, amplitude, conformal, point):
    """Direct full three-dimensional metric jets, Christoffel and Ricci."""
    k = mp.matrix(normal)
    khat = k / mp.sqrt((k.T * k)[0])
    e0 = mp.matrix([1, 0, 0])
    difference = e0 - khat
    O = (
        mp.eye(3)
        if mp.norm(difference) == 0
        else mp.eye(3) - 2 * difference * difference.T / (difference.T * difference)[0]
    )
    ell = mp.matrix([2, -1, 3])
    x = mp.matrix(point)
    angle, phase = (k.T * x)[0], (ell.T * x)[0]
    v = conformal * mp.sin(phase)
    vi = [conformal * mp.cos(phase) * ell[i] for i in range(3)]
    vij = [
        [-conformal * mp.sin(phase) * ell[i] * ell[j] for j in range(3)]
        for i in range(3)
    ]
    f = amplitude * mp.cos(angle)
    fi = [-amplitude * mp.sin(angle) * k[i] for i in range(3)]
    fij = [
        [-amplitude * mp.cos(angle) * k[i] * k[j] for j in range(3)] for i in range(3)
    ]
    eigen = (0, -1, 1)

    def diagonal(values):
        return mp.diag(values)

    D = diagonal([mp.exp(sign * f) for sign in eigen])
    Qinv = O * D * O.T
    dQ = [
        O * diagonal([sign * mp.exp(sign * f) * fi[i] for sign in eigen]) * O.T
        for i in range(3)
    ]
    ddQ = [
        [
            O
            * diagonal(
                [
                    mp.exp(sign * f) * (sign * fij[i][j] + sign**2 * fi[i] * fi[j])
                    for sign in eigen
                ]
            )
            * O.T
            for j in range(3)
        ]
        for i in range(3)
    ]
    scale = mp.exp(2 * v)
    g = scale * Qinv
    dg = [scale * (2 * vi[i] * Qinv + dQ[i]) for i in range(3)]
    ddg = [
        [
            scale
            * (
                (4 * vi[i] * vi[j] + 2 * vij[i][j]) * Qinv
                + 2 * vi[i] * dQ[j]
                + 2 * vi[j] * dQ[i]
                + ddQ[i][j]
            )
            for j in range(3)
        ]
        for i in range(3)
    ]
    inv = g**-1
    dinv = [-inv * dg[i] * inv for i in range(3)]
    Gamma = [
        [
            [
                mp.fsum(
                    inv[a, l] * (dg[i][j, l] + dg[j][i, l] - dg[l][i, j])
                    for l in range(3)
                )
                / 2
                for j in range(3)
            ]
            for i in range(3)
        ]
        for a in range(3)
    ]
    dGamma = [
        [
            [
                [
                    mp.fsum(
                        dinv[m][a, l] * (dg[i][j, l] + dg[j][i, l] - dg[l][i, j])
                        + inv[a, l]
                        * (ddg[m][i][j, l] + ddg[m][j][i, l] - ddg[m][l][i, j])
                        for l in range(3)
                    )
                    / 2
                    for j in range(3)
                ]
                for i in range(3)
            ]
            for a in range(3)
        ]
        for m in range(3)
    ]
    Ricci = mp.matrix(3, 3)
    for i, j in itertools.product(range(3), repeat=2):
        Ricci[i, j] = mp.fsum(dGamma[a][a][i][j] - dGamma[j][a][i][a] for a in range(3))
        Ricci[i, j] += mp.fsum(
            Gamma[a][a][b] * Gamma[b][i][j] - Gamma[a][j][b] * Gamma[b][i][a]
            for a in range(3)
            for b in range(3)
        )
    actual = mp.fsum(inv[i, j] * Ricci[i, j] for i in range(3) for j in range(3))
    # Independent exact warped/conformal formula. det Q=1 and div Q=0
    # hold for this whole exponential shape, including all harmonics.
    Q = Qinv**-1
    ell_square = (ell.T * Q * ell)[0]
    reference = mp.exp(-2 * v) * (
        -(amplitude**2) * (k.T * k)[0] * mp.sin(angle) ** 2 / 2
        + 4 * conformal * mp.sin(phase) * ell_square
        - 2 * conformal**2 * mp.cos(phase) ** 2 * ell_square
    )
    return actual, reference, g, Q


@pytest.mark.parametrize("normal", ((1, 0, 0), (1, 2, 2), (2, -1, 2)))
@pytest.mark.parametrize(
    "amplitude,conformal",
    (("0.003", "0.002"), ("-0.004", "0.001"), ("0.005", "-0.003"), ("0.002", "0")),
)
def test_literal_rotated_full_metric_curvature(normal, amplitude, conformal):
    with mp.workdps(65):
        actual, expected, g, Q = numeric_full_curvature(
            normal,
            mp.mpf(amplitude),
            mp.mpf(conformal),
            [mp.mpf("0.2"), mp.mpf("-0.3"), mp.mpf("0.1")],
        )
        assert abs(actual - expected) < mp.mpf("1e-52")
        assert abs(mp.det(Q) - 1) < mp.mpf("1e-55")
        assert mp.det(g) > 0
        if conformal == "0":
            assert actual < 0  # Nonlinear tensor curvature cannot be dropped.


@pytest.mark.parametrize("degree", range(1, 9))
def test_full_additive_sublattice_support_without_deleting_means(degree):
    basis = {
        tuple(sign * int(i == axis) for i in range(3))
        for axis in range(3)
        for sign in (-1, 1)
    }
    support = {(0, 0, 0)}
    for _ in range(degree):
        support |= {
            tuple(x[i] + y[i] for i in range(3)) for x in support for y in basis
        }
    assert (0, 0, 0) in support
    for n in support:
        k = tuple(moving.P * x for x in n)
        assert all(x % moving.P == 0 for x in k)
        if any(k):
            assert sum(x * x for x in k) >= moving.P**2


@pytest.mark.parametrize(
    "time,error,has_gap",
    (
        (source.TIME, quantum.OPERATOR_ERROR, True),
        (source.TIME, s.Rational(1, 10**200), False),
        (s.Rational(1, 10**2000), s.Rational(1, 10**199), False),
        (source.TIME, s.Rational(1, 10**320), True),
    ),
)
def test_endpoint_gap_requires_error_smaller_than_time_squared(time, error, has_gap):
    gap = (1 + time * time) ** 6 * (1 - error) - (1 + error)
    assert bool(gap > 0) is has_gap


@pytest.mark.parametrize("phase", (0, s.pi / 2, s.pi, 3 * s.pi / 2))
def test_noncommuting_unitary_need_not_be_near_identity(phase):
    X = s.Matrix([[0, 1], [1, 0]])
    U = s.cos(phase) * s.eye(2) - s.I * s.sin(phase) * X
    psi = s.Matrix([1, 0])
    F = s.eye(2) + quantum.OPERATOR_ERROR * s.diag(1, -1)
    state = U * psi
    assert U.conjugate().T * U == s.eye(2)
    mean = s.simplify((state.conjugate().T * F * state)[0])
    assert 1 - quantum.OPERATOR_ERROR <= mean <= 1 + quantum.OPERATOR_ERROR
    assert (1 + source.TIME**2) ** 6 * mean > 1 + quantum.OPERATOR_ERROR


@pytest.mark.parametrize(
    "name,call,args",
    audit.bad_cases(),
    ids=lambda value: value if isinstance(value, str) else None,
)
def test_rejected_scope_and_unproved_input(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_long_interval_does_not_extend_old_state_leakage_bounds():
    b = quantum.bounds()
    assert source.TIME * b["HC"] > 2 and source.TIME * b["HW"] > 2
    assert "EXTERNAL" in audit.observable()["not_established"]
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 128


def test_new_time_guard_accepts_only_the_stated_exact_interval():
    for value in (-source.TIME, 0, source.TIME):
        assert audit.require_time(value) == value
    with pytest.raises(ValueError):
        audit.require_time(2 * source.TIME)
