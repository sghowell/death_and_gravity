"""Independent coefficient, finite-q saddle, dynamical chart and energy controls."""

from functools import cache

import pytest
import sympy as s
from p8_affine import verify as serializer
from p8_vacuum_affine_scalar_tame_propagator import (
    audit,
    coefficients,
    energy,
    majorants,
)
from p8_vacuum_affine_scalar_tame_propagator import charts as c

PACKETS = (coefficients.data, c.data, majorants.data, energy.data)


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_integrated_exact_residual(name, value):
    row = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(s.factor(v) == 0 for v in row), name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_integrated_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_and_false_closure_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_complete_chart_cover_and_unchanged_frontier():
    for t in (-s.Rational(1, 2), 0, s.Rational(1, 2)):
        assert audit.require_scope(t)[0] == t
        assert audit.require_chart("low", t, 0) == ("low", t, 0)
        assert audit.require_chart("low", t, 100) == ("low", t, 100)
    for t in (-s.Rational(3, 16), s.Rational(3, 16)):
        for chart in ("outer", "central"):
            assert audit.require_chart(chart, t, 100) == (chart, t, 100)
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.frontier()) == 9
    assert audit.validate_scope(audit.frontier(), audit.matching())


@pytest.mark.parametrize("packet", PACKETS)
def test_exact_science_packet_and_serialization(packet):
    data = packet()
    for key, value in data["checks"].items():
        rows = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.factor(v) == 0 for v in rows), key
    assert all(bool(v) for v in data["gates"].values())
    serializer.serialize(
        {k: v for k, v in data.items() if k not in ("checks", "gates")}
    )


def test_independent_actual_principal_gradient_and_characteristics():
    t = s.Symbol("t", real=True)
    a = (1 + t * t) ** 2
    th = s.diff(a, t) / a - t / (1 + t * t) ** 4
    ell = 1 / (10 * (1 + t * t) ** 6)
    E = 1 - 3 / (2 * (1 + t * t) ** 3)
    F = s.factor(
        th * s.diff(E, t)
        - E * s.diff(th, t)
        + s.diff(a, t) / a * E * th
        - th * th
        - ell * ell * E * E / 2
    )
    actual = coefficients.clock()
    assert s.factor(F - actual["F"].subs(coefficients.t, t)) == 0
    assert (
        s.factor(actual["J"].subs(coefficients.t, t) - F - 1 / (50 * (1 + t * t) ** 6))
        == 0
    )
    for time in (
        -s.Rational(1, 2),
        -s.Rational(1, 4),
        0,
        s.Rational(1, 4),
        s.Rational(1, 2),
    ):
        f = F.subs(t, time)
        j = f + 1 / (50 * (1 + time * time) ** 6)
        assert 0 < f < j
    # The principal numerator remains nonzero at Theta=0.
    assert F.subs(t, 0) == s.Rational(1199, 800)
    assert th.subs(t, 0) == 0


@cache
def independent_fixture():
    """Independent smooth test family and original weighted canonical matrix.

    Nonzero rational retuning probes the complete algebra, not a numerical
    approximation to the extremely small physical QG1 expectation values.
    """
    t, p = s.symbols("clock transfer", real=True)
    a = (1 + t * t) ** 2
    H = s.diff(a, t) / a
    ell = 1 / (10 * (1 + t * t) ** 6)
    delta = 1 / (2 * (1 + t * t) ** 3)
    th = H - t / (1 + t * t) ** 4
    E = 1 - 3 * delta
    w = -ell * E
    F = th * s.diff(E, t) - E * s.diff(th, t) + H * E * th - th * th - w * w / 2
    Jbare = F + 1 / (50 * (1 + t * t) ** 6)
    rho = s.Rational(1, 10**8) * (1 + t + t * t)
    pressure = s.Rational(2, 10**8) * (1 - t + 3 * t * t)
    A = -pressure
    B = -(rho + pressure) / 2
    J = Jbare + (21 * delta * delta - 3 * delta) * A / 2 + (1 - 6 * delta) * B
    T = rho - 3 * delta * pressure
    q = p * p / a**2
    v, sigma, pv, ps = s.symbols("field matter momentum matter_momentum", real=True)
    Z = s.Matrix([v, sigma, pv, ps])
    lapse = (th * pv - w * ps + 3 * th * ell * sigma - (2 * E * q + 3 * T) * v) / (
        2 * J
    )
    hamiltonian = (
        ps * ps / 2 - ell * pv * sigma / 2 + (q / 2 - 3 * ell * ell / 4) * sigma * sigma
    )
    hamiltonian -= q * v * v + s.Rational(9, 2) * A * v * v
    hamiltonian += J * lapse * lapse
    flow = s.Matrix(
        [
            s.diff(hamiltonian, pv),
            s.diff(hamiltonian, ps),
            -s.diff(hamiltonian, v) - 3 * H * pv,
            -s.diff(hamiltonian, sigma) - 3 * H * ps,
        ]
    )
    matrix = flow.jacobian(Z)
    profiles = (th, E, ell, J, A, T, H)
    mapping = dict(zip(c.jetvars, profiles))
    mapping.update({key: s.diff(expr, t) for key, expr in zip(c.jets, profiles)})
    mapping.update({key: s.diff(expr, t, 2) for key, expr in zip(c.jet2, profiles)})
    mapping[c.q] = q
    return t, p, Z, matrix, mapping, lapse, hamiltonian


@pytest.mark.parametrize(
    "which,time",
    (
        ("central", s.Rational(-3, 16)),
        ("central", 0),
        ("central", s.Rational(3, 16)),
        ("outer", s.Rational(-2, 5)),
        ("outer", s.Rational(-3, 16)),
        ("outer", s.Rational(3, 16)),
        ("outer", s.Rational(2, 5)),
    ),
)
@pytest.mark.parametrize("momentum", (100, 300))
def test_independent_canonical_flow_equals_complete_chart_Euler(which, time, momentum):
    t, p, _Z, L, mapping, _, _ = independent_fixture()
    at = {t: time, p: momentum}
    actual = {key: s.factor(value.subs(at)) for key, value in mapping.items()}
    q = mapping[c.q]
    # Derive the coordinate map independently, including its q-dependent time derivative.
    Y = (
        s.Matrix([[0, 0, -1 / (2 * q), 0], [0, 1, 0, 0]])
        if which == "central"
        else s.Matrix([[1, 0, 0, 0], [0, 1, 0, 0]])
    )
    Y0 = Y.subs(at)
    Y1 = Y.diff(t).subs(at)
    Y2 = Y.diff(t, 2).subs(at)
    L0 = L.subs(at)
    L1 = L.diff(t).subs(at)
    velocity = Y1 + Y0 * L0
    acceleration = Y2 + 2 * Y1 * L0 + Y0 * (L1 + L0 * L0)
    chart = c.central() if which == "central" else c.outer()
    K = chart["K"].subs(actual)
    Kdot = chart["K"].applyfunc(c.dtime).subs(actual)
    gyro = chart["gyro"].subs(actual)
    potential = (mapping[c.q].subs(at) * chart["G"] + chart["lower"]).subs(actual)
    residual = (
        K * acceleration
        + (Kdot + 3 * actual[c.H] * K + gyro) * velocity
        + potential * Y0
    )
    assert all(s.cancel(value) == 0 for value in residual)


@pytest.mark.parametrize("time", (-s.Rational(3, 16), 0, s.Rational(3, 16)))
@pytest.mark.parametrize("momentum", (100, 300))
def test_independent_finite_q_auxiliary_saddle_and_phase_reconstruction(time, momentum):
    t, p, Z, flow, mapping, _, hamiltonian = independent_fixture()
    at = {t: time, p: momentum}
    actual = {key: s.factor(value.subs(at)) for key, value in mapping.items()}
    q = actual[c.q]
    v, sigma, pv, ps = Z
    b, bd, sd = s.symbols("independent_b independent_bd independent_sd", real=True)
    # Literal weighted canonical action, then solve its auxiliary stationarity.
    L = (
        2 * q * v * bd
        + 2 * actual[c.H] * q * v * b
        + ps * sd
        - hamiltonian.subs(at).subs(pv, -2 * q * b)
    )
    sol = s.solve([s.diff(L, v), s.diff(L, ps)], (v, ps), dict=True)[0]
    reduced = s.expand(L.subs(sol, simultaneous=True))
    yi = s.Matrix([b, sigma])
    dyi = s.Matrix([bd, sd])
    K = s.hessian(reduced, dyi)
    B = s.Matrix([[s.diff(reduced, d, y) for y in yi] for d in dyi])
    D = s.hessian(reduced, yi)
    for name, observed in (("K", K), ("B", B), ("D", D)):
        assert all(s.cancel(v) == 0 for v in observed - c.central()[name].subs(actual))
    for values in (
        (s.Rational(1, 7), s.Rational(-2, 9), s.Rational(3, 11), s.Rational(-5, 13)),
        (0, 1, 2, 0),
    ):
        zdot = flow.subs(at) * s.Matrix(values)
        bval = -values[2] / (2 * q)
        bdval = -zdot[2] / (2 * q) + 2 * actual[c.H] * bval
        fields = {b: bval, bd: bdval, sd: zdot[1], sigma: values[1]}
        assert s.cancel(sol[v].subs(fields) - values[0]) == 0
        assert s.cancel(sol[ps].subs(fields) - values[3]) == 0


def test_finite_q_chart_failure_does_not_singularize_original_Hamiltonian():
    _, _, Z, _, _, _, _ = independent_fixture()
    # Independent bare bounce Hessian has auxiliary determinant zero at152/25.
    q = s.Symbol("q", positive=True)
    v, sg, pv, ps = Z
    J = s.Rational(243, 160)
    ell = s.Rational(1, 10)
    E = s.Rational(-1, 2)
    w = -ell * E
    b = s.Symbol("b", real=True)
    H = (
        ps * ps / 2
        - ell * pv * sg / 2
        + (q / 2 - 3 * ell * ell / 4) * sg * sg
        - q * v * v
        + (-w * ps - 2 * E * q * v) ** 2 / (4 * J)
    )
    aux = s.hessian(H.subs(pv, -2 * q * b), (v, ps))
    qstar = s.Rational(152, 25)
    assert s.factor(aux.det().subs(q, qstar)) == 0
    original = s.hessian(H, Z).subs(q, qstar)
    assert all(value.is_finite for value in original)
    assert s.factor(aux.det().subs(q, 4096)) > 0


def test_deleting_weighted_transport_changes_complete_dynamics():
    chart = c.central()
    y = s.Matrix([c.b, c.sigma])
    yd = s.Matrix([c.bd, c.sd])
    rhs = chart["M"] * yd - chart["XY"] * y
    wrong = (rhs.T * chart["inverse"] * rhs)[0] / 2 - (y.T * chart["YY"] * y)[0] / 2
    good = (
        (yd.T * chart["K"] * yd)[0] / 2
        + (yd.T * chart["B"] * y)[0]
        + (y.T * chart["D"] * y)[0] / 2
    )
    probe = {
        c.th: s.Rational(1, 2),
        c.E: s.Rational(-1, 3),
        c.l: s.Rational(1, 10),
        c.J: 2,
        c.A: s.Rational(1, 10**7),
        c.T: s.Rational(-1, 10**7),
        c.H: s.Rational(1, 4),
        c.q: 5000,
        c.b: s.Rational(1, 3),
        c.sigma: s.Rational(1, 5),
        c.bd: s.Rational(2, 7),
        c.sd: s.Rational(3, 11),
    }
    assert s.factor((good - wrong).subs(probe)) != 0


def test_independent_energy_differentiation_on_nonsymmetric_lower_term():
    t = s.Symbol("test_time", real=True)
    y = s.Matrix([1 + t + t * t, 2 - t + 3 * t**3])
    dy = y.diff(t)
    ddy = dy.diff(t)
    q = s.exp(-2 * t)
    K = s.Matrix([[2 + t * t, t / 10], [t / 10, 3 + t * t]])
    G = s.Matrix([[4 + t * t, t / 7], [t / 7, 2 + t * t]])
    gyro = s.Matrix([[0, 1 + t], [-1 - t, 0]])
    # Choose forcing so this arbitrary curve solves the exact Euler equation.
    R = s.Matrix([[1, 2], [-3, 4]])
    forcing = K * ddy + (K.diff(t) + 3 * K + gyro) * dy + (q * G + R) * y
    E = (dy.T * K * dy + q * y.T * G * y)[0] / 2
    rhs = -(dy.T * K.diff(t) * dy)[0] / 2 - 3 * (dy.T * K * dy)[0] - (dy.T * R * y)[0]
    rhs += q * (y.T * (G.diff(t) - 2 * G) * y)[0] / 2 + (dy.T * forcing)[0]
    assert s.simplify(s.diff(E, t) - rhs) == 0
    assert R != R.T


def test_continuum_weight_and_no_same_space_inference():
    r = s.Symbol("r", real=True)
    # lambda^6 in phase amplitude requires twelve, not six, Sobolev derivatives.
    assert (r + 12) - r == 2 * 6
    assert energy.SPATIAL_LOSS == 12
    assert energy.LOG_PROPAGATOR == s.Integer(10) ** 29
    d = energy.data()
    assert "not small relative to kappa" in d["strict_scope"]
    assert "not an unproved map" in d["Duhamel_Sobolev_statement"]
