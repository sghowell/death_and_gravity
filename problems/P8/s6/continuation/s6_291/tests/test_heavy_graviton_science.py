"""Independent production tensors, full phase normalization and threshold limits."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_heavy_graviton_production_threshold import (
    audit,
    production,
    source,
    threshold,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not value.atoms(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_every_unsupported_input(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("dimension", (4, 5, 6))
@pytest.mark.parametrize("choice", (1, 2, 3))
def test_independent_all_components_four_graph_Ward_and_TT(dimension, choice):
    eta = s.diag(1, *([-1] * (dimension - 1)))
    E, a, b, w = (
        s.Integer(choice + 3),
        s.Rational(choice, 3),
        s.Rational(1, 2),
        s.Rational(1, 50),
    )
    p = s.zeros(dimension, 1)
    p[0] = E
    p[1] = a
    p[-1] = b
    r = s.zeros(dimension, 1)
    r[0] = E
    r[1] = -a
    r[-1] = -b
    q = s.zeros(dimension, 1)
    q[0] = w
    q[-1] = w
    dot = lambda x, y: (x.T * eta * y)[0]
    P = p + r
    H = P - q
    mu = dot(p, p)
    n = dot(H, H)
    energy = dot(P, P)
    assert mu > 0 and n > 4 * mu and energy > n
    vertex = lambda x, y, m: x * y.T + y * x.T - eta * (dot(x, y) - m)
    graphs = (
        vertex(p, p - q, mu) / (dot(p - q, p - q) - mu),
        vertex(r, r - q, mu) / (dot(r - q, r - q) - mu),
        vertex(P, H, n) / (energy - n),
        eta,
    )
    full = sum(graphs, s.zeros(dimension))
    assert full * eta * q == s.zeros(dimension, 1)
    assert (full - graphs[2]) * eta * q == -H
    assert (full - graphs[3]) * eta * q == -q
    d = dimension - 2
    unit = s.eye(d)
    projector = s.Matrix(
        d * d,
        d * d,
        lambda i, j: (
            (
                unit[i // d, j // d] * unit[i % d, j % d]
                + unit[i // d, j % d] * unit[i % d, j // d]
            )
            / 2
            - unit[i // d, i % d] * unit[j // d, j % d] / d
        ),
    )
    assert projector**2 == projector
    assert s.trace(projector) == s.Rational(dimension * (dimension - 3), 2)
    vector = s.Matrix(list(full[1:-1, 1:-1]))
    actual = (vector.T * projector * vector)[0]
    angle = b / s.sqrt(a * a + b * b)
    expected = production.polarization_sum(energy, angle, mu, n, 1, 1, dimension)
    assert s.simplify(actual - expected) == 0
    assert actual > 0
    if dimension == 4:
        for sign in (-1, 1):
            eps = s.Matrix([0, 1, sign * s.I, 0]) / s.sqrt(2)
            hel = (eps.T * eta * full * eta * eps)[0]
            assert (
                s.simplify(hel - production.helicity(energy, angle, mu, n, 1, 1)) == 0
            )
        assert actual == 2 * hel * hel


@pytest.mark.parametrize(
    "beta2", (s.S.Zero, s.Rational(1, 4), s.Rational(3, 4), s.S.One)
)
@pytest.mark.parametrize("epsilon", (s.S.Zero, s.Rational(1, 16), s.Rational(1, 8)))
def test_independent_whole_angular_integral_hypergeometric(beta2, epsilon):
    with mp.workdps(60):
        b = mp.mpf(str(beta2.p)) / int(beta2.q)
        e = mp.mpf(str(epsilon.p)) / int(epsilon.q)
        numerator = mp.quad(
            lambda x: (1 - x * x) ** (2 + e) / (1 - b * x * x) ** 2, [0, 1]
        )
        measure = mp.quad(lambda x: (1 - x * x) ** e, [0, 1])
        pref = 4 * (e + 1) * (e + 2) / ((2 * e + 3) * (2 * e + 5))
        expected = pref * mp.hyp2f1(2, mp.mpf("0.5"), e + mp.mpf("3.5"), b)
        assert abs(numerator / measure - expected) < mp.mpf("1e-48")
        assert pref - mp.mpf("1e-48") <= expected <= 1 + mp.mpf("1e-48")


@pytest.mark.parametrize(
    "mu,n,energy", ((1, 9, 10), (1, 9, 40), (2, 17, 20), (2, 17, 90))
)
@pytest.mark.parametrize("epsilon", (s.Rational(1, 16), s.Rational(1, 8)))
def test_independent_full_D_radial_phase_vs_B0_and_compact_cut(mu, n, energy, epsilon):
    with mp.workdps(60):
        mu, n, energy = map(mp.mpf, (mu, n, energy))
        e = mp.mpf(str(epsilon.p)) / int(epsilon.q)
        D = 4 + 2 * e
        nu2 = mp.mpf(3)
        beta = (energy - n) / energy
        p = (energy - n) / (2 * mp.sqrt(energy))
        sphere = 2 * mp.pi ** ((D - 1) / 2) / mp.gamma((D - 1) / 2)
        radial = (
            (2 * mp.pi) ** (2 - D)
            * sphere
            * p ** (D - 3)
            / (4 * mp.sqrt(energy) * nu2**e)
        )
        imB = (
            -mp.gamma(-e)
            * mp.sin(mp.pi * e)
            * (4 * mp.pi * nu2) ** (-e)
            * mp.quad(lambda x: (energy * x * (beta - x)) ** e, [0, beta])
        )
        assert abs(radial / 2 - imB / (16 * mp.pi**2)) < mp.mpf("1e-48")
        b = 1 - 4 * mu / energy
        numerator = mp.quad(
            lambda x: (1 - x * x) ** (e + 2) / (1 - b * x * x) ** 2, [0, 1]
        )
        measure = mp.quad(lambda x: (1 - x * x) ** e, [0, 1])
        literal = (
            radial
            / 2
            * 4
            * (D - 3)
            / (D - 2)
            * (energy - 4 * mu) ** 2
            / (energy - n) ** 2
            * numerator
            / measure
        )
        compact = (
            (energy - n) ** (-1 + 2 * e)
            * (energy - 4 * mu) ** 2
            / (4 * mp.pi * energy)
            * (1 + 2 * e)
            / (2 + 2 * e)
            * (16 * mp.pi * nu2 * energy) ** (-e)
            / mp.gamma(1 + e)
            * numerator
        )
        assert abs(literal - compact) < mp.mpf("1e-48")
        assert abs(literal - compact / 2) > mp.mpf("1e-6")


@pytest.mark.parametrize("mu,n,nu2", ((1, 9, 1), (1, 30, 2), (2, 17, 3), (2, 100, 1)))
def test_independent_whole_evanescent_derivative(mu, n, nu2):
    with mp.workdps(60):
        mu, n, nu2 = map(mp.mpf, (mu, n, nu2))
        b = 1 - 4 * mu / n
        N = lambda e: mp.quad(
            lambda x: (1 - x * x) ** (e + 2) / (1 - b * x * x) ** 2, [0, 1]
        )
        K = lambda e: (
            (n - 4 * mu) ** 2
            / (4 * mp.pi * n)
            * (1 + 2 * e)
            / (2 + 2 * e)
            * (16 * mp.pi * nu2 * n) ** (-e)
            / mp.gamma(1 + e)
            * N(e)
        )
        Llog = mp.quad(
            lambda x: (1 - x * x) ** 2 * mp.log(1 - x * x) / (1 - b * x * x) ** 2,
            [0, 1],
        )
        expected = (
            (n - 4 * mu) ** 2
            / (8 * mp.pi * n)
            * ((mp.euler + 1 - mp.log(16 * mp.pi * nu2 * n)) * N(0) + Llog)
        )
        assert abs(mp.diff(K, 0) - expected) < mp.mpf("1e-45")
        assert Llog < 0 and -2 < Llog
        assert abs(mp.diff(K, 0)) < 1000 * K(0)


@pytest.mark.parametrize(
    "beta2", (s.S.Zero, s.Rational(1, 3), s.Rational(8, 9), s.S.One)
)
@pytest.mark.parametrize("angle", (s.S.Zero, s.Rational(1, 3), s.Rational(3, 4)))
def test_independent_full_pointwise_kernel_bounds(beta2, angle):
    kernel = production.angular_kernel(beta2, angle)
    assert (1 - angle * angle) ** 2 <= kernel <= 1
    assert not kernel.has(s.Float)


@pytest.mark.parametrize(
    "mu,n,window",
    ((1, 9, s.Rational(1, 4)), (2, 17, s.Rational(1, 2)), (1, 100, s.Rational(1, 10))),
)
def test_exact_public_window_domain_and_positive_threshold(mu, n, window):
    assert threshold.require_window(mu, n, window) == (mu, n, window, 1)
    coefficient = threshold.leading_coefficient(n, mu, 1, 1)
    assert coefficient.is_positive
    assert not coefficient.has(s.Float)
    assert production.forward_cut(n + 1, mu, n, 1, 1).is_positive


@pytest.mark.parametrize(
    "args",
    (
        (1, 4, s.Rational(1, 4)),
        (1, 3, s.Rational(1, 4)),
        (1, 9, 1),
        (0, 9, s.Rational(1, 4)),
        (1, 9, 0),
        (1.0, 9, s.Rational(1, 4)),
        (True, 9, s.Rational(1, 4)),
    ),
)
def test_invalid_threshold_domains_rejected(args):
    with pytest.raises((TypeError, ValueError)):
        threshold.require_window(*args)


@pytest.mark.parametrize("dimension", (3, 0, -1, 4.0, True, s.Symbol("D")))
def test_nonphysical_matrix_dimension_rejected(dimension):
    with pytest.raises((TypeError, ValueError)):
        source.require_dimension(dimension)


@pytest.mark.parametrize("choice", (1, 2, 3))
def test_independent_integrated_smooth_polynomial_Laurent_model(choice):
    e, L = s.symbols("epsilon window", positive=True)
    y = s.Symbol("dimensionless_threshold", positive=True)
    A = s.Integer(choice)
    B = s.Rational(choice, 7)
    C = s.Rational(choice, 11)
    # A complete smooth, nonconstant weight with an independent epsilon derivative.
    H = A + e * B + e * e * C / 2 + (2 * A + e * B) * y + 3 * C * y * y
    whole = (
        (A + e * B + e * e * C / 2) * L ** (2 * e) / (2 * e)
        + (2 * A + e * B) * L ** (1 + 2 * e) / (1 + 2 * e)
        + 3 * C * L ** (2 + 2 * e) / (2 + 2 * e)
    )
    finite = A * s.log(L) + B / 2 + 2 * A * L + 3 * C * L * L / 2
    assert s.simplify(s.limit(whole - A / (2 * e), e, 0) - finite) == 0
    assert s.factor((H - H.subs(y, 0)) / y - (2 * A + e * B + 3 * C * y)) == 0
    # A D0-only numerator loses the B/2 finite term.
    assert B / 2 != 0


def test_every_scientific_packet_is_recursively_exact_and_ancestry_unchanged():
    def exact(value):
        if isinstance(value, float):
            pytest.fail("Python float in scientific payload")
        if isinstance(value, (s.Basic, s.MatrixBase)):
            assert not value.has(s.Float)
        elif isinstance(value, dict):
            for item in value.values():
                exact(item)
        elif isinstance(value, (tuple, list)):
            for item in value:
                exact(item)

    for packet in audit.packets().values():
        exact(packet)
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 147
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    assert source.data()["checks"] is not source.previous.data()["checks"]
    assert len(source.previous.data()["checks"]) == 100
