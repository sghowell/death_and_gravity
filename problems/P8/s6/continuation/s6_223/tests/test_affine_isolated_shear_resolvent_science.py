"""Independent curvature, cut, enclosure and pole-retention checks."""

import mpmath as mp
import pytest
import sympy as s
from p8_affine import verify as serializer
from p8_vacuum_affine_isolated_shear_resolvent import (
    audit,
    kernel,
    normalization,
    spectral,
)

PACKETS = (normalization.data, spectral.data, kernel.data)


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_integrated_exact_residual(name, value):
    assert s.factor(value) == 0, name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_integrated_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_and_false_completion_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize("packet", PACKETS)
def test_exact_packets_and_serialization(packet):
    data = packet()
    assert all(s.factor(v) == 0 for v in data["checks"].values())
    assert all(bool(v) for v in data["gates"].values())
    serializer.serialize(
        {k: v for k, v in data.items() if k not in ("checks", "gates")}
    )


def test_scoped_stages_and_all_previous_rows_retained():
    for stage in ("isolated_factor", "pole_plus_cut", "finite_window_inverse"):
        assert audit.require_stage(stage) == stage
    for t in (-s.Rational(1, 2), 0, s.Rational(1, 2)):
        assert audit.require_scope(t)[0] == t
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.validate_scope(audit.frontier(), audit.matching())


@pytest.mark.parametrize(
    "entries", ((2, -1, -1, 0, 0, 0), (1, -1, 0, 2, 0, 0), (3, -2, -1, 1, -2, 3))
)
def test_independent_literal_linear_curvature_and_finite_TF_Hessian(entries):
    aa, bb, cc, ab, ac, bc = map(s.Integer, entries)
    Q = s.Matrix([[aa, ab, ac], [ab, bb, bc], [ac, bc, cc]])
    assert s.trace(Q) == 0
    norm = s.trace(Q * Q)
    q = s.Symbol("independent_second_time_jet", real=True)
    h = s.zeros(4)
    h[1:4, 1:4] = -q * Q
    eta = (1, -1, -1, -1)

    def dd(a, b, c, d):
        return h[a, b] if c == d == 0 else 0

    curvature = {}
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    curvature[a, b, c, d] = (
                        dd(a, d, b, c)
                        + dd(b, c, a, d)
                        - dd(a, c, b, d)
                        - dd(b, d, a, c)
                    ) / 2
    Ric = s.Matrix(
        4, 4, lambda b, d: sum(eta[a] * curvature[a, b, a, d] for a in range(4))
    )
    scalar = sum(eta[a] * Ric[a, a] for a in range(4))
    Ric2 = sum(eta[a] * eta[b] * Ric[a, b] ** 2 for a in range(4) for b in range(4))
    Riem2 = sum(
        eta[a] * eta[b] * eta[c] * eta[d] * value**2
        for (a, b, c, d), value in curvature.items()
    )
    Weyl2 = s.expand(Riem2 - 2 * Ric2 + scalar * scalar / 3)
    Euler = s.expand(Riem2 - 4 * Ric2 + scalar * scalar)
    assert scalar == 0 and Euler == 0
    assert s.factor(Weyl2 - q * q * norm / 2) == 0
    literal = -Weyl2 / 30 + Euler / 90 - scalar * scalar / 18
    assert s.factor(s.diff(literal, q, 2) / norm + s.Rational(1, 30)) == 0
    assert s.diff(literal, q, 2) / norm != literal.coeff(q, 2) / norm


@pytest.mark.parametrize("spin,factor", ((2, 64), (0, 768)))
def test_independent_direct_spectral_mass_substitution(spin, factor):
    y, m, p = s.symbols("independent_y independent_mass independent_p", positive=True)
    sig = 4 * m * m / (1 - y * y)
    poly = (
        (13 * sig**2 + 56 * m * m * sig + 48 * m**4) / 3840
        if spin == 2
        else (sig**2 - 4 * m * m * sig + 12 * m**4) / 384
    )
    transformed = s.factor(
        factor * p * y * poly * s.diff(sig, y) / (sig**3 * (sig + p))
    )
    W = normalization.W2 if spin == 2 else normalization.W0
    expected = p * W.subs(normalization.y, y) / (4 * m * m + p * (1 - y * y))
    assert s.factor(transformed - expected) == 0
    if spin == 0:
        assert s.factor(transformed / 12 - expected) != 0


def test_independent_binomial_moments_and_shorter_rigorous_enclosure():
    N = 12
    M = [
        sum(
            (-1) ** k
            * s.binomial(j, k)
            * (
                s.Rational(30, 3 + 2 * k)
                - s.Rational(20, 5 + 2 * k)
                + s.Rational(3, 7 + 2 * k)
            )
            for k in range(j + 1)
        )
        for j in range(N + 1)
    ]
    assert tuple(M) == spectral.moments()[: N + 1]
    lo, hi = s.Rational(2899, 5000), s.Rational(5799, 10000)
    intervals = []
    derivatives = []
    for r in (lo, hi):
        x = r / 4
        partial = s.Rational(1, 30) - r * sum(M[j] * x**j for j in range(N + 1)) / 120
        remainder = r * M[0] * x ** (N + 1) / (120 * (1 - x))
        intervals.append((partial - remainder, partial))
        der = sum((j + 1) * M[j] * x**j for j in range(N + 1)) / 120
        tail = M[0] * x ** (N + 1) * ((N + 2) - (N + 1) * x) / (120 * (1 - x) ** 2)
        derivatives.append((der, der + tail))
    assert intervals[0][0] > 0 > intervals[1][1]
    assert 16 < 1 / derivatives[1][1] < 1 / derivatives[0][0] < 17


def test_independent_complete_pole_plus_cut_high_precision_reconstruction():
    # Direct radial and cut integrals, not a lambdified production formula.
    # Numerical comparisons do not prove zero counts or continuum integrability.
    with mp.workdps(65):

        def radial(p):
            return mp.mpf(1) / 30 + mp.quad(
                lambda y: (
                    p
                    * y
                    * y
                    * (30 - 20 * y * y + 3 * y**4)
                    / (30 * (4 + p * (1 - y * y)))
                ),
                [0, mp.mpf(".4"), mp.mpf(".7"), mp.mpf(".9"), 1],
            )

        def closed(p):
            d = 1 + 4 / p
            return (
                -mp.mpf(172) / 225
                + 19 * d / 30
                - d * d / 10
                + mp.sqrt(d) * (30 - 20 * d + 3 * d * d) * mp.atanh(1 / mp.sqrt(d)) / 30
            )

        r = -mp.findroot(radial, (-mp.mpf(".58"), -mp.mpf(".5798")))
        R = 1 / mp.diff(radial, -r)

        def rho(u):
            z = -mp.expm1(-u)
            beta = mp.sqrt(z)
            P = 30 - 20 * z + 3 * z * z
            D = (
                -mp.mpf(172) / 225
                + 19 * z / 30
                - z * z / 10
                + beta * P * (u / 2 + mp.log1p(beta)) / 30
            )
            U = beta * P / 60
            return U / (D * D + mp.pi**2 * U * U)

        regions = [0, mp.mpf(".1"), mp.mpf(".5"), 1, 2, 4, 10, 30, mp.inf]
        for p in (
            mp.mpf(1),
            mp.mpc(1, 2),
            mp.mpc(-3, mp.mpf(".4")),
            mp.mpc(-9, mp.mpf(".3")),
        ):
            assert abs(radial(p) - closed(p)) < mp.mpf("1e-55")
            if mp.im(p) > 0:
                assert mp.im(radial(p)) > 0
        for p in (0, mp.mpf(1), mp.mpc(1, 2), -mp.mpf(".2")):
            cut = mp.quad(
                lambda u, point=p: rho(u) / (1 + (point / 4) * mp.exp(-u)), regions
            )
            assert abs(R / (p + r) + cut - 1 / radial(p)) < mp.mpf("1e-50")
            assert abs(cut - 1 / radial(p)) > 1
            assert abs(-R / (p + r) + cut - 1 / radial(p)) > 1
        static = mp.quad(rho, regions)
        nextmoment = mp.quad(lambda u: rho(u) * mp.exp(-u) / 4, regions)
        assert abs(R / r + static - 30) < mp.mpf("1e-50")
        assert abs(R / r**2 + nextmoment - mp.mpf(675) / 14) < mp.mpf("1e-50")
        assert 0 < static < 3 and 0 < nextmoment < mp.mpf(".1")


def test_independent_bounded_forward_kernel_beta_integrals():
    x = s.Symbol("angle", real=True)
    integral = s.integrate(
        s.sin(x) ** 2 - s.Rational(2, 3) * s.sin(x) ** 4 + s.sin(x) ** 6 / 10,
        (x, 0, s.pi / 2),
    )
    assert integral == 9 * s.pi / 64


@pytest.mark.parametrize("mass", (1, 2, 1000))
def test_independent_mass_scaling_and_prepared_pole_convolution(mass):
    mass = s.Integer(mass)
    t, lam = s.symbols("time lambda", positive=True)
    r = s.Rational(3, 5)
    weight = s.Rational(33, 2) * mass**2
    omega = mass * s.sqrt(r)
    # Rational illustrative frequency/weight fixture, not the actual enclosed root.
    pole = -weight * s.sin(omega * t) / omega
    K1 = -s.Rational(33, 2) * s.sin(s.sqrt(r) * t) / s.sqrt(r)
    assert s.simplify(pole - mass * K1.subs(t, mass * t)) == 0
    primitive = -weight * (1 - s.cos(omega * t)) / omega**2
    # For f(t)=t, J*f' has this independently integrated closed expression.
    solution = -weight / omega**2 * (t - s.sin(omega * t) / omega)
    assert s.simplify(s.diff(solution, t) - primitive) == 0
    assert solution.subs(t, 0) == 0
    assert (
        s.simplify(
            s.laplace_transform(solution, t, lam, noconds=True)
            + weight / (lam**2 * (lam**2 + omega**2))
        )
        == 0
    )
