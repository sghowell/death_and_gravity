"""Full-dimensional and time-domain independent quotient-reference checks."""

import mpmath as mp
import pytest
import sympy as s
from p8_affine import verify as serializer
from p8_vacuum_affine_flat_scalar_quotient_inverse import audit, estimates
from p8_vacuum_affine_flat_scalar_quotient_inverse import geometry as g
from p8_vacuum_affine_flat_scalar_quotient_inverse import resolvent as r


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_integrated_exact_residual(name, value):
    rows = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(s.factor(v) == 0 for v in rows), name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_integrated_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_and_false_completion_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_scoped_stages_and_previous_frontier_unchanged():
    for stage in (
        "flat_gauge_quotient",
        "two_channel_reference_inverse",
        "uniform_spatial_bound",
    ):
        assert audit.require_stage(stage) == stage
    for t in (-s.Rational(1, 2), 0, s.Rational(1, 2)):
        assert audit.require_scope(t)[0] == t
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.validate_scope(audit.frontier(), audit.matching())


def test_supported_exact_momenta_and_continuous_endpoint():
    for value in (0, s.Rational(1, 10**12), 1, 10000):
        assert r.transfer(value) == value
        assert r.coordinate_kernel(value).shape == (2, 2)


@pytest.mark.parametrize(
    "direction",
    (
        (1, 0, 0),
        (s.Rational(3, 5), s.Rational(4, 5), 0),
        (s.Rational(2, 3), s.Rational(-2, 3), s.Rational(1, 3)),
    ),
)
def test_independent_full_four_metric_projectors_and_curvature(direction):
    lam, k = s.symbols("independent_lambda transfer", positive=True)
    q = k * k
    p = lam * lam + q
    P = k * s.Matrix(direction)
    assert P.dot(P) == q
    nd, zd, bd, ng, zg, bg = s.symbols("nD zD bD nG zG bG")
    eta = s.diag(1, -1, -1, -1)
    derivative = s.Matrix([lam, *[s.I * x for x in P]])
    T = p * eta - derivative * derivative.T

    def metric(n, z, b, sign):
        out = s.diag(2 * n, -2 * z, -2 * z, -2 * z)
        for j in range(3):
            out[0, j + 1] = out[j + 1, 0] = sign * s.I * P[j] * b / q
        return out

    HD, HG = metric(nd, zd, bd, -1), metric(ng, zg, bg, 1)
    UD, UG = eta * HD * eta, eta * HG * eta
    spin0 = s.trace(T * UD) * s.trace(T * UG) / 3
    spin2 = s.trace(T * UD * T * UG) - spin0
    SD = q * nd / 3 + (lam * lam + 2 * q / 3) * zd + lam * bd / 3
    SG = q * ng / 3 + (lam * lam + 2 * q / 3) * zg - lam * bg / 3
    WD = q * (nd - zd) + lam * bd
    WG = q * (ng - zg) - lam * bg
    assert s.expand(spin0 - 12 * SD * SG) == 0
    assert s.expand(spin2 - s.Rational(8, 3) * WD * WG) == 0

    def curvature(H):
        return {
            (a, b, c, d): s.expand(
                (
                    derivative[b] * derivative[c] * H[a, d]
                    + derivative[a] * derivative[d] * H[b, c]
                    - derivative[b] * derivative[d] * H[a, c]
                    - derivative[a] * derivative[c] * H[b, d]
                )
                / 2
            )
            for a in range(4)
            for b in range(4)
            for c in range(4)
            for d in range(4)
        }

    RD, RG = curvature(HD), curvature(HG)

    def ricci(curv):
        return s.Matrix(
            4, 4, lambda b, d: sum(eta[a, a] * curv[a, b, a, d] for a in range(4))
        )

    RicD, RicG = ricci(RD), ricci(RG)
    scalarD, scalarG = s.trace(eta * RicD), s.trace(eta * RicG)
    riem = sum(
        eta[a, a] * eta[b, b] * eta[c, c] * eta[d, d] * RD[a, b, c, d] * RG[a, b, c, d]
        for a in range(4)
        for b in range(4)
        for c in range(4)
        for d in range(4)
    )
    ric = s.trace(eta * RicD * eta * RicG)
    rs = scalarD * scalarG
    Euler = s.expand(riem - 4 * ric + rs)
    Weyl = s.expand(riem - 2 * ric + rs / 3)
    assert Euler == 0
    assert s.expand(rs - 36 * SD * SG) == 0
    assert s.expand(Weyl - s.Rational(4, 3) * WD * WG) == 0
    literal = -Weyl / 15 + Euler / 45 - rs / 9
    assert s.expand(literal + 4 * SD * SG + s.Rational(4, 45) * WD * WG) == 0


@pytest.mark.parametrize("q", (s.Rational(1, 10**20), 1, 10000))
@pytest.mark.parametrize("lam", (s.Rational(1, 3), 2))
def test_independent_forced_four_metric_gauge_quotient_and_inverse(q, lam):
    q, lam = map(s.Rational, (q, lam))
    eta, z, c = s.symbols("independent_eta independent_zeta independent_c")
    # Both time legs are separately substituted, not assumed to have identical maps.
    for sign in (-1, 1):
        n = sign * lam * eta
        b = sign * lam * c + q * eta
        S = q * n / 3 + (lam * lam + 2 * q / 3) * z - sign * lam * b / 3
        W = q * (n - z) - sign * lam * b
        actual = s.Matrix([S, W])
        B = g.B.subs({g.q: q, g.lam: lam})
        assert (actual - B * s.Matrix([z, c])).applyfunc(s.expand) == s.zeros(2, 1)
        wrong = s.Matrix(
            [
                q * n / 3 + (lam * lam + 2 * q / 3) * z + sign * lam * b / 3,
                q * (n - z) + sign * lam * b,
            ]
        )
        assert wrong != actual
    f0 = -s.Rational(19, 4)
    f2 = -s.Rational(7, 30)
    B = g.B.subs({g.q: q, g.lam: lam})
    matrix = B.T * s.diag(f0, s.Rational(8, 3) * f2) * B
    force = s.Matrix([s.Rational(7, 11), s.Rational(-13, 17)])
    direct = matrix.inv() * force
    inverse = r.BINV.subs({g.q: q, g.lam: lam})
    predicted = inverse * s.diag(1 / f0, s.Rational(3, 8) / f2) * inverse.T * force
    assert direct == predicted
    deleted = inverse * s.diag(1 / f0, 0) * inverse.T * force
    assert direct != deleted


@pytest.mark.parametrize("qval", (0, s.Rational(1, 10**12), 1, 10000))
def test_independent_causal_differential_coordinate_inverse(qval):
    qval = s.Rational(qval)
    t = s.Symbol("test_time", real=True)
    f = s.Matrix([t**4 / 12, t**5 / 60])
    # Solve the two forced second-order equations directly.
    drive = f[0] - f[1] / 3
    z = s.dsolve(
        s.diff(s.Function("z")(t), t, 2) + qval * s.Function("z")(t) - drive,
        ics={s.Function("z")(0): 0, s.Subs(s.diff(s.Function("z")(t), t), t, 0): 0},
    ).rhs
    c = s.integrate(s.integrate(-qval * z - f[1], t), t)
    c = c - c.subs(t, 0) - t * s.diff(c, t).subs(t, 0)
    observed = s.Matrix(
        [
            s.diff(z, t, 2) + 2 * qval * z / 3 - s.diff(c, t, 2) / 3,
            -qval * z - s.diff(c, t, 2),
        ]
    )
    assert (observed - f).applyfunc(s.simplify) == s.zeros(2, 1)
    # Explicit B^-1 cancellation avoids division byq in the continuous endpoint.
    x = s.Symbol("integration_time", real=True)
    primitive = s.integrate((t - x) * (f[0] + 2 * f[1] / 3).subs(t, x), (x, 0, t))
    assert s.simplify(c - z + primitive) == 0


def test_independent_coordinate_bounds_and_initial_primitive():
    c = s.Symbol("cosine", real=True)
    M = s.Matrix([[c, -c / 3], [c - 1, -c / 3 - s.Rational(2, 3)]])
    norm = s.expand(sum(x * x for x in M))
    assert norm == (20 * c * c - 14 * c + 13) / 9
    # Convex quadratic attains its maximum at an endpoint on[-1,1].
    assert s.diff(norm, c, 2) > 0
    assert (
        max(norm.subs(c, -1), norm.subs(c, 1)) == s.Rational(47, 9) < s.Rational(25, 4)
    )
    assert r.coordinate_kernel(0) == s.Matrix([[r.t, -r.t / 3], [0, -r.t]])
    G = s.Matrix([[12, -4], [-4, 4]])
    assert (G - 2 * s.eye(2)).is_positive_definite
    assert (14 * s.eye(2) - G).is_positive_definite


def test_independent_shifted_full_channel_dispersion():
    # Both actual continuum densities are entered independently. The checks
    # validate shifts and normalization, not continuum existence by quadrature.
    with mp.workdps(60):

        def A2(p):
            return mp.mpf(1) / 30 + mp.quad(
                lambda y: (
                    p
                    * y
                    * y
                    * (30 - 20 * y * y + 3 * y**4)
                    / (30 * (4 + p * (1 - y * y)))
                ),
                [0, mp.mpf(".5"), mp.mpf(".9"), 1],
            )

        def A0(p):
            return 4 + mp.quad(
                lambda y: (
                    p * y * y * (3 - 2 * y * y + 3 * y**4) / (4 + p * (1 - y * y))
                ),
                [0, mp.mpf(".5"), mp.mpf(".9"), 1],
            )

        root = -mp.findroot(A2, (-mp.mpf(".58"), -mp.mpf(".5798")))
        R = 1 / mp.diff(A2, -root)

        def density(u, spin):
            z = -mp.expm1(-u)
            b = mp.sqrt(z)
            atanh = u / 2 + mp.log1p(b)
            if spin == 2:
                poly = 30 - 20 * z + 3 * z * z
                D = (
                    -mp.mpf(172) / 225
                    + 19 * z / 30
                    - z * z / 10
                    + b * poly * atanh / 30
                )
                U = b * poly / 60
            else:
                poly = 3 - 2 * z + 3 * z * z
                D = mp.mpf(16) / 15 + z - 3 * z * z + b * poly * atanh
                U = b * poly / 2
            return U / (D * D + mp.pi**2 * U * U)

        regions = [0, mp.mpf(".1"), mp.mpf(".5"), 1, 2, 4, 10, 30, mp.inf]
        for q in (mp.mpf(0), mp.mpf(1), mp.mpf(10000)):
            for lam in (mp.mpf(1), mp.mpc(1, 2)):
                p = lam * lam + q
                cut0 = mp.quad(
                    lambda u, point=p: density(u, 0) / (1 + point * mp.exp(-u) / 4),
                    regions,
                )
                cut2 = mp.quad(
                    lambda u, point=p: density(u, 2) / (1 + point * mp.exp(-u) / 4),
                    regions,
                )
                assert abs(cut0 - 1 / A0(p)) < mp.mpf("1e-45")
                assert abs(R / (p + root) + cut2 - 1 / A2(p)) < mp.mpf("1e-45")
                B = mp.matrix(
                    [[lam * lam + 2 * q / 3, -lam * lam / 3], [-q, -lam * lam]]
                )
                coordinate = mp.matrix(
                    [
                        [1 / p, -1 / (3 * p)],
                        [1 / p - 1 / (lam * lam), -1 / (3 * p) - 2 / (3 * lam * lam)],
                    ]
                )
                D = mp.diag([-A0(p), -8 * A2(p) / 3])
                cut_inverse = mp.diag([-cut0, -3 * (R / (p + root) + cut2) / 8])
                recovered = coordinate * cut_inverse * coordinate.T
                assert mp.norm((B.T * D * B) * recovered - mp.eye(2)) < mp.mpf("1e-38")
            assert 0 < 1 / A0(q) <= mp.mpf(1) / 4
            assert 0 < 1 / A2(q) <= 30


def test_independent_ordered_time_convolution_constants_and_units():
    t, a, b = s.symbols("time a b", positive=True)
    C = s.Rational(5, 2)
    J = s.Rational(45, 2)
    kernel = s.integrate(s.integrate(C * J * C * (t - a - b), (b, 0, t - a)), (a, 0, t))
    operator = s.integrate(kernel, (t, 0, 1))
    assert operator == s.Rational(375, 64) < 6
    kap = s.Integer(10) ** 800
    physical = 64 * s.pi**2 * kap * operator
    assert s.simplify(physical - 375 * s.pi**2 * kap) == 0
    assert physical > 1


@pytest.mark.parametrize("module", (g, r, estimates))
def test_exact_module_checks_and_gates(module):
    data = module.data()
    for key, value in data["checks"].items():
        values = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.factor(v) == 0 for v in values), key
    assert all(bool(v) for v in data["gates"].values())
    serializer.serialize(
        {k: v for k, v in data.items() if k not in ("checks", "gates")}
    )
