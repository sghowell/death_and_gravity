"""Independent threshold, full-parameter, Dirac and cut checks."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_two_loop_elastic_cut import (
    audit,
    bounds,
    calibration,
    cut,
    fermion,
    scalar,
)


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_exact_identity(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[v[0] for v in audit.bad_cases()]
)
def test_invalid_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def bubble_real(v):
    if v == 0:
        return mp.mpf(0)
    if v == 4:
        return mp.mpf(2)
    if v < 0:
        beta = mp.sqrt(1 - 4 / v)
        return 2 - beta * mp.log((beta + 1) / (beta - 1))
    if v < 4:
        r = mp.sqrt(4 / v - 1)
        return 2 - 2 * r * mp.atan(1 / r)
    beta = mp.sqrt(1 - 4 / v)
    return 2 - beta * mp.log((1 + beta) / (1 - beta))


@pytest.mark.parametrize("v", (-2, -1, 0, 1, 3, 4, 5, 6))
def test_logarithmic_threshold_integral_and_bound(v):
    with mp.workdps(45):
        v = mp.mpf(v)
        splits = [mp.mpf(0), mp.mpf(1)]
        if v >= 4:
            beta = mp.sqrt(1 - 4 / v)
            r1, r2 = (1 - beta) / 2, (1 + beta) / 2
            splits = sorted({mp.mpf(0), r1, r2, mp.mpf(1)})

            def delta(x):
                return v * (x - r1) * (x - r2)
        else:

            def delta(x):
                return 1 - v * x * (1 - x)

        def integrand(x):
            d = delta(x)
            return -mp.log(abs(d)) if d else mp.mpf(0)

        val = mp.quad(integrand, splits)
        assert abs(val - bubble_real(v)) < mp.mpf("1e-35")
        assert abs(val) <= 2 + mp.mpf("1e-35")


def box_radial(delta, M, z, y):
    C = M - z * y * (1 - y)
    d = mp.sqrt(M * M - 4 * C * delta)
    log = mp.log(abs((M + d) / (M - d)))
    if delta < 0:
        log += 1j * mp.pi
    return (M * log - 2 * d) / d**3


@cache
def pole_subtracted_box_functions():
    H, r, a, b = s.symbols("H r a b")
    d = b + 2 * a * H
    k = b + a * (H + r)
    f = H * (1 - H)
    numerator = (1 - 2 * H) / d - 2 * a * f / d**2
    u = numerator / k
    u0 = s.factor(u.subs(H, r))
    quotient = s.factor(s.cancel((u - u0) / (H - r)))
    return s.lambdify((r, a, b), u0, "mpmath"), s.lambdify(
        (H, r, a, b), quotient, "mpmath"
    )


def physical_box_subtracted(delta, M, z, y):
    a = delta - z * y * (1 - y)
    b = M - 2 * delta
    root = -2 * delta / (b + mp.sqrt(b * b - 4 * a * delta))
    u0fun, qfun = pole_subtracted_box_functions()
    u0 = u0fun(root, a, b)
    value = u0 * mp.log(abs((1 - root) / root))
    if 0 < root < 1:
        value += 1j * mp.pi * u0
    value += mp.quad(lambda H: qfun(H, root, a, b), [0, 1])
    return value


def physical_triangle_subtracted(delta, M):
    a, b = delta, M - 2 * delta
    d = mp.sqrt(b * b - 4 * a * delta)
    root = -2 * delta / (b + d)
    u0 = (1 - root) / d
    value = u0 * mp.log(abs((1 - root) / root))
    if 0 < root < 1:
        value += 1j * mp.pi * u0
    remainder = lambda H: -(b + a * (1 + root)) / (d * (b + a * (H + root)))
    return value + mp.quad(remainder, [0, 1])


@pytest.mark.parametrize("M", (32, 100))
@pytest.mark.parametrize("delta", ("-.5", "-.1", ".01", ".5", "1.5"))
def test_physical_box_boundary_value_matches_independent_radial_form(M, delta):
    with mp.workdps(46):
        M, delta = mp.mpf(M), mp.mpf(delta)
        for z in (-2, 0, 6):
            y = mp.mpf(".3")
            actual = physical_box_subtracted(delta, M, z, y)
            radial = box_radial(delta, M, z, y)
            assert abs(actual - radial) < mp.mpf("1e-35")
            allowance = (
                6 * (mp.log(M) + 1 + abs(mp.log(abs(delta))) + mp.pi) + 10
            ) / M**2
            assert abs(actual) < allowance
            if delta < 0:
                assert actual.imag > 0


@pytest.mark.parametrize("M", (32, 100))
@pytest.mark.parametrize("delta", (".01", ".5", "1.5"))
def test_direct_positive_parameter_triangle_and_box(M, delta):
    with mp.workdps(42):
        M, delta = mp.mpf(M), mp.mpf(delta)
        triangle = mp.quad(lambda H: (1 - H) / (delta * (1 - H) ** 2 + M * H), [0, 1])
        assert abs(triangle - physical_triangle_subtracted(delta, M)) < mp.mpf("1e-32")
        for z in (-2, 0, 6):
            y = mp.mpf(".3")
            box = mp.quad(
                lambda H, z=z, y=y: (
                    H
                    * (1 - H)
                    / (delta * (1 - H) ** 2 + M * H - z * H * H * y * (1 - y)) ** 2
                ),
                [0, 1],
            )
            assert abs(box - box_radial(delta, M, z, y)) < mp.mpf("1e-32")


@pytest.mark.parametrize("M", (32, 100))
@pytest.mark.parametrize("delta", ("-.5", "-.1", ".01", ".5", "1.5"))
def test_triangle_boundary_value_and_full_uniform_pointwise_bound(M, delta):
    with mp.workdps(42):
        M, delta = mp.mpf(M), mp.mpf(delta)
        value = physical_triangle_subtracted(delta, M)
        bound = (3 * (mp.log(M) + abs(mp.log(abs(delta))) + mp.pi) + 4) / M
        assert abs(value) < bound
        if delta < 0:
            d = mp.sqrt(M * M - 4 * M * delta)
            root = -2 * delta / (M - 2 * delta + d)
            assert abs(value.imag - mp.pi * (1 - root) / d) < mp.mpf("1e-36")


def test_unsplit_absolute_double_pole_is_not_the_boundary_value():
    with mp.workdps(35):
        delta, M, z, y = mp.mpf("-.3"), mp.mpf(32), 0, mp.mpf(".5")
        a, b = delta, M - 2 * delta
        root = -2 * delta / (b + mp.sqrt(b * b - 4 * a * delta))
        epsilon1, epsilon2 = mp.mpf("1e-4"), mp.mpf("1e-6")
        f = lambda H: H * (1 - H) / (delta * (1 - H) ** 2 + M * H) ** 2
        first = mp.quad(f, [0, root - epsilon1]) + mp.quad(f, [root + epsilon1, 1])
        second = mp.quad(f, [0, root - epsilon2]) + mp.quad(f, [root + epsilon2, 1])
        assert second > 50 * first
        assert abs(box_radial(delta, M, z, y)) < 1


def gamma_matrices():
    zero = mp.zeros(2)
    identity = mp.eye(2)
    pauli = [
        mp.matrix([[0, 1], [1, 0]]),
        mp.matrix([[0, -1j], [1j, 0]]),
        mp.matrix([[1, 0], [0, -1]]),
    ]

    def block(a, b, c, d):
        out = mp.zeros(4)
        for i in range(2):
            for j in range(2):
                out[i, j], out[i, j + 2] = a[i, j], b[i, j]
                out[i + 2, j], out[i + 2, j + 2] = c[i, j], d[i, j]
        return out

    return [block(identity, zero, zero, -identity)] + [
        block(zero, P, P, zero) for P in pauli
    ]


def opnorm(matrix):
    return max(mp.svd(matrix, compute_uv=False))


@pytest.mark.parametrize("mass", (24, 37))
@pytest.mark.parametrize("angle", ("-1", ".3", "1"))
def test_literal_Dirac_resolvent_product_difference(mass, angle):
    with mp.workdps(34):
        gammas = gamma_matrices()
        eye = mp.eye(4)
        for i in range(4):
            for j in range(4):
                assert (
                    mp.norm(
                        gammas[i] * gammas[j]
                        + gammas[j] * gammas[i]
                        - (2 * eye if i == j else mp.zeros(4))
                    )
                    == 0
                )
        E, p, z = mp.sqrt(mp.mpf(3) / 2), mp.sqrt(mp.mpf(1) / 2), mp.mpf(angle)
        p1 = [1j * E, 0, 0, p]
        p2 = [1j * E, 0, 0, -p]
        p4 = [-1j * E, p * mp.sqrt(1 - z * z), 0, p * z]
        routes = [[0] * 4, p1, [a + b for a, b in zip(p1, p2)], [-v for v in p4]]
        q = list(map(mp.mpf, (".3", ".2", ".1", ".4")))
        m = mp.mpf(mass)
        D0 = m * eye
        for i in range(4):
            D0 += 1j * gammas[i] * q[i]
        S0 = D0**-1
        n = mp.sqrt(m * m + sum(v * v for v in q))
        assert abs(opnorm(S0) - 1 / n) < mp.mpf("1e-30")
        product = eye
        for r in routes:
            R = mp.zeros(4)
            for i in range(4):
                R += 1j * gammas[i] * r[i]
            assert opnorm(R) < 6
            product = product * (D0 + R) ** -1
        difference = opnorm(product - S0**4)
        majorant = 24 / (n**5 * (1 - 6 / m) ** 4)
        assert difference < majorant


@pytest.mark.parametrize("m", (24, 100))
def test_finite_radial_norm_integral(m):
    with mp.workdps(38):
        m = mp.mpf(m)
        radial = mp.quad(
            lambda r: r**3 / (r * r + m * m) ** mp.mpf("2.5"), [0, m, mp.inf]
        )
        assert abs(radial - 2 / (3 * m)) < mp.mpf("1e-34")
        assert (1 - 6 / m) ** -4 < 4


def test_full_local_fermion_vertex_is_required():
    x = s.Symbol("x")
    series = s.series(
        (1 + x) ** 4 * (2 * s.log(1 + x) - s.Rational(3, 2)), x, 0, 5
    ).removeO()
    assert series.coeff(x, 4) == s.Rational(8, 3)
    assert 24 * series.coeff(x, 4) == 64
    assert fermion.data()["zero_momentum_MS_scattering_vertex"] != 0


@pytest.mark.parametrize("lam,B", (("1", ".2"), (".03", "2")))
def test_independent_optical_interference_integral_and_factor(lam, B):
    with mp.workdps(38):
        lam, B = mp.mpf(lam), mp.mpf(B)

        def density(v):
            # Explicit bounded test amplitudes; no assertion that these
            # polynomials equal the GY14 amplitude.
            angular = mp.quad(
                lambda z: lam * (40 + 10 * z * z + 5 * (v - 4)) * B * (2 * z * z - 1),
                [0, 1],
            )
            return mp.sqrt(1 - 4 / v) * angular / (16 * mp.pi)

        result = 2 / mp.pi * mp.quad(lambda v: density(v) / (v - 2) ** 3, [4, 5, 6])
        assert abs(result) < 73 * lam * B / 1280
        assert result != 0
        wrong = result / 2
        assert abs(result - wrong) > mp.mpf("1e-6") * lam * B


def test_exact_parameter_bounds_and_scope():
    d = calibration.data()
    assert all(d["bounds"].values())
    assert d["complete_through_two_loop_improved_formal_band"][
        "positive_formal_uniform_lower"
    ]
    assert len(cut.rows()) == 6
    assert len(audit.matching()) == 17
    assert sum(v["status"] == "OPEN" for v in audit.matching()) == 5
    assert scalar.data()["checks"]["box_weight_zero_at_first_endpoint"] == 0
    assert bounds.improved_band(1, 2, 1)["positive_formal_uniform_lower"] is False
    assert audit.controls()["original_P8_not_closed"] is True
