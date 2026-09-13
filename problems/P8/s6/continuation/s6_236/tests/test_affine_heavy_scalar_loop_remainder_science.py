"""Independent primitive, branch, Cauchy, full diagram and physical-window tests."""

import itertools

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_heavy_scalar_four_point_loop import audit as previous
from p8_vacuum_affine_heavy_scalar_loop_remainder import audit, bounds
from p8_vacuum_affine_heavy_scalar_tree_matching import model


def mpq(v):
    q = s.Rational(v)
    return mp.mpf(int(q.p)) / int(q.q)


def physical_log(L):
    if L < 0:
        return mp.log(-L) - mp.j * mp.pi
    return mp.log(L)


def closed_T_D(n, L, b):
    d = mp.sqrt(n * n - 4 * n * L + 4 * L * b)
    a = (n - 2 * L + d) / 2
    beta = (L - b) / a
    v = L / a
    H = (1 + beta) * mp.log1p(beta) / beta if beta else mp.mpf(1)
    Hp = (beta - mp.log1p(beta)) / (beta * beta) if beta else mp.mpf("0.5")
    log = mp.log(a + L) - physical_log(L)
    R = (1 + v) * log - H
    dn = (n - 2 * L) / d
    an = (1 + dn) / 2
    T = R / d
    D = R * dn / (d * d) - (an / a) * (1 - v * log + beta * Hp) / d
    return T, D


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_every_exact_entry(name, value):
    assert s.cancel(value) == 0, name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_every_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize(
    "n,L,b",
    [
        (100, s.Rational(2, 3), 0),
        (100, s.Rational(2, 3), s.Rational(2, 3)),
        (100, s.Rational(2, 3), -1),
        (100, 3, 2),
        (1000, 7, -2),
        (1000, s.Rational(1, 100), 3),
    ],
)
def test_entire_positive_parameter_integrals_and_independent_mass_derivative(n, L, b):
    with mp.workdps(70):
        nn, ll, bb = map(mpq, (n, L, b))
        den = lambda z: nn * z + (1 - z) ** 2 * ll - bb * z * z
        ts = mp.quad(lambda z: (1 - z) / den(z), [0, 1 / nn, mp.mpf("0.1"), 1])
        ds = mp.quad(lambda z: z * (1 - z) / den(z) ** 2, [0, 1 / nn, mp.mpf("0.1"), 1])
        T, D = closed_T_D(nn, ll, bb)
        assert abs(T - ts) < mp.mpf("1e-60")
        assert abs(D - ds) < mp.mpf("1e-60")
        numerical_derivative = -mp.diff(lambda n0: closed_T_D(n0, ll, bb)[0], nn)
        assert abs(numerical_derivative - D) < mp.mpf("1e-60")


@pytest.mark.parametrize("L,b", [(-s.Rational(1, 10), 0), (-1, -2), (-3, 2)])
def test_physical_light_cut_against_independent_delta_function_and_sign(L, b):
    with mp.workdps(80):
        n = mp.mpf(1000)
        ll, bb = map(mpq, (L, b))
        delta = mp.sqrt(n * n - 4 * n * ll + 4 * ll * bb)
        alpha = (n - 2 * ll + delta) / 2
        root = -ll / alpha
        jac = n - 2 * ll + 2 * (ll - bb) * root
        independent = mp.pi * (1 - root) / abs(jac)
        T, D = closed_T_D(n, ll, bb)
        assert 0 < root < 1
        assert abs(mp.im(T) - independent) < mp.mpf("1e-70")
        assert mp.im(T) > 0 and mp.im(D) > 0
        wrong = -mp.pi * (1 + ll / alpha) / delta
        assert abs(mp.im(T) - wrong) > mp.mpf("1e-4")


@pytest.mark.parametrize("lratio,bratio", list(itertools.product((-1, 0, 1), repeat=2)))
def test_complete_complex_coefficient_disk_independent_samples(lratio, bratio):
    with mp.workdps(70):
        K = mp.mpf(7)
        L, b = K * lratio, K * bratio
        for j in range(16):
            w = mp.exp(2 * mp.pi * mp.j * j / 16) / (16 * K)
            d = mp.sqrt(1 - 4 * L * w + 4 * L * b * w * w)
            a = (1 - 2 * L * w + d) / 2
            beta = (L - b) * w / a
            v = L * w / a
            h = (1 + beta) * mp.log1p(beta) / beta if beta else 1
            P = (1 + v) / d
            Q = ((1 + v) * mp.log((1 + d) / 2) - h) / d
            assert abs(d - 1) < mp.mpf(1) / 4
            assert abs(a - 1) < mp.mpf(3) / 16
            assert abs(P) < 2 and abs(Q) < 4


@pytest.mark.parametrize("n", [1000, 10**6, model.MASS2])
@pytest.mark.parametrize(
    "L,b", [(s.Rational(2, 3), s.Rational(1, 3)), (-1, -2), (3, 2), (-3, 0)]
)
def test_complete_T_and_D_remainders_including_actual_hierarchy(n, L, b):
    with mp.workdps(900 if n == model.MASS2 else 80):
        nn, ll, bb = map(mpq, (n, L, b))
        S = mp.mpf(16)
        w = 1 / nn
        H = mp.log(nn) - physical_log(ll)
        T, D = closed_T_D(nn, ll, bb)
        tj = w * (H - 1) + w * w * (3 * ll * H + (bb - 7 * ll) / 2)
        dj = w * w * (H - 2) + w**3 * (6 * ll * H - 10 * ll + bb)
        assert abs(T - tj) <= 64 * S * S * w**3 * (abs(H) + 2)
        assert abs(D - dj) <= 256 * S * S * w**4 * (abs(H) + 3)
        assert abs(T - tj) > 0 and abs(D - dj) > 0


@pytest.mark.parametrize(
    "r", [0, s.Rational(1, 10), s.Rational(1, 2), s.Rational(9, 10), 1]
)
def test_integrated_root_logarithm_bound_including_coincident_threshold_root(r):
    with mp.workdps(60):
        rr = mpq(r)
        endpoints = sorted({mp.mpf(0), rr, mp.mpf(1)})
        actual = mp.quad(lambda x: -mp.log(abs(x - rr)), endpoints)
        term = lambda x: 0 if x == 0 else x * mp.log(x)
        expected = 1 - term(rr) - term(1 - rr)
        assert abs(actual - expected) < mp.mpf("1e-50")
        assert actual < 2


@pytest.mark.parametrize(
    "S,cosine",
    [
        (4, -1),
        (4, 1),
        (5, s.Rational(1, 7)),
        (16, -1),
        (16, 0),
        (16, 1),
        (bounds.S_MAX, -1),
        (bounds.S_MAX, 0),
        (bounds.S_MAX, 1),
    ],
)
def test_full_physical_parameter_cubes_have_the_claimed_scale(S, cosine):
    name, svalue, tvalue, uvalue, order = audit.require_physical(
        "V2S-T1-OS4", S, cosine, 1
    )
    assert name == "V2S-T1-OS4" and order == 1
    assert svalue + tvalue + uvalue == 4
    K = s.Rational(S, 4)
    for a in (svalue, tvalue, uvalue):
        for x in (0, s.Rational(1, 7), s.Rational(1, 2), 1):
            assert abs(1 - a * x * (1 - x)) <= K
            assert abs(a * x * (1 - x)) <= K
    assert S / model.MASS2 < s.Rational(1, 8)


def test_independent_all_diagram_cancellation_with_distinct_light_functions():
    w, N = s.symbols("w log_n")
    x, y = s.symbols("s t")
    channels = (x, y, 4 - x - y)
    Bs = s.symbols("B0:3")
    Js = s.symbols("J0:3")
    pieces = []
    for i, a in enumerate(channels):
        others = [channels[j] for j in range(3) if j != i]
        A = -2 * w + (a - 4) * w * w
        C = w * (N - 1 + Bs[i]) + w * w * (
            3 * (1 - a / 6) * N - 3 * Js[i] - s.Rational(7, 2) * (1 - a / 6)
        )
        boxes = sum(
            w * w * (N - 2 + Bs[i])
            + w**3 * (6 * (1 - a / 6) * N - 6 * Js[i] - 10 * (1 - a / 6) + b / 6)
            for b in others
        )
        pieces.append(A * A * Bs[i] / 2 + 2 * A * C + boxes)
    full = s.series(sum(pieces), w, 0, 4).removeO().expand()
    assert s.expand(full - (-6 * N * w * w + (-16 * N + s.Rational(10, 3)) * w**3)) == 0
    for function in Bs + Js:
        assert s.diff(full, function) == 0
    # Omitting any channel bubble fails the nonanalytic cancellation.
    for i, a in enumerate(channels):
        wrong = s.series(
            sum(pieces) - (-2 * w + (a - 4) * w * w) ** 2 * Bs[i] / 2, w, 0, 4
        ).removeO()
        assert s.expand(s.diff(wrong, Bs[i])) != 0


def test_full_mass_derivative_keeps_the_logarithm_derivative():
    w, H, L, b = s.symbols("w H L b")
    T = w * (H - 1) + w * w * (3 * L * H + (b - 7 * L) / 2)
    correct = s.expand(w * w * (s.diff(T, w) - s.diff(T, H) / w))
    wrong = s.expand(w * w * s.diff(T, w))
    assert s.expand(wrong - correct - w * w - 3 * L * w**3) == 0
    assert correct != wrong


def test_independent_analytic_tail_sums_and_entire_constant_budget():
    x = s.Symbol("x")
    remainder = x * x / (1 - x)
    weighted = s.cancel(remainder + x * s.diff(remainder, x))
    assert s.cancel(weighted - x * x * (3 - 2 * x) / (1 - x) ** 2) == 0
    assert weighted.subs(x, s.Rational(1, 2)) == 2
    assert 2 + 12 + s.Rational(8, 64) < 16
    assert 384 + 12 + 14 < 512
    assert 16 + 512 + 2 * 256 < 2048
    assert 3 * 2048 * 1000 < 10**7
    assert 2 * 462 + 8 + 3 < 1000


def test_actual_relative_bound_and_scope_not_exact_quantum_amplitude():
    n = model.MASS2
    g2 = model.G2
    assert model.LAMBDA == g2 / (2 * (n - 2) ** 3)
    upper = s.Rational(10**7, 9) * g2 / n
    assert upper < s.Rational(1, 10**199)
    assert s.Rational(1, 60) + upper < s.Rational(1, 59)
    assert audit.frontier() == previous.frontier()
    assert audit.matching()[:-1] == previous.matching()
    assert len(audit.matching()) == 92
    assert "truncation" in audit.observable()["specified_truncation"]
    assert "omitted loop" in audit.observable()["remaining"]
    assert len({name for name, _, _ in audit.bad_cases()}) == len(audit.bad_cases())


def full_forward_loop_diagnostic(n, S, nodes=12):
    """Numerical diagnostic only: complete ordered diagrams, with root-split logs."""
    nn, ss = mp.mpf(n), mp.mpf(S)
    gx, gw = mp.gauss_quadrature(nodes, "legendre")
    grid = [((gx[j] + 1) / 2, gw[j] / 2) for j in range(nodes)]

    def integrate_light(a, fun):
        if a > 4:
            beta = mp.sqrt(1 - 4 / a)
            r1, r2 = (1 - beta) / 2, (1 + beta) / 2
            points = [0, r1, r2, 1]
            getL = lambda x: a * (x - r1) * (x - r2)
        elif a == 4:
            points = [0, mp.mpf("0.5"), 1]
            getL = lambda x: (2 * x - 1) ** 2
        else:
            points = [0, mp.mpf("0.5"), 1]
            getL = lambda x: 1 - a * x * (1 - x)
        return mp.quad(lambda x: fun(getL(x)) if getL(x) != 0 else 0, points)

    channels = (ss, mp.mpf(0), 4 - ss)
    Dgap = nn - 2
    contact = -3 / Dgap + 2 / Dgap**2
    fixed = lambda a: contact + 1 / (nn - a)

    def db(a, b):
        if b == 0:
            return integrate_light(a, lambda L: closed_T_D(nn, L, 0)[1])
        if a == 0:
            return sum(
                weight * closed_T_D(nn, 1, b * x * (1 - x))[1] for x, weight in grid
            )
        return integrate_light(
            a,
            lambda L: sum(
                weight * closed_T_D(nn, L, b * x * (1 - x))[1] for x, weight in grid
            ),
        )

    result = 0
    for i, a in enumerate(channels):
        B = integrate_light(a, lambda L: -physical_log(L))
        C = integrate_light(a, lambda L: closed_T_D(nn, L, 0)[0])
        result += fixed(a) ** 2 * B / 2 + 2 * fixed(a) * C
        for j, b in enumerate(channels):
            if i != j:
                result += db(a, b)
    return result / (16 * mp.pi**2)


def test_full_forward_loop_numerics_and_complete_optical_cut():
    with mp.workdps(40):
        n, S = mp.mpf(256), mp.mpf(5)
        complete = full_forward_loop_diagnostic(n, S)
        # Independent angle-space optical theorem for the entire tree.
        Dgap = n - 2
        contact = -3 / Dgap + 2 / Dgap**2
        beta = mp.sqrt(1 - 4 / S)
        tree = lambda x: (
            contact
            + 1 / (n - S)
            + 1 / (n + (S - 4) * (1 - x) / 2)
            + 1 / (n + (S - 4) * (1 + x) / 2)
        )
        optical = beta / (64 * mp.pi) * mp.quad(lambda x: tree(x) ** 2, [-1, 0, 1])
        assert abs(mp.im(complete) - optical) < mp.mpf("1e-30")
        assert optical > 0
        universal = (
            -6 * mp.log(n) / n**2 + (-16 * mp.log(n) + mp.mpf(10) / 3) / n**3
        ) / (16 * mp.pi**2)
        remainder = complete - universal
        assert abs(remainder) < mp.mpf(10) ** 7 * S * S / (16 * mp.pi**2 * n**4)
        assert abs(mp.re(remainder)) > 0
        # Complete quadrature, not just its series, changes at the second mass.
        n2 = 2 * n
        complete2 = full_forward_loop_diagnostic(n2, S)
        universal2 = (
            -6 * mp.log(n2) / n2**2 + (-16 * mp.log(n2) + mp.mpf(10) / 3) / n2**3
        ) / (16 * mp.pi**2)
        assert abs(complete2 - universal2) < abs(remainder)
