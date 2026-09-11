"""Independent full-function, tensor, source and labelled-vertex controls."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_canonical_affine_decoupling import (
    amplitude,
    audit,
    family,
    gravity,
    source,
)

ETA = s.diag(1, -1, -1, -1)


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_exact_residual(name):
    assert audit.residuals()[name] == 0


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_unsupported_scope_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@cache
def physical_functions(ratio):
    u, X = family.base.u, family.base.X
    R, F = family.at_base_coordinates(s.sqrt(ratio) * u, ratio * X, ratio)
    return R, F, s.diff(R, X), s.diff(R, X, 2)


@pytest.mark.parametrize(
    "ratio", [s.S.One, s.Integer(2), s.Integer(9), s.Integer(10) ** 8]
)
@pytest.mark.parametrize(
    "v,x",
    [(0, 0), (s.Rational(1, 3), s.Rational(1, 8)), (-2, 1), (2, -s.Rational(1, 8192))],
)
def test_full_functions_at_fixed_canonical_points(ratio, v, x):
    u, X = family.base.u, family.base.X
    R, F, RX, RXX = physical_functions(ratio)
    physical = {u: v / s.sqrt(ratio), X: x / ratio}
    baseline = {u: v, X: x}
    d = family.base.data()
    assert s.simplify(ratio * (R.subs(physical) - 1) - (d["R"].subs(baseline) - 1)) == 0
    assert s.simplify(ratio * F.subs(physical) - d["F"].subs(baseline)) == 0
    assert s.simplify(RX.subs(physical) - d["RX"].subs(baseline)) == 0
    if x:
        canonical_a3 = RX.subs(physical) / x
        assert s.simplify(canonical_a3 - d["A3"].subs(baseline)) == 0
    else:
        assert (
            s.simplify(RXX.subs(physical) / ratio - s.diff(d["R"], X, 2).subs(baseline))
            == 0
        )
    actual = s.N(R.subs(physical), 50)
    assert actual > s.Rational(1, 2) + s.Rational(1, 8 * family.N)
    assert actual < s.Rational(6, 5)


@pytest.mark.parametrize("ratio", [s.Integer(2), s.Integer(9)])
def test_complete_higher_scalar_couplings_not_dropped(ratio):
    u, X = family.base.u, family.base.X
    R, F, _, _ = physical_functions(ratio)
    # Units kappa0^2 stripped off, so enormous normalization does not hide a zero.
    fY3 = s.simplify(s.diff(F, X, 3).subs({u: 0, X: 0}) / (6 * ratio**2))
    a3Y = s.simplify(s.diff(R, X, 3).subs({u: 0, X: 0}) / (2 * ratio**2))
    assert fY3 == -s.Rational(109184, 25)
    assert a3Y == 3 * family.N
    assert fY3 != 0 and a3Y != 0


def test_wrong_fixed_physical_coordinates_change_the_family():
    u, X = family.base.u, family.base.X
    R, _, _, _ = physical_functions(s.Integer(2))
    at = {u: 0, X: s.Rational(1, 32)}
    wrong = s.N(2 * (R.subs(at) - 1), 60)
    fixed = s.N(family.base.data()["R"].subs(at) - 1, 60)
    assert abs(wrong - fixed) > s.Rational(1, 1000)


@pytest.mark.parametrize("ratio", [s.S.One, s.Integer(2), s.Integer(9)])
def test_nonzero_null_source_from_actual_base_germ(ratio):
    eps = s.Symbol("epsilon", positive=True)
    v = s.Rational(1, 3)
    a = -family.N / (1 + v * v) ** 3
    ap = 6 * family.N * v / (1 + v * v) ** 4
    X = eps / ratio
    R = 1 + a * ratio * X * X
    RX = 2 * a * ratio * X
    Ru = ap * ratio ** s.Rational(3, 2) * X * X
    gradient = s.Matrix([s.sqrt(1 + eps), 1, 0, 0]) / s.sqrt(ratio)
    H = s.diag(2, 1, 0, 0) / s.sqrt(ratio)
    Q = s.trace(ETA * H)
    Z = (gradient.T * ETA * H * ETA * gradient)[0]
    scalar = (R - 1) * (
        -3 * s.Rational(4, 5)
        + 3 * Ru / (4 * R)
        + Q / X
        + (-1 / X**2 + 3 * RX / (2 * R * X)) * Z
    )
    null = (gradient * scalar).applyfunc(lambda e: s.simplify(s.limit(e, eps, 0)))
    expected = -3 * a * s.Matrix([1, 1, 0, 0]) / ratio
    assert null == expected
    assert null != s.zeros(4, 1)


@pytest.mark.parametrize("u0", [0, 1])
@pytest.mark.parametrize("x0", [mp.mpf("-0.25"), mp.mpf("0.1"), mp.mpf("0.5")])
def test_regular_q_integral_from_both_sides(u0, x0):
    with mp.workdps(65):
        u = mp.mpf(u0)
        a = (1 + u * u) / 10
        ap = u / 5

        def q(x):
            R = 1 + a * x * x
            integrand = lambda t: (
                mp.sqrt(t)
                * (2 * a * t * x)
                * (ap * (t * x) ** 2)
                * (1 + a * (t * x) ** 2) ** (-mp.mpf(7) / 4)
            )
            return 3 * x / 4 * R ** (mp.mpf(3) / 4) * mp.quad(integrand, [0, 1])

        x = mp.mpf(x0)
        R = 1 + a * x * x
        residual = (
            mp.diff(q, x)
            + (1 / (2 * x) - 3 * (2 * a * x) / (4 * R)) * q(x)
            - 3 * (2 * a * x) * (ap * x * x) / (4 * R)
        )
        assert abs(residual) < mp.mpf("1e-55")
        assert mp.isfinite(q(x))


def coefficient_all_labels(expr, labels):
    out = s.expand(expr)
    for label in labels:
        out = out.coeff(label, 1)
    return s.expand(out)


@pytest.mark.parametrize(
    "momenta",
    [
        [(2, 1, 0, 0), (1, -1, 1, 0), (-1, 0, -1, 2)],
        [(1, 1, 0, 0), (1, 0, 1, 0), (0, 0, 0, 1)],
        [(3, 0, 1, 1), (-2, 1, 0, 1), (1, -2, 1, 0)],
    ],
)
def test_literal_off_shell_tensor_quartic_vertices(momenta):
    ps = [s.Matrix(v) for v in momenta]
    ps.append(-sum(ps, s.zeros(4, 1)))
    labels = s.symbols("e0:4")
    gradient = s.I * sum([labels[i] * ps[i] for i in range(4)], s.zeros(4, 1))
    H = -sum([labels[i] * ps[i] * ps[i].T for i in range(4)], s.zeros(4))
    Y = (gradient.T * ETA * gradient)[0]
    Q = s.trace(ETA * H)
    L3 = Q * (gradient.T * ETA * H * ETA * gradient)[0]
    L4 = (gradient.T * ETA * H * ETA * H * ETA * gradient)[0]
    gal = Y * (Q * Q - s.trace(ETA * H * ETA * H))
    dots = s.Matrix(4, 4, lambda i, j: (ps[i].T * ETA * ps[j])[0])
    lower, vertex, boundary, potential = amplitude.labelled_vertices(dots, 2, 3)
    assert coefficient_all_labels(2 * Y * Y, labels) == lower
    assert coefficient_all_labels(6 * (L4 - L3), labels) == vertex
    assert coefficient_all_labels(3 * gal, labels) == boundary
    assert vertex == boundary
    assert coefficient_all_labels(-(sum(labels) ** 4), labels) == potential


@pytest.mark.parametrize(
    "energy,cosine",
    [(s.Rational(5, 4), 0), (s.Rational(5, 3), s.Rational(3, 5)), (2, -1)],
)
def test_exact_on_shell_amplitude_from_physical_momenta(energy, cosine):
    p = s.sqrt(energy * energy - 1)
    sn = s.sqrt(1 - cosine * cosine)
    ps = [
        s.Matrix([energy, 0, 0, p]),
        s.Matrix([energy, 0, 0, -p]),
        s.Matrix([-energy, -p * sn, 0, -p * cosine]),
        s.Matrix([-energy, p * sn, 0, p * cosine]),
    ]
    dots = s.Matrix(4, 4, lambda i, j: (ps[i].T * ETA * ps[j])[0])
    assert all(s.simplify(dots[i, i] - 1) == 0 for i in range(4))
    S = (ps[0] + ps[1]).dot(ETA * (ps[0] + ps[1]))
    T = (ps[0] + ps[2]).dot(ETA * (ps[0] + ps[2]))
    U = (ps[0] + ps[3]).dot(ETA * (ps[0] + ps[3]))
    lam, gam = s.Rational(2, 7), s.Rational(3, 11)
    lo, de, _, po = amplitude.labelled_vertices(dots, lam, gam)
    expected = (
        2 * lam * sum((z - 2) ** 2 for z in (S, T, U)) + 3 * gam * S * T * U - 8 * gam
    )
    assert s.simplify(lo + de + po - expected) == 0


@pytest.mark.parametrize("index", [0, 1, 2])
def test_literal_metric_variation_gives_the_canonical_scalar_vertex(index):
    eps, e0, e1 = s.symbols("epsilon e0 e1", real=True)
    h = s.Matrix(
        4, 4, lambda i, j: s.Rational((min(i, j) + 1) * (max(i, j) + 2) + index, 7)
    )
    metric = ETA + 2 * eps * h
    inv = metric.inv()
    p = s.Matrix([2 + index, 1, 0, -1])
    q = s.Matrix([-1, 0, 1, index])
    grad = s.I * (p * e0 + q * e1)
    density = s.sqrt(-metric.det()) * ((grad.T * inv * grad)[0] - (e0 + e1) ** 2) / 2
    actual = s.diff(density, eps).subs(eps, 0).expand().coeff(e0).coeff(e1)
    dot = (p.T * ETA * q)[0]
    V = p * q.T + q * p.T - ETA * (dot + 1)
    expected = s.trace(ETA * h * ETA * V)
    assert s.simplify(actual - expected) == 0


@pytest.mark.parametrize("index", [0, 1, 2])
def test_dense_scalar_stress_projector_and_forward_boundary(index):
    E = s.Rational(5 + index, 2)
    cosine = s.Rational(1, 3 + index)
    p = s.sqrt(E * E - 1)
    sn = s.sqrt(1 - cosine * cosine)
    ps = [
        s.Matrix([E, 0, 0, p]),
        s.Matrix([E, 0, 0, -p]),
        s.Matrix([-E, -p * sn, 0, -p * cosine]),
        s.Matrix([-E, p * sn, 0, p * cosine]),
    ]

    def V(i, j):
        dot = (ps[i].T * ETA * ps[j])[0]
        return ps[i] * ps[j].T + ps[j] * ps[i].T - ETA * (dot + 1)

    A, B = V(0, 2), V(1, 3)
    assert s.simplify((ps[0] + ps[2]).T * ETA * A) == s.zeros(1, 4)
    actual = s.trace(ETA * A * ETA * B) - s.trace(ETA * A) * s.trace(ETA * B) / 2
    S = 4 * E * E
    T = (ps[0] + ps[2]).dot(ETA * (ps[0] + ps[2]))
    U = 4 - S - T
    assert s.simplify(actual - (2 - 2 * T - S * U)) == 0
    # Nonuniformity is real already in this tree observable.
    K = s.Symbol("kappa", positive=True)
    t = s.Symbol("t", negative=True)
    pole = (S * (4 - S - t) + 2 * t - 2) / (K * t)
    assert s.limit(pole, K, s.oo) == 0
    assert s.limit(pole.subs(t, -1 / K), K, s.oo) == S * S - 4 * S + 2


@pytest.mark.parametrize("scale", [10, 100, 10000])
def test_full_vector_source_and_contact_norm_majorants(scale):
    V = s.Matrix([1, -2, 3, 1])
    B = s.Matrix([2, 1, -1, 3])
    C = s.Matrix([-1, 3, 2, 1])
    k0 = s.Integer(10)
    K = s.Integer(scale)
    r = -s.Rational(2, 5) * k0
    R = 1 + r / K
    U = V + B / s.sqrt(K) + C / (K * R)
    upper = s.sqrt(V.dot(V)) + s.sqrt(B.dot(B)) / s.sqrt(k0) + 2 * s.sqrt(C.dot(C)) / k0
    assert s.N(s.sqrt(U.dot(U)), 60) < s.N(upper, 60)
    contact = (U.T * ETA * U)[0] / (2 * K)
    assert s.N(abs(contact), 60) <= s.N(upper**2 / (2 * K), 60)


@pytest.mark.parametrize(
    "k,x",
    [
        (family.K0, 0),
        (2 * family.K0, -s.Rational(1, 8192)),
        (10**6 * family.K0, s.Rational(11, 10)),
    ],
)
def test_allowed_full_family_points(k, x):
    assert family.require_scope(k, x) == (k, x, family.ZETA)


def test_all_fixed_higher_coefficients_nonzero_and_frontier_open():
    assert all(v != 0 for v in family.germs()["fixed_higher_coefficients"].values())
    assert len(audit.frontier()) == 9
    assert audit.matching()[-1]["status"].endswith(
        "NOT_QUANTUM_CONTOUR_TRUNCATION_OR_V_CLOSURE"
    )
    assert "interacting-state" in source.data()["not_inferred"].lower()


@pytest.mark.parametrize("index", [0, 1, 2])
def test_full_metric_Einstein_Gamma_density_converges_to_canonical_tensor(index):
    h = s.Matrix(
        4, 4, lambda i, j: s.Rational((min(i, j) + 1) * (max(i, j) + 1) + index, 20)
    )
    d = [
        [
            [
                s.Rational((a + 1) * (min(i, j) + 2) - (max(i, j) + 1) + index, 11)
                for j in range(4)
            ]
            for i in range(4)
        ]
        for a in range(4)
    ]
    expected = gravity.gamma_density(d)

    def full(eps):
        metric = ETA + 2 * eps * h
        inv = metric.inv()
        # Full Gamma/eps from actual inverse metric and derivative g=2eps*d.
        G = [
            [
                [
                    sum(
                        inv[r, k] * (d[m][k][n] + d[n][k][m] - d[k][m][n])
                        for k in range(4)
                    )
                    for n in range(4)
                ]
                for m in range(4)
            ]
            for r in range(4)
        ]
        value = (
            -s.sqrt(-metric.det())
            * sum(
                inv[m, n] * (G[r][m][q] * G[q][n][r] - G[r][m][n] * G[q][q][r])
                for m in range(4)
                for n in range(4)
                for r in range(4)
                for q in range(4)
            )
            / 2
        )
        return value

    assert s.simplify(full(s.S.Zero) - expected) == 0
    errors = []
    for den in (128, 256, 512):
        eps = s.Rational(1, den)
        error = abs(s.N(full(eps) - expected, 70))
        assert error / eps < 10**5
        errors.append(error)
    assert errors[1] < s.Rational(3, 5) * errors[0]
    assert errors[2] < s.Rational(3, 5) * errors[1]


def test_exact_cosmological_Einstein_boundary_before_limit():
    a, N, ad, add, Nd, K = s.symbols("a N a_dot a_ddot N_dot kappa", positive=True)
    R = -6 * (add / (a * N * N) + ad * ad / (a * a * N * N) - ad * Nd / (a * N**3))
    actual = -K * N * a**3 * R / 2
    exact_divergence = (
        3 * K * (2 * a * ad * ad / N + a * a * add / N - a * a * ad * Nd / (N * N))
    )
    assert s.factor(actual - exact_divergence + 3 * K * a * ad * ad / N) == 0
    eps, b, bd, n = s.symbols("epsilon b b_dot n", real=True)
    reduced = (-3 * K * a * ad * ad / N).subs(
        {K: 1 / eps**2, a: 1 + eps * b, ad: eps * bd, N: 1 + eps * n}, simultaneous=True
    )
    assert s.limit(reduced, eps, 0) == -3 * bd * bd
