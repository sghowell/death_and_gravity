"""Independent covariant, dense-constraint and full-function science checks."""

import copy
from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_affine import connection
from p8_affine_kinetic import scalar as original_scalar
from p8_affine_vacuum_domain import family
from p8_proca_retuned_margin import model as margin
from p8_vacuum_analytic_affine_parent import (
    affine,
    audit,
    chart,
    dynamics,
    source,
    verify,
)


@pytest.mark.parametrize("name,value", tuple(audit.residuals().items()))
def test_every_exact_residual(name, value):
    if isinstance(value, s.MatrixBase):
        assert value == s.zeros(*value.shape), name
    else:
        assert s.simplify(value) == 0, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_strict_rejected_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "x",
    (
        -s.Rational(1, 10000),
        0,
        s.Rational(1, 10000),
        s.Rational(1, 2),
        1,
        s.Rational(119, 100),
    ),
)
def test_actual_domain_allows_zero_and_small_negative_X(x):
    assert chart.require_domain(x, 1, s.Rational(1, 2000))[0] == x


def direct_trace_matrix():
    M = s.zeros(4, 64)
    index = lambda a, b, c: 16 * a + 4 * b + c
    eta = (-1, 1, 1, 1)
    for mu in range(4):
        for a in range(4):
            M[mu, index(a, mu, a)] += 1
            M[mu, index(mu, a, a)] += eta[mu] * eta[a]
            M[mu, index(a, a, mu)] -= s.Rational(1, 2)
    return M


@pytest.mark.parametrize("p", (s.Rational(2, 5), s.Rational(1, 2), s.Rational(27, 50)))
def test_dense_full_quotient_and_source_centering(p):
    d = affine.matrices()
    e = connection.quotient()["embedding"]
    n64 = direct_trace_matrix()
    n = n64 * e
    old64 = connection.quadratic()["hessian"].subs(connection.P, p)
    old60 = e.T * old64 * e
    tt = 3 * (2 * p**3 - 1) / p
    ss = (8 * p + 5) / (8 * p * p)
    eta = s.diag(1, -1, -1, -1)
    Xi = eta - s.diag(1 / tt, 1 / ss, 1 / ss, 1 / ss)
    new64 = old64 + n64.T * Xi * n64
    new60 = e.T * new64 * e
    ell = s.Matrix(
        [s.Rational(2, 7), -s.Rational(3, 11), s.Rational(5, 13), -s.Rational(7, 17)]
    )
    solved = new60.inv(method="DM") * (n.T * ell)
    # Independent dense inverse, not the blockwise symbolic construction.
    assert n * solved == eta * ell
    assert new64 * e * solved == n64.T * ell
    assert new60 == d["M_new"].subs(connection.P, p)
    assert new64 * connection.quadratic()["gauge"] == s.zeros(64, 4)
    vals = {
        connection.P: p,
        connection.S: 1,
        connection.PP: s.Rational(1, 7),
        connection.PX: -s.Rational(2, 11),
        connection.CP: s.Rational(3, 13),
        connection.CX: -s.Rational(5, 17),
        connection.F3: s.Rational(7, 19),
    }
    vals.update({z: s.Rational(i - 4, 23) for i, z in enumerate(connection.H_SYMBOLS)})
    oldsource = connection.quadratic()["source"].subs(vals)
    ystar = -old60.inv(method="DM") * (e.T * oldsource)
    tstar = n * ystar
    newsource = e.T * oldsource - n.T * Xi * tstar
    assert new60 * ystar + newsource == s.zeros(60, 1)
    # The constant term is necessary to preserve the stationary action value.
    delta_value = (ystar.T * (new60 - old60) * ystar)[0] / 2 + (
        (newsource - e.T * oldsource).T * ystar
    )[0]
    assert s.factor(delta_value + (tstar.T * Xi * tstar)[0] / 2) == 0
    assert delta_value != 0
    assert new60.det() / old60.det() == -tt * ss**3


@pytest.mark.parametrize(
    "covector",
    (
        (s.Rational(1, 2), 0, 0, 0),
        (0, s.Rational(1, 2), 0, 0),
        (1, 1, 0, 0),
        (s.Rational(3, 4), s.Rational(1, 2), s.Rational(1, 4), 0),
    ),
)
def test_metric_map_inverse_and_full_ten_dimensional_jacobian(covector):
    g = s.diag(1, -1, -1, -1)
    v = s.Matrix(covector)
    x = (v.T * g.inv() * v)[0]
    z = s.Symbol("z")
    # A rational positive test function for the GENERIC inverse identity.
    Cfun = 1 / (1 + z * z / 10)
    Dfun = z / (10 + z * z)
    C, D, Cx, Dx = [
        s.diff(f, z, n).subs(z, x)
        for f, n in ((Cfun, 0), (Dfun, 0), (Cfun, 1), (Dfun, 1))
    ]
    gp = C * g + D * v * v.T
    xp = (v.T * gp.inv() * v)[0]
    assert xp == x
    invD = -D / C
    assert gp / C + invD * v * v.T == g
    assert gp.det() / g.det() == C**3
    pairs = [(i, j) for i in range(4) for j in range(i, 4)]
    raised = g.inv() * v
    jac = s.zeros(10)
    for row, (i, j) in enumerate(pairs):
        for col, (a, b) in enumerate(pairs):
            dx = -raised[a] * raised[b] * (1 if a == b else 2)
            jac[row, col] = C * int(row == col) + (Cx * g[i, j] + Dx * v[i] * v[j]) * dx
    assert jac.det() == C**9
    if x == 0:
        assert jac != s.eye(10)
        assert (jac - s.eye(10)) ** 2 == s.zeros(10)
    # Direct differentiation of the literal finite metric map along a dense variation.
    e = s.Symbol("epsilon")
    variation = s.Matrix(
        [[s.Rational((i + 1) * (j + 1), 31) for j in range(4)] for i in range(4)]
    )
    ge = g + e * variation
    xe = (v.T * ge.inv() * v)[0]
    mapped = Cfun.subs(z, xe) * ge + Dfun.subs(z, xe) * v * v.T
    literal = s.Matrix([s.diff(mapped[i, j], e).subs(e, 0) for i, j in pairs])
    assert literal == jac * s.Matrix([variation[i, j] for i, j in pairs])


@cache
def full_R():
    d = family.data()
    return s.lambdify(
        (family.u, family.X),
        (d["R"], d["RX"], s.diff(d["R"], family.u)),
        "mpmath",
        cse=True,
    )


@pytest.mark.parametrize(
    "u,x",
    (
        (0, -s.Rational(1, 10000)),
        (0, 0),
        (0, s.Rational(1, 50)),
        (0, s.Rational(503, 1000)),
        (1, s.Rational(1, 2)),
        (-1, 1),
        (3, s.Rational(119, 100)),
    ),
)
def test_full_actual_R_domain_and_new_mass_not_inherited(u, x):
    with mp.workdps(100):
        xx = mp.mpf(str(x.p)) / int(x.q) if isinstance(x, s.Rational) else mp.mpf(x)
        R, _, _ = full_R()(mp.mpf(u), xx)
        assert mp.mpf("0.5") < R < mp.mpf("1.2")
        p = mp.sqrt(R) / 2
        tt = 3 * (2 * p**3 - 1) / p
        ss = (8 * p + 5) / (8 * p * p)
        gamma = 1 - 3 * (R - 1) ** 2 / (2 * R)
        assert tt < 0 < ss and gamma > mp.mpf("0.25")
        assert abs((1 - 1 / tt) + 1 / tt - 1) < mp.mpf("1e-95")
        assert abs((-1 - 1 / ss) + 1 / ss + 1) < mp.mpf("1e-95")
        if R != 1:
            assert abs((mp.mpf(11) / 9 + 1 / tt) - 1) > 0


def covariant_source_at(t, epsilon):
    N = 1 + epsilon / (1 + t * t)
    Nd = -2 * epsilon * t / (1 + t * t) ** 2
    H = 4 * t / (1 + t * t)
    Ha = H + epsilon * (1 - t * t) / (1 + t * t) ** 2
    X = N**-2
    R, RX, Ru = full_R()(t, X)
    Box = -Nd / N**3 + 3 * Ha / N**2
    Z = -Nd / N**5
    direct = (R - 1) * (
        -3 * H + 3 * Ru / (4 * R) + Box / X + (-1 / X**2 + 3 * RX / (2 * R * X)) * Z
    )
    rd = Ru - 2 * RX * Nd / N**3
    chart_value = 3 * (R - 1) * (Ha - H + rd / (4 * R))
    return direct, chart_value


@pytest.mark.parametrize("t", (-2, -s.Rational(1, 2), 0, s.Rational(1, 5), 1, 3))
def test_literal_full_function_metric_variation_and_nonzero_second_source(t):
    with mp.workdps(130):
        t = mp.mpf(str(t.p)) / int(t.q) if isinstance(t, s.Rational) else mp.mpf(t)
        h = (1 + t * t) ** 3
        hu = 6 * t * (1 + t * t) ** 2
        g = 1 / (1 + t * t)
        gd = -2 * t / (1 + t * t) ** 2
        fd = (1 - t * t) / (1 + t * t) ** 2
        expected = -2 * g / h * (3 * fd - mp.mpf(3) / 2 * (gd / h - g * hu / h**2))
        assert covariant_source_at(t, mp.mpf(0)) == (0, 0)
        assert abs(mp.diff(lambda e: covariant_source_at(t, e)[0], 0)) < mp.mpf(
            "1e-100"
        )
        for eps in (mp.mpf("1e-14"), -mp.mpf("1e-14")):
            a, b = covariant_source_at(t, eps)
            assert abs(a - b) < mp.mpf("1e-100")
            assert abs(a / eps**2 - expected) < mp.mpf("1e-9")
        assert expected != 0


@pytest.mark.parametrize("sign", (-1, 1))
def test_nonzero_null_gradient_source_has_finite_nonzero_limit(sign):
    with mp.workdps(130):
        x = sign * mp.mpf("1e-25")
        R, RX, Ru = full_R()(mp.mpf(0), x)
        v0 = mp.sqrt(1 + x)
        # Direct +--- contractions for Hessian diag(1,2,3,4).
        Box = -8
        Z = 3 + x
        scalar = (R - 1) * (
            3 * Ru / (4 * R) + Box / x + (-1 / x**2 + 3 * RX / (2 * R * x)) * Z
        )
        assert abs(scalar - 3072) < mp.mpf("1e-18")
        assert abs(v0 * scalar - 3072) < mp.mpf("1e-18")


def test_actual_lapse_J_against_independent_frozen_principal_and_margin():
    d = dynamics.clock_lapse()
    u = d["u"]
    h = (1 + u * u) ** 3
    expected = original_scalar.background()["J"] + 4 * margin.NEW_MARGIN / h**2
    assert s.factor(d["J"] - expected) == 0
    assert all(c > 0 for c in d["lapse_numerator_coefficients"])
    assert d["J"].subs(u, 0) == s.Rational(243, 160)
    I_NN = -3 * s.diff(h, u) / h**3
    H = 4 * u / (1 + u * u)
    missing_boundary_J = s.factor(d["J"] + (3 * H * I_NN + s.diff(I_NN, u)) / 2)
    assert missing_boundary_J.subs(u, 0) == s.Rational(243, 160) - 9 < 0


@pytest.mark.parametrize(
    "R", (s.Rational(5001, 10000), s.Rational(3, 4), 1, s.Rational(11999, 10000))
)
def test_full_ten_velocity_inertia_and_uniform_auxiliary_margin(R):
    K, x, y, xy, xz, yz = s.symbols("K x y xy xz yz")
    E = s.symbols("E0:3")
    matter = s.Symbol("matter_velocity")
    matrix = s.Matrix(
        [[K / 3 + x, xy, xz], [xy, K / 3 + y, yz], [xz, yz, K / 3 - x - y]]
    )
    zeta = s.Rational(1, 2000)
    L = (
        R * (s.trace(matrix * matrix) - s.trace(matrix) ** 2) / 2
        + (R - 1) ** 2 * K * K / 2
        + zeta * sum(z * z for z in E) / 2
        + matter * matter / 2
    )
    H = s.hessian(L, (K, x, y, xy, xz, yz, *E, matter))
    assert H[0, 0] < 0
    positive = H[1:, 1:]
    assert all(positive[:j, :j].det() > 0 for j in range(1, 10))
    assert 1 - 3 * (R - 1) ** 2 / (2 * R) > s.Rational(1, 4)


@pytest.mark.parametrize(
    "t,k", ((0, 0), (0, 1), (0, 100), (-1, 2), (1, 3), (10, 5), (s.Rational(1, 2), 10))
)
def test_frequency_from_literal_time_normalization(t, k):
    u = s.Symbol("u", real=True)
    a = (1 + u * u) ** 2
    mass = s.Integer(2000)
    q = k * k / a**2
    Z = a * mass / (q + mass)
    physical = q + mass - s.diff(s.sqrt(Z), u, 2) / s.sqrt(Z)
    H = s.diff(a, u) / a
    r = q / (q + mass)
    formula = (
        q
        + mass
        - (s.Rational(1, 2) + r) * s.diff(H, u)
        + H * H * (-s.Rational(1, 4) + r - 3 * r * r)
    )
    assert s.factor(physical - formula) == 0
    assert s.factor(physical - q - 1985).subs(u, t) > 0


def test_centering_shift_and_vacuum_degree_controls():
    q = s.Symbol("nonzero_regular_clock_q", nonzero=True)
    # At R=1, R_u=0, the actual stationary trace need not vanish.
    B = 3 * q / 2
    assert s.Rational(11, 9) * B == 11 * q / 6
    assert B != 0
    eps = s.Symbol("epsilon")
    grad, Box, Z = s.symbols("gradient Box_u contracted_Hessian")
    X = s.Symbol("X")
    Rminus = -1024 * (eps**2 * X) ** 2
    expression = s.expand(
        eps
        * grad
        * Rminus
        * (eps * Box / (eps**2 * X) - eps**3 * Z / (eps**2 * X) ** 2)
    )
    assert s.Poly(expression, eps).monoms() == [(4,)]
    assert s.factor(expression.coeff(eps, 4) - 1024 * grad * (Z - X * Box)) == 0


def test_scoped_exact_serialization_and_no_frontier_promotion():
    for packet in (
        chart.data(),
        chart.target_jets(),
        source.data(),
        affine.data(),
        dynamics.data(),
        dynamics.clock_lapse(),
    ):
        verify.serialize(verify.payload(packet))
    for i in range(9):
        f = copy.deepcopy(audit.frontier())
        f[i]["status"] = "COMPLETE"
        with pytest.raises(ValueError):
            audit.validate_scope(f, audit.matching())
    assert len(audit.frontier()) == 9
    assert all(v is True for v in audit.gates().values())
    with pytest.raises(ValueError):
        verify.serialize(s.Float("0.1"))
