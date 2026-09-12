"""Independent full flat covariance, original form factor and causal identities."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_full_flat_tensor_response import (
    analytic,
    audit,
    normalization,
    response,
)


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_all_exact_residuals(name, value):
    assert value == (
        s.zeros(*value.shape) if isinstance(value, s.MatrixBase) else s.S.Zero
    ), name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_all_proof_gates(name, value):
    assert value is True, name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_every_unsupported_input_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_original_frontier_and_actual_retained_flat_scope():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 86
    assert audit.frontier() == audit.previous.frontier()
    assert all(row["status"] != "COMPLETE" for row in audit.frontier())
    assert "GROWING_POLES" in audit.ITEM["status"]
    assert "NOT_CURVED_INSTABILITY_PHYSICAL_UV_NO_GO_OR_P8" in audit.ITEM["status"]


def test_independent_full_six_by_six_covariance_linear_solve():
    omega, lam = s.Integer(5), s.Integer(3)
    J = s.zeros(3).row_join(s.eye(3)).col_join((-s.eye(3)).row_join(s.zeros(3)))
    X = s.Matrix([[2, 1, -1], [1, 3, 2], [-1, 2, 5]])
    Y = s.Matrix([[3, -2, 1], [-2, 7, -1], [1, -1, 11]])
    XD = s.Matrix([[5, 3, 1], [3, 2, -2], [1, -2, 7]])
    YD = s.Matrix([[2, 1, -3], [1, 5, 2], [-3, 2, 4]])
    v = s.symbols("v0:21")
    tangent = s.zeros(6)
    n = 0
    for i in range(6):
        for j in range(i, 6):
            tangent[i, j] = tangent[j, i] = v[n]
            n += 1
    G = s.diag(X, Y)
    eq = lam * tangent - omega * (J * tangent - tangent * J) - (J * G - G * J) / 2
    solved = s.linsolve([eq[i, j] for i in range(6) for j in range(i, 6)], v)
    literal = tangent.subs(dict(zip(v, next(iter(solved)))))
    detector = -s.trace(s.diag(XD, YD) * literal) / 2
    wanted = omega * s.trace((YD - XD) * (Y - X)) / (2 * (lam**2 + 4 * omega**2))
    assert detector == wanted and detector != 0
    assert (literal[:3, 3:] - lam * (Y - X) / (2 * (lam**2 + 4 * omega**2))) == s.zeros(
        3
    )
    assert detector != 2 * wanted


@pytest.mark.parametrize("S", [s.Integer(4), s.Integer(17), s.Rational(71, 3)])
def test_original_second_vertex_and_three_derivative_marker_subtractions(S):
    p, marker, contact = s.symbols("p marker original_second_vertex")
    raw = contact + 7 / (S + marker**2 * p)
    sub = s.series(raw, marker, 0, 6).removeO().subs(marker, 1)
    rest = s.cancel(raw.subs(marker, 1) - sub)
    assert s.cancel(rest + 7 * p**3 / (S**3 * (S + p))) == 0
    assert s.expand(sub).coeff(contact) == 1
    wrong = sub - contact
    assert s.cancel(raw.subs(marker, 1) - wrong - rest) == contact


def test_independent_finite_action_force_sign_and_full_tree_composition():
    t = s.Symbol("t", real=True)
    h = s.Function("h")(t)
    m, k = s.symbols("m kappa", positive=True)
    lag = 5 * m * m * s.diff(h, t) ** 2 / 12 - s.diff(h, t, 2) ** 2 / 60
    euler = -s.diff(s.diff(lag, s.diff(h, t)), t) + s.diff(
        s.diff(lag, s.diff(h, t, 2)), t, 2
    )
    assert s.expand(euler + 5 * m * m * s.diff(h, t, 2) / 6 + s.diff(h, t, 4) / 30) == 0
    p, A = s.symbols("p original_A2")
    C0 = 16 * s.pi**2 * k
    C = C0 + 5 * m * m / 6
    O = p * (C + p * A) / C0
    G = C0 / (p * (C + p * A))
    assert s.cancel(O * G - 1) == 0
    assert s.cancel(G.subs(A, 0) - C0 / (C * p)) == 0
    assert s.cancel(C0 / (C + p * A) - 1 + (5 * m * m / 6 + p * A) / (C + p * A)) == 0
    assert s.cancel(O.subs(A, 0) - p) != 0


@pytest.mark.parametrize(
    "radius,angle", [(r, a) for r in (8, 32, 10000) for a in (0, 1, 3)]
)
def test_independent_massive_integral_exterior_branch_and_error_diagnostic(
    radius, angle
):
    # Finite precision diagnostic only; the all-sheet count and bound are written proofs.
    with mp.workdps(55):
        p = mp.mpf(radius) * mp.exp(1j * mp.pi * angle / 4)
        W = lambda y: y * y * (30 - 20 * y * y + 3 * y**4) / 30
        direct = mp.mpf(1) / 30 + mp.quad(
            lambda y: p * W(y) / (4 + p * (1 - y * y)), [0, 0.5, 0.9, 0.99, 1]
        )
        d = 1 + 4 / p
        root = mp.sqrt(d)
        B = -mp.mpf(172) / 225 + mp.mpf(19) * d / 30 - d * d / 10
        T = root * (30 - 20 * d + 3 * d * d) / 60
        stable = B + T * (mp.log(p) + 2 * mp.log((1 + root) / 2))
        assert abs(direct - stable) < mp.mpf("1e-45")
        error = abs(stable - mp.mpf(13) * mp.log(p) / 60 + mp.mpf(52) / 225)
        assert error <= 8 * (1 + mp.log(radius) + mp.pi) / radius


def test_independent_analytic_coefficient_derivative_identity():
    p, m = s.symbols("p m", positive=True)
    d = s.Symbol("d", positive=True)
    B = -s.Rational(172, 225) + 19 * d / 30 - d * d / 10
    T = s.sqrt(d) * (30 - 20 * d + 3 * d * d) / 60
    Lc = 2 * s.log((1 + s.sqrt(d)) / 2)
    expr = B + T * (s.log(p / (m * m)) + Lc)
    target = T - (d - 1) * (
        s.diff(B, d) + s.diff(T, d) * (s.log(p / (m * m)) + Lc) + T * s.diff(Lc, d)
    )
    actual = p * s.diff(expr.subs(d, 1 + 4 * m * m / p), p)
    assert s.simplify(actual - target.subs(d, 1 + 4 * m * m / p)) == 0


def test_independent_full_cut_sign_both_moment_coefficients_and_required_pair():
    p, tau, C = s.symbols("p tau C", positive=True)
    a, b, u, v, DA, U = s.symbols("a b u v DA U", real=True)
    z = a + s.I * b
    R = u + s.I * v
    Ebank = 1 / (-tau * (C - tau * (DA + s.I * s.pi * U)))
    assert (
        s.factor(
            -s.im(Ebank) / s.pi - U / ((C - tau * DA) ** 2 + s.pi**2 * tau**2 * U**2)
        )
        == 0
    )
    M0 = -1 / C - 2 * u
    M1 = 2 * (u * a - v * b)
    assert 1 / C + R + s.conjugate(R) + M0 == 0
    assert s.expand(R * z + s.conjugate(R * z) - M1) == 0
    assert (
        s.cancel(1 / (p + tau) - 1 / p + tau / p**2 - tau**2 / (p**2 * (p + tau))) == 0
    )
    assert s.cancel(1 / (p - z) - 1 / p - z / p**2 - z * z / (p**2 * (p - z))) == 0
    assert s.expand(2 * s.re(R * z) - M1) == 0
    assert (1 / C + s.Symbol("positive_cut_mass", positive=True)).is_positive


@pytest.mark.parametrize("q", [0, 7, 1004])
def test_complex_pole_geometry_and_prepared_source_annihilator(q):
    z = -3 + 4 * s.I
    plus = s.sqrt(z - q)
    minus = s.conjugate(plus)
    assert s.simplify(plus**2 + q - z) == 0
    assert s.re(plus).is_positive and s.im(plus).is_positive
    assert s.simplify(minus**2 + q - s.conjugate(z)) == 0
    lam = s.Symbol("lambda")
    source = lam - minus
    assert source.subs(lam, minus) == 0
    assert s.simplify(source.subs(lam, plus)) != 0
    assert s.simplify(s.re(plus) ** 2 - (s.Abs(z - q) + s.re(z - q)) / 2) == 0


def test_independent_forward_beta_moment_and_complex_wave_Duhamel_coefficients():
    y = s.Symbol("y", real=True)
    W = y * y * (30 - 20 * y * y + 3 * y**4) / 30
    assert s.integrate(W, (y, 0, 1)) == s.Rational(3, 14)
    # y=sin(theta), evaluated through elementary even sine beta moments.
    beta = [s.pi / 4, 3 * s.pi / 16, 5 * s.pi / 32]
    assert (30 * beta[0] - 20 * beta[1] + 3 * beta[2]) / 30 == 9 * s.pi / 64
    t, v = s.symbols("t v", nonnegative=True)
    q, z = s.symbols("q z")
    real = sum((-q) ** j * t ** (2 * j + 1) / s.factorial(2 * j + 1) for j in range(6))
    complex_wave = sum(
        (z - q) ** j * t ** (2 * j + 1) / s.factorial(2 * j + 1) for j in range(6)
    )
    duhamel = real + z * s.integrate(
        real.subs(t, t - v) * complex_wave.subs(t, v), (v, 0, t)
    )
    assert s.expand(s.series(duhamel - complex_wave, t, 0, 13).removeO()) == 0
    assert complex_wave.subs(t, 0) == 0 and s.diff(complex_wave, t).subs(t, 0) == 1


@pytest.mark.parametrize("call", [normalization.data, analytic.data, response.data])
def test_each_source_packet_keeps_exact_checks_and_true_gates(call):
    data = call()
    for value in data["checks"].values():
        vv = (
            value.applyfunc(s.cancel)
            if isinstance(value, s.MatrixBase)
            else s.cancel(value)
        )
        assert vv == (s.zeros(*vv.shape) if isinstance(vv, s.MatrixBase) else 0)
    assert all(bool(v) for v in data["gates"].values())


def test_all_sixteen_coarse_margins_and_physical_low_disk_scale():
    margins = analytic.data()["all_exact_positive_margins"]
    assert len(margins) == 16 and all(v > 0 for v in margins.values())
    assert analytic.data()["low_disk_A_bound"] == s.Rational(11, 105)
    assert analytic.data()["physical_low_disk_relative_error_bound"] < s.Rational(
        1, 10**796
    )
    assert analytic.RESIDUE_BOUND < 2
    assert normalization.C - normalization.C0 == 5 * normalization.MASS**2 / 6
