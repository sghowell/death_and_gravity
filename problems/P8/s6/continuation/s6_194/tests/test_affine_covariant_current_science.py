"""Independent complete covariant Euler tensors and mixed metric bounds."""

from functools import cache
from math import factorial

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_covariant_current import (
    audit,
    canonical,
    current,
    metrics,
    tensors,
)


def mp_matrix(M):
    return mp.matrix([[mp.mpf(str(s.N(x, 90))) for x in row] for row in M.tolist()])


def mp_trace(M):
    return sum(M[j, j] for j in range(M.rows))


def fixture(which):
    F = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 17
    D = s.Matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]]) / 19
    spatial = [
        s.eye(3) if which == 0 else s.diag(1, 4, 9),
        F,
        D,
        (F + D) / 5,
        (2 * F - D) / 7,
    ]
    raw = []
    for j in range(5):
        metric = s.zeros(4)
        metric[0, 0] = -1 if j == 0 else 0
        metric[1:, 1:] = spatial[j]
        raw.append(metric)
    return spatial, raw, F, D


@cache
def geometric_fixture(which):
    return tensors.geometry(fixture(which)[1])


def literal_spatial_density(H, Hd, Hdd, invariant):
    inv = H**-1
    B = inv * Hd / 2
    Bd = inv * Hdd / 2 - inv * Hd * inv * Hd / 2
    X = Bd + B * B
    Y = Bd + mp_trace(B) * B
    R = mp_trace(X) + mp_trace(Y)
    Ric2 = mp_trace(X) ** 2 + mp_trace(Y * Y)
    value = {"R": R, "R2": R * R, "Ricci2": Ric2}[invariant]
    return mp.sqrt(mp.det(H)) * value


@pytest.mark.parametrize("which", (0, 1))
@pytest.mark.parametrize("invariant", ("R", "R2", "Ricci2"))
@pytest.mark.parametrize("direction", ("F", "D"))
def test_independent_full_local_metric_Euler_tensor(which, invariant, direction):
    with mp.workdps(85):
        spatial, _raw, F, D = fixture(which)
        direction_matrix = F if direction == "F" else D
        jets = list(map(mp_matrix, spatial))
        delta = mp_matrix(direction_matrix)

        def varied(t, b, slot):
            fields = [
                sum(
                    (t**j * jets[j + n] / factorial(j) for j in range(5 - n)),
                    mp.zeros(3),
                )
                for n in range(3)
            ]
            fields[slot] += b * delta
            return literal_spatial_density(*fields, invariant)

        actual = mp.mpf(0)
        for slot in range(3):
            actual += (-1) ** slot * mp.diff(
                lambda t, b, slot=slot: varied(t, b, slot),
                (mp.mpf(0), mp.mpf(0)),
                (slot, 1),
            )
        geometry = geometric_fixture(which)
        covariant = geometry["Euler_" + invariant]
        inv = spatial[0].inv()
        expected = s.sqrt(spatial[0].det()) * s.trace(
            inv * covariant[1:, 1:] * inv * direction_matrix
        )
        expected = mp.mpf(str(s.N(expected, 80)))
        assert abs(actual - expected) < mp.mpf("1e-55") * max(1, abs(expected))


def test_homogeneous_scalar_keeps_spatial_covariant_connections():
    t = s.Symbol("t", real=True)
    H, H1, H2, H3 = (
        s.Rational(1, 5),
        s.Rational(1, 7),
        s.Rational(1, 11),
        s.Rational(1, 13),
    )
    scale = s.exp(2 * H * t + H1 * t * t + H2 * t**3 / 3 + H3 * t**4 / 12)
    raw = []
    for n in range(5):
        metric = s.zeros(4)
        metric[0, 0] = -1 if n == 0 else 0
        metric[1:, 1:] = s.eye(3) * s.diff(scale, t, n).subs(t, 0)
        raw.append(metric)
    g = tensors.geometry(raw)
    assert g["scalar_Hessian"][1:, 1:] != s.zeros(3)
    assert (g["scalar_Hessian"][1:, 1:] + H * g["scalar"][1] * s.eye(3)).applyfunc(
        s.cancel
    ) == s.zeros(3)
    assert s.cancel(g["scalar_box"] + g["scalar"][2] + 3 * H * g["scalar"][1]) == 0
    assert g["scalar_box"] != -g["scalar"][2]


@pytest.mark.parametrize("which", (0, 1))
def test_all_covariant_Euler_tensors_are_symmetric(which):
    g = geometric_fixture(which)
    for name in ("Euler_R", "Euler_R2", "Euler_Ricci2", "scalar_Hessian", "Ricci_box"):
        assert (g[name] - g[name].T).applyfunc(s.cancel) == s.zeros(4)


@pytest.mark.parametrize("j", range(5))
@pytest.mark.parametrize("a", range(4))
def test_independent_absolute_mixed_exponential_majorant(j, a):
    t = s.Symbol("t", real=True)
    generating = s.exp(s.Rational(1, 100) * (s.exp(t) - 1) + a * t)
    expected = 2 * s.diff(generating, t, j).subs(t, 0)
    assert s.simplify(metrics.exponential(j, a) - expected) == 0


def test_relative_zero_jet_is_not_an_absolute_metric_bound():
    with mp.workdps(50):
        assert mp.exp(mp.mpf(1) / 100) > 1
        assert metrics.exponential(0, 0) == 2


@pytest.mark.parametrize(
    "j,a",
    (
        (0, 0),
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 1),
        (1, 3),
        (2, 0),
        (2, 2),
        (2, 3),
        (3, 1),
        (3, 3),
        (4, 0),
        (4, 2),
        (4, 3),
    ),
)
@pytest.mark.parametrize("inverse", (False, True))
def test_independent_actual_noncommuting_mixed_metric_jets(j, a, inverse):
    with mp.workdps(75):
        A = mp.matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 100
        B = mp.matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]]) / 100
        C = mp.matrix([[0, 1, 2], [1, 1, -1], [2, -1, -1]]) / 100
        jets = [A, B, C, (A + B) / 2, (B - C) / 3]

        def metric(t, e):
            direction = sum(
                (t**n * jets[n] / factorial(n) for n in range(5)), mp.zeros(3)
            )
            scale = (1 + t * t) ** 4
            return (1 / scale if inverse else scale) * mp.expm(
                (-1 if inverse else 1) * e * direction
            )

        value = mp.diff(metric, (mp.mpf(1) / 3, mp.mpf(1) / 200), (j, a))
        # Entry-l1 bound uses exactly the proved matrix-op to component conversion.
        l1 = sum(abs(value[r, c]) for r in range(3) for c in range(3))
        bound = 6 * metrics.constants()["raw_metric_operator_bounds"][inverse, j, a]
        assert l1 < mp.mpf(str(s.N(bound, 70)))
        assert mp.norm(A * B - B * A) > 0


@pytest.mark.parametrize("inverse", (False, True))
@pytest.mark.parametrize("shift", range(5))
def test_every_time_shifted_metric_seminorm_uses_all_required_jets(inverse, shift):
    c = metrics.constants()
    direct = sum(
        6
        * c["raw_metric_operator_bounds"][inverse, j + shift, a]
        / (s.factorial(j) * s.factorial(a))
        for j in range(5 - shift)
        for a in range(4)
    ) + (1 if shift == 0 else 0)
    assert (
        s.simplify(direct - c["all_time_shifted_component_jet_bounds"][inverse, shift])
        == 0
    )
    assert direct < metrics.G


def test_underestimated_metric_constant_is_rejected_by_actual_majorants():
    c = metrics.constants()["all_time_shifted_component_jet_bounds"]
    assert c[True, 3] > 100000
    assert max(c.values()) < 200000


@pytest.mark.parametrize("name", ("metrics", "tensors", "current", "canonical"))
def test_all_private_quantitative_packet_gates(name):
    module = {
        "metrics": metrics,
        "tensors": tensors,
        "current": current,
        "canonical": canonical,
    }[name]
    p = module.data()
    assert all(s.cancel(x) == 0 for x in p["checks"].values())
    assert all(p["gates"].values())


@pytest.mark.parametrize(
    "name,value",
    audit.residuals().items(),
    ids=lambda x: x if isinstance(x, str) else None,
)
def test_all_exact_residuals(name, value):
    assert s.cancel(value) == 0, name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=lambda x: x if isinstance(x, str) else None
)
def test_all_unsupported_scope_mutations(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_complete_counts_and_unchanged_frontier():
    assert len(audit.residuals()) == 22 and audit.scalar_entry_count() == 22
    assert len(audit.gates()) == 38 and all(audit.gates().values())
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 116
    assert len(audit.frontier()) == 9 and audit.matching()[-1] == audit.ITEM


@pytest.mark.parametrize("which", (0, 1, 2))
@pytest.mark.parametrize("momentum", (0, 1000, 10**16))
@pytest.mark.parametrize("detector", (0, 1))
def test_independent_constrained_clock_vertex_and_longitudinal_covariance(
    which, momentum, detector
):
    raw = [
        s.Matrix([1, 0, 0]),
        s.Matrix([1, 2, 3]) / s.sqrt(14),
        s.Matrix([-2, 1, 2]) / 3,
    ][which]
    n = raw
    k = momentum * n
    P = n * n.T
    D = [
        s.Matrix([[2, 1, -1], [1, -1, 2], [-1, 2, -1]]) / 10,
        s.Matrix([[1, -2, 1], [-2, 3, 0], [1, 0, 2]]) / 10,
    ][detector]
    a = s.Rational(25, 16)
    m = s.Integer(1000)
    cross = s.Matrix([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    VD = -a * m * m * D + cross.T * D * cross / a
    KD = D / a
    cT, cL, pT, pL = (
        s.Rational(2, 7),
        s.Rational(3, 11),
        s.Rational(5, 13),
        s.Rational(7, 17),
    )
    CA = cT * (s.eye(3) - P) + cL * P
    CP = pT * (s.eye(3) - P) + pL * P
    actual = s.trace(VD * CA) + s.trace(KD * CP)
    expression = canonical.data()["modewise_full_metric_vertex_contraction"]
    values = {
        "a": a,
        "mass": m,
        "k2": momentum**2,
        "nDn": (n.T * D * n)[0],
        "trD": s.trace(D),
        "cT": cT,
        "cL": cL,
        "pT": pT,
        "pL": pL,
    }
    expected = expression.subs(
        {symbol: values[str(symbol)] for symbol in expression.free_symbols}
    )
    assert s.simplify(actual - expected) == 0
    assert cT != cL and pT != pL


@pytest.mark.parametrize("which", (0, 1, 2))
def test_exact_rotational_quadrature_not_unregulated_divergent_cancellation(which):
    D = [
        s.diag(1, -1, 0),
        s.Matrix([[0, 1, 2], [1, 2, -1], [2, -1, -2]]),
        s.diag(1, 2, 3),
    ][which]
    # The six axes integrate every quadratic angular polynomial exactly.
    average = (
        sum(
            ((sign * s.eye(3)[:, j]).T * D * (sign * s.eye(3)[:, j]))[0]
            for j in range(3)
            for sign in (-1, 1)
        )
        / 6
    )
    assert average == s.trace(D) / 3


@pytest.mark.parametrize("epsilon", (-s.Rational(1, 100), 0, s.Rational(1, 100)))
@pytest.mark.parametrize("order", range(3))
def test_independent_noncommuting_detector_metric_jacobian(epsilon, order):
    with mp.workdps(70):
        F = mp.matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]]) / 100
        D = mp.matrix([[2, -1, 1], [-1, 0, 2], [1, 2, -2]]) / 100
        a = mp.mpf(25) / 16
        ep = mp.mpf(str(s.N(epsilon, 65)))
        actual = mp.diff(
            lambda e, b: a * a * mp.expm(e * F + b * D), (ep, mp.mpf(0)), (order, 1)
        )
        l1 = sum(abs(actual[i, j]) for i in range(3) for j in range(3))
        assert l1 < 12 * a * a
        assert mp.norm(F * D - D * F) > 0


@pytest.mark.parametrize(
    "epsilon",
    (
        -s.Rational(1, 100),
        -s.Rational(1, 200),
        0,
        s.Rational(1, 200),
        s.Rational(1, 100),
    ),
)
def test_Taylor_integral_both_signs_and_exact_zero(epsilon):
    e, x = s.symbols("e x", real=True)
    J = 3 + 5 * e + 7 * e**2 - 11 * e**3 + 13 * e**4
    remainder = J.subs(e, epsilon) - J.subs(e, 0) - epsilon * s.diff(J, e).subs(e, 0)
    integral = epsilon**2 * s.integrate(
        (1 - x) * s.diff(J, e, 2).subs(e, x * epsilon), (x, 0, 1)
    )
    assert s.factor(remainder - integral) == 0
    assert abs(remainder) <= epsilon**2 * 16 / 2
    assert (remainder == 0) == (epsilon == 0)


@pytest.mark.parametrize("kappa", (4, 100, 10**800))
def test_both_canonical_variations_from_literal_quadratic_functional(kappa):
    h, R, a = s.symbols("h R a", positive=True)
    gamma = 2 * h / s.sqrt(kappa)
    functional = R * gamma**2 / 2
    force = s.diff(functional, h) / a**3
    assert s.simplify(s.diff(force, h) - 4 * R / (a**3 * kappa)) == 0
    assert (
        s.simplify(s.diff(force, h) - 2 * R / (a**3 * s.sqrt(kappa))) != 0 or kappa == 4
    )


@pytest.mark.parametrize("raw", ([], [s.eye(4)] * 4, [s.eye(3)] * 5))
def test_full_coordinate_tensor_input_shape_required(raw):
    with pytest.raises(ValueError, match="five"):
        tensors.geometry(raw)


def test_scope_does_not_turn_derivative_loss_into_contraction():
    assert "C12-to-C0" in canonical.data()["norm_boundary"]
    assert "same-space" in canonical.data()["feedback_boundary"]
    assert "not a full spatial/mixed inverse" in current.data()["not_a_solution"]
    assert current.data()["physical_current_divided_by_kappa"] == tuple(
        x / s.Integer(10) ** 800 for x in current.CURRENT
    )
