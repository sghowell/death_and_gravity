"""Independent literal 4D curvature, full time-jet and normalization checks."""

from functools import cache
from itertools import product

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_local_spatial_hessian import (
    audit,
    bounds,
    helicities,
    local,
    proper,
)

Z = (s.Integer(0),) * 3
ONE = (s.Integer(1), s.Integer(0), s.Integer(0))
indices = tuple(product(range(4), repeat=4))
a = s.Symbol("a", positive=True)
a1, a2 = s.symbols("a1 a2", real=True)
x0, x1, x2 = s.symbols("x0 x1 x2")
p = s.Symbol("p", real=True)


def plus(*args):
    return tuple(s.expand(sum(arg[j] for arg in args)) for j in range(3))


def scale(c, value):
    return tuple(s.expand(c * v) for v in value)


def times(left, right):
    if left == Z or right == Z:
        return Z
    return tuple(
        s.expand(sum(left[k] * right[j - k] for k in range(j + 1))) for j in range(3)
    )


def constant(value):
    return (s.sympify(value), s.Integer(0), s.Integer(0))


def contract_matrix(left, right):
    return [
        [plus(*(times(left[i][k], right[k][j]) for k in range(4))) for j in range(4)]
        for i in range(4)
    ]


def geometric_density(T, cosine, sine, jets=None, momentum=None):
    T = s.Matrix(T)
    P = s.Matrix([0, 0, p]) if momentum is None else s.Matrix(momentum)
    M0, M1, M2 = (
        (x0 * T, x1 * T, x2 * T) if jets is None else tuple(map(s.Matrix, jets))
    )
    square = M0 * M0
    mixed_jet = M0 * M1 + M1 * M0
    second_jet = (M0 * M2 + M2 * M0) / 2 + M1 * M1

    def spatial_e(sign):
        return [
            [
                (
                    s.KroneckerDelta(i, j),
                    sign * cosine * M0[i, j],
                    cosine * cosine * square[i, j] / 2,
                )
                for j in range(3)
            ]
            for i in range(3)
        ]

    E = spatial_e(1)
    Ei = spatial_e(-1)
    Et = [
        [
            (0, cosine * M1[i, j], cosine * cosine * mixed_jet[i, j] / 2)
            for j in range(3)
        ]
        for i in range(3)
    ]
    Ett = [
        [(0, cosine * M2[i, j], cosine * cosine * second_jet[i, j]) for j in range(3)]
        for i in range(3)
    ]

    def Es(direction, sign=1):
        return [
            [
                (
                    0,
                    -sign * sine * P[direction] * M0[i, j],
                    -cosine * sine * P[direction] * square[i, j],
                )
                for j in range(3)
            ]
            for i in range(3)
        ]

    def Ets(direction):
        return [
            [
                (
                    0,
                    -sine * P[direction] * M1[i, j],
                    -cosine * sine * P[direction] * mixed_jet[i, j],
                )
                for j in range(3)
            ]
            for i in range(3)
        ]

    def Ess(left, right):
        return [
            [
                (
                    0,
                    -cosine * P[left] * P[right] * M0[i, j],
                    (sine * sine - cosine * cosine) * P[left] * P[right] * square[i, j],
                )
                for j in range(3)
            ]
            for i in range(3)
        ]

    metric = [[Z for _ in range(4)] for _ in range(4)]
    inverse = [[Z for _ in range(4)] for _ in range(4)]
    first = [[[Z for _ in range(4)] for _ in range(4)] for _ in range(4)]
    inverse_first = [[[Z for _ in range(4)] for _ in range(4)] for _ in range(4)]
    second = [
        [[[Z for _ in range(4)] for _ in range(4)] for _ in range(4)] for _ in range(4)
    ]
    metric[0][0] = constant(-a * a)
    inverse[0][0] = constant(-(a**-2))
    first[0][0][0] = constant(-2 * a * a1)
    inverse_first[0][0][0] = constant(2 * a1 / a**3)
    second[0][0][0][0] = constant(-2 * (a1 * a1 + a * a2))
    for i, j in product(range(3), repeat=2):
        I, J = i + 1, j + 1
        metric[I][J] = scale(a * a, E[i][j])
        inverse[I][J] = scale(a**-2, Ei[i][j])
        first[0][I][J] = plus(scale(2 * a * a1, E[i][j]), scale(a * a, Et[i][j]))
        inverse_first[0][I][J] = plus(
            scale(-2 * a1 / a**3, Ei[i][j]),
            scale(
                a**-2, (0, -cosine * M1[i, j], cosine * cosine * mixed_jet[i, j] / 2)
            ),
        )
        second[0][0][I][J] = plus(
            scale(2 * (a1 * a1 + a * a2), E[i][j]),
            scale(4 * a * a1, Et[i][j]),
            scale(a * a, Ett[i][j]),
        )
        for l in range(3):
            first[l + 1][I][J] = scale(a * a, Es(l)[i][j])
            inverse_first[l + 1][I][J] = scale(a**-2, Es(l, -1)[i][j])
            mixed = plus(scale(2 * a * a1, Es(l)[i][j]), scale(a * a, Ets(l)[i][j]))
            second[0][l + 1][I][J] = second[l + 1][0][I][J] = mixed
            for k in range(3):
                second[l + 1][k + 1][I][J] = scale(a * a, Ess(l, k)[i][j])
    identity = contract_matrix(metric, inverse)
    assert all(
        identity[i][j] == (ONE if i == j else Z) for i, j in product(range(4), repeat=2)
    )
    gamma = {}
    gamma_first = {}
    for i, j, k in product(range(4), repeat=3):
        gamma[i, j, k] = scale(
            s.Rational(1, 2),
            plus(
                *(
                    times(
                        inverse[i][e],
                        plus(first[j][e][k], first[k][e][j], scale(-1, first[e][j][k])),
                    )
                    for e in range(4)
                )
            ),
        )
        for l in range(4):
            gamma_first[l, i, j, k] = scale(
                s.Rational(1, 2),
                plus(
                    *(
                        plus(
                            times(
                                inverse_first[l][i][e],
                                plus(
                                    first[j][e][k],
                                    first[k][e][j],
                                    scale(-1, first[e][j][k]),
                                ),
                            ),
                            times(
                                inverse[i][e],
                                plus(
                                    second[l][j][e][k],
                                    second[l][k][e][j],
                                    scale(-1, second[l][e][j][k]),
                                ),
                            ),
                        )
                        for e in range(4)
                    )
                ),
            )
    Riem = {}
    for i, j, k, l in indices:
        Riem[i, j, k, l] = plus(
            gamma_first[k, i, l, j],
            scale(-1, gamma_first[l, i, k, j]),
            *(
                plus(
                    times(gamma[i, k, e], gamma[e, l, j]),
                    scale(-1, times(gamma[i, l, e], gamma[e, k, j])),
                )
                for e in range(4)
            ),
        )
    Ric = [
        [plus(*(Riem[e, i, e, j] for e in range(4))) for j in range(4)]
        for i in range(4)
    ]
    R = plus(*(times(inverse[i][j], Ric[i][j]) for i, j in product(range(4), repeat=2)))
    mixed = contract_matrix(inverse, Ric)
    Ric2 = plus(
        *(times(mixed[i][j], mixed[j][i]) for i, j in product(range(4), repeat=2))
    )
    cov = {
        key: plus(*(times(metric[key[0]][e], Riem[(e,) + key[1:]]) for e in range(4)))
        for key in indices
    }
    raised = cov
    for slot in range(4):
        old = raised
        raised = {
            key: plus(
                *(
                    times(
                        inverse[key[slot]][e], old[key[:slot] + (e,) + key[slot + 1 :]]
                    )
                    for e in range(4)
                )
            )
            for key in indices
        }
    Riem2 = plus(*(times(cov[key], raised[key]) for key in indices))
    R2 = times(R, R)
    Weyl = plus(Riem2, scale(-2, Ric2), scale(s.Rational(1, 3), R2))
    result = {
        name: tuple(s.factor(a**4 * v) for v in value)
        for name, value in {
            "R_old": R,
            "R_old_squared": R2,
            "Weyl_squared": Weyl,
        }.items()
    }
    assert s.factor(result["R_old"][0] - 6 * a * a2) == 0
    assert s.factor(result["R_old_squared"][0] - 36 * a2 * a2 / a**2) == 0
    assert result["Weyl_squared"][0] == result["Weyl_squared"][1] == 0
    assert s.factor(result["R_old"][1] + a * a * (P.T * M0 * P)[0] * cosine) == 0
    assert (
        s.factor(result["R_old_squared"][1] + 12 * a2 * (P.T * M0 * P)[0] * cosine / a)
        == 0
    )
    return result


@cache
def averaged_quadratic(kind):
    T = physical_test_basis()[kind]
    left = geometric_density(T, s.Integer(1), s.Integer(0))
    right = geometric_density(T, s.Integer(0), s.Integer(1))
    return {name: s.factor((left[name][2] + right[name][2]) / 2) for name in left}


def physical_test_basis():
    return {
        "tensor": s.diag(1, -1, 0),
        "tensor_cross": s.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
        "vector": s.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
        "vector_cross": s.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
        "scalar": s.diag(-1, -1, 2),
    }


A3, A4 = s.symbols("a3 a4")
X3, X4 = s.symbols("x3 x4")


def clock_derivative(expr):
    substitutions = (
        (a, a1),
        (a1, a2),
        (a2, A3),
        (A3, A4),
        (x0, x1),
        (x1, x2),
        (x2, X3),
        (X3, X4),
    )
    return s.expand(sum(s.diff(expr, left) * right for left, right in substitutions))


def euler(expr):
    return s.expand(
        s.diff(expr, x0)
        - clock_derivative(s.diff(expr, x1))
        + clock_derivative(clock_derivative(s.diff(expr, x2)))
    )


@pytest.mark.parametrize("kind", tuple(physical_test_basis()))
@pytest.mark.parametrize("curvature", ("R_old", "R_old_squared", "Weyl_squared"))
def test_literal_full_four_dimensional_compact_action_Hessian(kind, curvature):
    T = physical_test_basis()[kind]
    actual = averaged_quadratic(kind)[curvature]
    jets = [x0 * T, x1 * T, x2 * T, X3 * T, X4 * T]
    H = local.operators(s.Matrix([0, 0, p]), jets, a, a1, a2, A3)[curvature]
    assert s.factor(euler(actual) - s.trace(T * H) / 2) == 0


@pytest.mark.parametrize("kind", tuple(physical_test_basis()))
def test_literal_curved_Weyl_density_has_no_conformal_scale_jets(kind):
    value = averaged_quadratic(kind)["Weyl_squared"]
    assert not value.has(a, a1, a2)


@pytest.mark.parametrize("kind", tuple(physical_test_basis()))
def test_exact_compact_Einstein_and_R_squared_quadratic_densities(kind):
    T = physical_test_basis()[kind]
    P = s.Matrix([0, 0, p])
    q = p * p
    norm = s.trace(T * T)
    longitudinal = (P.T * T * T * P)[0]
    form = norm * x1 * x1 - (q * norm - 2 * longitudinal) * x0 * x0
    expected_R = a * a * form / 8
    expected_R2 = (P.T * T * P)[0] ** 2 * x0 * x0 / 2 + 3 * a2 * form / (2 * a)
    actual = averaged_quadratic(kind)
    assert s.factor(actual["R_old"] - expected_R) == 0
    assert s.factor(actual["R_old_squared"] - expected_R2) == 0


@cache
def noncommuting_full_geometry(case):
    D = s.Matrix([[1, 2, -1], [2, -3, 1], [-1, 1, 2]])
    G = s.Matrix([[2, -1, 3], [-1, 1, 2], [3, 2, -3]])
    assert D * G != G * D
    M0, M1, M2 = D, G, D + 2 * G
    P = p * s.Matrix(((2, -1, 3), (0, 0, 1), (0, 0, 0))[case])
    left = geometric_density(D, 1, 0, (M0, M1, M2), P)
    right = geometric_density(D, 0, 1, (M0, M1, M2), P)
    actual = {name: s.factor((left[name][2] + right[name][2]) / 2) for name in left}
    q = (P.T * P)[0]
    norm = lambda M: s.trace(M * M)
    dot = lambda M, N: s.trace(M * N)
    v = lambda M, N: (P.T * M * N * P)[0]
    form = norm(M1) - q * norm(M0) + 2 * v(M0, M0)
    predicted = {
        "R_old": a * a * form / 8,
        "R_old_squared": (P.T * M0 * P)[0] ** 2 / 2 + 3 * a2 * form / (2 * a),
        "Weyl_squared": norm(M2) / 4
        - q * norm(M1)
        + s.Rational(3, 2) * v(M1, M1)
        - q * dot(M0, M2) / 2
        + v(M0, M2)
        + q * q * norm(M0) / 4
        - q * v(M0, M0) / 2
        + (P.T * M0 * P)[0] ** 2 / 6,
    }
    return actual, predicted


@pytest.mark.parametrize("case", range(3))
@pytest.mark.parametrize("curvature", ("R_old", "R_old_squared", "Weyl_squared"))
def test_literal_noncommuting_matrix_time_jets_and_generic_spatial_direction(
    case, curvature
):
    actual, predicted = noncommuting_full_geometry(case)
    assert s.factor(actual[curvature] - predicted[curvature]) == 0


@pytest.mark.parametrize("kind", tuple(physical_test_basis()))
def test_Weyl_compact_boundary_term_not_dropped_pointwise(kind):
    T = physical_test_basis()[kind]
    P = s.Matrix([0, 0, p])
    q = p * p
    M0 = x0 * T
    M1 = x1 * T
    M2 = x2 * T
    B0 = helicities.B(P, M0)
    C0 = helicities.C(P, M0)
    compact = (
        s.trace(M2 * M2) / 4
        - s.trace(M1 * (2 * q * M1 - 2 * helicities.B(P, M1))) / 4
        + s.trace(M0 * (q * q * M0 - 2 * q * B0 + s.Rational(2, 3) * C0)) / 4
    )
    boundary_term = s.trace(M0 * (-q * M1 / 2 + helicities.B(P, M1)))
    assert (
        s.factor(
            averaged_quadratic(kind)["Weyl_squared"]
            - compact
            - clock_derivative(boundary_term)
        )
        == 0
    )


@pytest.mark.parametrize("module", (helicities, local, proper, bounds))
def test_all_exact_core_packets(module):
    d = module.data()
    for value in d["checks"].values():
        entries = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.cancel(x) == 0 for x in entries)
    assert all(d["gates"].values())


TIME = s.Symbol("time", real=True)
SCALE = (1 + TIME * TIME) ** 2
HUBBLE = s.diff(SCALE, TIME) / SCALE
PROFILE = (1 - 4 * TIME * TIME) ** 6


def actual_clock_substitution(profile, qvalue):
    return {
        proper.a: SCALE,
        proper.H: HUBBLE,
        proper.H1: s.diff(HUBBLE, TIME),
        proper.H2: s.diff(HUBBLE, TIME, 2),
        proper.H3: s.diff(HUBBLE, TIME, 3),
        proper.q: s.sympify(qvalue),
        **{proper.X[j]: s.diff(profile, TIME, j) for j in range(7)},
    }


@pytest.mark.parametrize("kind", ("tensor", "vector", "scalar"))
@pytest.mark.parametrize("curvature", ("R_old", "R_old_squared", "Weyl_squared"))
def test_actual_CD_clock_measure_conversion_and_wrong_measure_control(kind, curvature):
    c, b, d, e = proper.HELICITY[kind]
    qvalue = s.Integer(7)
    de = lambda value: s.expand(SCALE * s.diff(value, TIME))
    r = de(de(SCALE)) / SCALE
    x = PROFILE * (1 + TIME)
    conformal = {
        "R_old": -(de(SCALE * SCALE * de(x)) + SCALE * SCALE * c * qvalue * x) / 2,
        "R_old_squared": -6 * (de(r * de(x)) + r * c * qvalue * x)
        + 2 * e * qvalue * qvalue * x,
        "Weyl_squared": de(de(de(de(x))))
        + b * qvalue * de(de(x))
        + d * qvalue * qvalue * x,
    }[curvature]
    actual = proper.operators(kind)[curvature].subs(
        actual_clock_substitution(x, qvalue), simultaneous=True
    )
    assert s.factor(actual - conformal / SCALE) == 0
    assert s.factor((actual - conformal).subs(TIME, s.Rational(1, 4))) != 0


@cache
def numeric_fixed_operator(kind, detector=False):
    profile = (
        PROFILE * (1 + TIME) if detector else PROFILE * (1 - 2 * TIME + 3 * TIME * TIME)
    )
    ops = proper.operators(kind)
    expr = (
        s.Rational(5, 3) * 1000**2 * ops["R_old"]
        - ops["Weyl_squared"] / 30
        - ops["R_old_squared"] / 18
    ) / (64 * s.pi**2)
    raw = s.factor(
        expr.subs(actual_clock_substitution(profile, proper.q), simultaneous=True)
    )
    return s.lambdify((TIME, proper.q), raw, "mpmath", cse=True), s.lambdify(
        TIME, profile, "mpmath"
    )


@pytest.mark.parametrize("kind", ("tensor", "vector", "scalar"))
@pytest.mark.parametrize("qvalue", ("0", "7", "1e12"))
def test_full_fixed_local_proper_time_self_adjoint_pairing(kind, qvalue):
    with mp.workdps(80):
        qvalue = mp.mpf(qvalue)
        HG, G = numeric_fixed_operator(kind)
        HD, D = numeric_fixed_operator(kind, True)
        lhs = mp.quad(
            lambda t: D(t) * HG(t, qvalue), [mp.mpf("-0.5"), 0, mp.mpf("0.5")]
        )
        rhs = mp.quad(
            lambda t: HD(t, qvalue) * G(t), [mp.mpf("-0.5"), 0, mp.mpf("0.5")]
        )
        assert abs(lhs - rhs) < mp.mpf("1e-60") * max(mp.mpf(1), abs(lhs), abs(rhs))


@cache
def derivative_profile_functions():
    profile = PROFILE * (1 - 2 * TIME + 3 * TIME * TIME)
    return tuple(s.lambdify(TIME, s.diff(profile, TIME, j), "mpmath") for j in range(5))


@pytest.mark.parametrize("kind", ("tensor", "vector", "scalar"))
@pytest.mark.parametrize("qvalue", ("0", "7", "1e12"))
def test_complete_local_Hessian_L2_bound_with_actual_derivative_norm(kind, qvalue):
    with mp.workdps(80):
        qvalue = mp.mpf(qvalue)
        HG, _ = numeric_fixed_operator(kind)
        actual = mp.sqrt(
            mp.quad(lambda t: HG(t, qvalue) ** 2, [mp.mpf("-0.5"), 0, mp.mpf("0.5")])
        )
        square = mp.fsum(
            mp.quad(lambda t, f=f: f(t) ** 2, [mp.mpf("-0.5"), 0, mp.mpf("0.5")])
            for f in derivative_profile_functions()
        )
        Phi = (1 + qvalue) ** 2 * mp.sqrt(square)
        assert 0 < actual < mp.mpf(50000) * Phi
        assert 4 * actual / mp.mpf("1e800") < mp.mpf("2e-795") * Phi


@pytest.mark.parametrize("axis", ((1, 2, 3), (2, -1, 1)))
@pytest.mark.parametrize("momentum", ((0, 0, 0), (0, 0, 2), (2, -1, 3)))
def test_full_polynomial_Hessians_under_nontrivial_spatial_rotations(axis, momentum):
    v = s.Matrix(axis) / 7
    skew = s.Matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation = (s.eye(3) - skew).inv() * (s.eye(3) + skew)
    assert rotation.T * rotation == s.eye(3)
    D = s.Matrix([[1, 2, -1], [2, -3, 1], [-1, 1, 2]])
    G = s.Matrix([[2, -1, 3], [-1, 1, 2], [3, 2, -3]])
    jets = [(j + 1) * D + j * j * G for j in range(5)]
    P = s.Matrix(momentum)
    values = (s.Rational(5, 4), s.Rational(1, 3), s.Rational(-1, 2), s.Rational(2, 7))
    original = local.operators(P, jets, *values)
    transformed = local.operators(
        rotation * P, [rotation * M * rotation.T for M in jets], *values
    )
    for name in original:
        assert (transformed[name] - rotation * original[name] * rotation.T).applyfunc(
            s.factor
        ) == s.zeros(3)


@pytest.mark.parametrize("time", tuple(s.Rational(j, 16) for j in range(-8, 9)))
def test_actual_CD_clock_bounds_and_all_helicity_coefficient_sums(time):
    sub = {TIME: time}
    assert 1 <= SCALE.subs(sub) <= s.Rational(25, 16)
    assert abs(HUBBLE.subs(sub)) <= 2
    assert abs(s.diff(HUBBLE, TIME).subs(sub)) <= 4
    assert abs(s.diff(HUBBLE, TIME, 2).subs(sub)) <= 12
    assert bounds.constants()["fixed_local_Hessian_norm_upper"] < 50000


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_every_exact_residual(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(s.cancel(x) == 0 for x in entries)


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_each_invalid_scope_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_complete_counts_and_true_gates():
    assert len(audit.residuals()) == 37 and audit.scalar_entry_count() == 177
    assert len(audit.gates()) == 36 and all(
        value is True for value in audit.gates().values()
    )
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 126


def test_original_primitive_and_matching_statuses_retained():
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert "NOT_QUANTUM_MATCHING_MIXED_INVERSE_OR_V_G_B" in audit.ITEM["status"]
    audit.validate_scope(audit.frontier(), audit.matching())
