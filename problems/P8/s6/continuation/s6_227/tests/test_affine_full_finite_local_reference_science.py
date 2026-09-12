"""Independent full local Hessian, adjoint reciprocity and ordered inverse tests."""

import pytest
import sympy as s
from p8_affine import verify as serializer
from p8_vacuum_affine_full_finite_local_reference import (
    audit,
    coefficients,
    factorization,
    resolvent,
)


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_integrated_exact_residual(name, value):
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(s.factor(v) == 0 for v in entries), name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_integrated_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_unsupported_input_and_false_completion_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_scoped_stages_and_unchanged_frontier():
    for stage in (
        "actual_full_finite_local_factorization",
        "source_time_adjoint_bound",
        "full_finite_reference_inverse",
    ):
        assert audit.require_stage(stage) == stage
    for t in (-s.Rational(1, 2), 0, s.Rational(1, 2)):
        assert audit.require_scope(t)[0] == t
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.validate_scope(audit.frontier(), audit.matching())


def test_independent_complete_local_Euler_factorization_boundary_and_actual_bounds():
    import sympy as s

    t = s.Symbol("conformal_time", real=True)
    q = s.Symbol("squared_transfer", nonnegative=True)
    alpha = s.Function("alpha")(t)
    gamma = s.Function("gamma")(t)
    h = s.Function("h")(t)
    wd, cd, wg, cg = [s.Function(name)(t) for name in ("wD", "cD", "wG", "cG")]
    YD, YG = s.Matrix([wd, cd]), s.Matrix([wg, cg])
    D = lambda value: s.diff(value, t)
    L = s.Matrix([[1, -s.Rational(1, 3)], [0, -1]])
    Q = q * s.Matrix([[s.Rational(2, 3), 0], [-1, 0]])
    C = s.Matrix([[3 * h, -h], [0, 0]])

    def B(Y):
        return L * Y.diff(t, 2) + Q * Y + C * Y.diff(t)

    def Bstar(Y):
        return L.T * Y.diff(t, 2) + Q.T * Y - (C.T * Y).diff(t)

    SD, _WD = B(YD)
    SG, _WG = B(YG)
    TD, TG = 3 * wd - cd, 3 * wg - cg
    rDG = (
        24 * D(wd) * D(wg)
        - 8 * (D(wd) * D(cg) + D(cd) * D(wg))
        + 4 * D(cd) * D(cg)
        - 20 * q * wd * wg
        + 4 * q * (cd * wg + cg * wd)
    )
    remainder = alpha * rDG + 6 * alpha * (TD * SG + TG * SD) + gamma * TD * TG
    observed = s.Matrix(
        [
            sum(
                (-1) ** j * s.diff(s.diff(remainder, s.diff(field, t, j)), t, j)
                for j in range(3)
            )
            for field in YD
        ]
    )
    A = alpha * s.Matrix([[18, -6], [10, -4]])
    H = alpha * s.Matrix([[24, -18], [-18, 12]])
    J = -D(alpha) * s.Matrix([[0, 0], [10, -4]])
    tv = s.Matrix([3, -1])
    V = gamma * tv * tv.T
    predicted = (
        Bstar(A * YG)
        + A.T * B(YG)
        - (H * YG.diff(t)).diff(t)
        - (J * YG).diff(t)
        + J.T * YG.diff(t)
        + V * YG
    )
    assert (observed - predicted).applyfunc(s.expand) == s.zeros(2, 1)
    bilinear = (
        B(YD).T * A * YG
        + YD.T * A.T * B(YG)
        + YD.diff(t).T * H * YG.diff(t)
        + YD.diff(t).T * J * YG
        + YD.T * J.T * YG.diff(t)
        + YD.T * V * YG
    )[0]
    boundary = alpha * (D(cd) * (10 * wg - 4 * cg) + (10 * wd - 4 * cd) * D(cg))
    assert s.expand(remainder - bilinear - D(boundary)) == 0
    wrong = predicted + (J * YG).diff(t) - J.T * YG.diff(t)
    assert (observed - wrong).applyfunc(s.expand) != s.zeros(2, 1)
    assert (observed - predicted).diff(q).applyfunc(s.expand) == s.zeros(2, 1)
    m = s.Integer(1000)
    u = s.Symbol("proper_time", real=True)
    a = (1 + u * u) ** 2
    U = 4 * (1 + u * u) ** 2 * (1 + 7 * u * u)
    Up = s.factor(a * s.diff(U, u))
    assert s.factor(Up - 24 * u * (1 + u * u) ** 3 * (3 + 7 * u * u)) == 0
    assert Up.subs(u, s.Rational(1, 2)) == s.Rational(7125, 64) < 112
    assert s.Rational(5, 3) * m * m * s.Rational(625, 256) < 5 * m * m
    assert 5 * m * m - 2 * 18 > 0
    assert 25 * m * m + s.Rational(224, 3) < 26 * m * m
    assert 2 * 18**2 + 10 * m * m * 3 * 18 + s.Rational(5, 2) * m**4 * 6 < 16 * m**4
    assert sum(v * v for v in A / alpha) == 476 < 22**2
    assert sum(v * v for v in H / alpha) == 1368 < 37**2
    assert sum(v * v for v in J / D(alpha)) == 116 < 11**2
    assert tv.dot(tv) == 10
    C0 = s.Rational(5, 2)
    assert 2 * 110 * C0 + 185 * C0 * C0 == s.Rational(6825, 4)
    assert 2 * 286 * C0 * C0 == 3575
    assert 160 * C0 * C0 == 1000
    upper = (
        s.Rational(6825, 4) / 100**2
        + s.Rational(3575, 100**3) / m
        + s.Rational(1000, 100**4)
    )
    assert upper < s.Rational(1, 4)


def test_independent_advanced_reciprocity_source_derivative_and_initial_atom():
    import sympy as s

    # Exactly solvable zero-transfer variable-coefficient B with h constant:
    # diagonalize the fixed leading matrix, but keep the first-order nonsymmetric C.
    t, u = s.symbols("time source_time", real=True)
    h = s.Rational(2, 5)
    L = s.Matrix([[1, -s.Rational(1, 3)], [0, -1]])
    C = s.Matrix([[3 * h, -h], [0, 0]])
    N = L.inv() * C
    assert N == s.Matrix([[3 * h, -h], [0, 0]])
    I = s.eye(2)

    def velocity(x, sign):
        return s.exp(-sign * N * x) * L.inv()

    def kernel(x, sign):
        v = s.Symbol("integration_time", real=True)
        return s.integrate(velocity(v, sign), (v, 0, x))

    Y = kernel(t - u, 1)
    # The reversed-time B has -C, and its retarded kernel gives Y_adv(u,t).
    Z = kernel(t - u, -1).T
    assert (L * Y.diff(t, 2) + C * Y.diff(t)).applyfunc(s.simplify) == s.zeros(2)
    assert (L.T * Z.diff(t, 2) - C.T * Z.diff(t)).applyfunc(s.simplify) == s.zeros(2)
    assert Y.subs(t, u) == s.zeros(2) and Z.subs(t, u) == s.zeros(2)
    assert (L * Y.diff(t).subs(t, u)).applyfunc(s.simplify) == I
    assert (L.T * Z.diff(t).subs(t, u)).applyfunc(s.simplify) == I
    assert (Z - Y.T).applyfunc(s.simplify) != s.zeros(2)
    # The source-time derivative, not a naively transposed same-time primal,
    # represents ZD*. Include the original derivative's initial delta explicitly.
    f = s.Matrix([1 + u + u * u, s.Rational(2, 3) - u + u**3])
    ordinary = -s.integrate(Z * f.diff(u), (u, 0, t)) - Z.subs(u, 0) * f.subs(u, 0)
    by_source = s.integrate(Z.diff(u) * f, (u, 0, t))
    assert (ordinary - by_source).applyfunc(s.simplify) == s.zeros(2, 1)
    omitted = -s.integrate(Z * f.diff(u), (u, 0, t))
    assert (omitted - by_source).applyfunc(s.simplify) != s.zeros(2, 1)
    # Independently nonconstant h on a finite polynomial Taylor jet:
    # causal Green kernel normalized at source u is solved by recurrences.
    x = s.Symbol("elapsed", real=True)
    hvar = 1 + u + x + (u + x) ** 2 / 3
    Cx = s.Matrix([[3 * hvar, -hvar], [0, 0]])
    Cp = Cx.diff(x)
    q = s.Rational(7, 5)
    Q = q * s.Matrix([[s.Rational(2, 3), 0], [-1, 0]])
    order = 6

    def jet_inverse(adjoint):
        lead = L.T if adjoint else L
        coeff = {0: s.zeros(2), 1: lead.inv()}
        for n in range(order - 1):
            poly = sum((coeff[j] * x**j for j in coeff), s.zeros(2))
            if adjoint:
                rest = Q.T * poly - Cx.T * poly.diff(x) - Cp.T * poly
            else:
                rest = Q * poly + Cx * poly.diff(x)
            coefficient = rest.applyfunc(lambda v, n=n: s.expand(v).coeff(x, n))
            coeff[n + 2] = -lead.inv() * coefficient / ((n + 2) * (n + 1))
        return sum((coeff[j] * x**j for j in coeff), s.zeros(2))

    Yjet = jet_inverse(False)
    Zjet = jet_inverse(True)
    # For a general smooth fundamental kernel, the advanced B kernel in the
    # forward signed displacement is the negative analytic continuation ofYret.
    # Therefore Zret(u+x,u)=-Yret(u,u+x)^T. Expand the exchanged source point.
    exchanged = -Yjet.subs({x: -x, u: u + x}, simultaneous=True).T
    diff = (Zjet - exchanged).applyfunc(s.expand)
    for n in range(order + 1):
        assert diff.applyfunc(lambda v, n=n: v.coeff(x, n)) == s.zeros(2), n


@pytest.mark.parametrize("q", (s.Rational(0), s.Rational(7, 5), s.Rational(10000)))
def test_independent_full_noncommuting_causal_block_inverse(q):
    # Finite operator algebra diagnostics, not a continuum discretization.
    I, O = s.eye(2), s.zeros(2)
    D = s.Matrix([[1, 0], [-1, 1]])
    Db = s.diag(D, D)
    D2 = D * D
    h = s.diag(s.Rational(1, 5), s.Rational(2, 7))
    C = s.BlockMatrix([[3 * h, -h], [O, O]]).as_explicit()
    CT = s.BlockMatrix([[3 * h, O], [-h, O]]).as_explicit()
    B0 = s.BlockMatrix([[D2 + 2 * q * I / 3, -D2 / 3], [-q * I, -D2]]).as_explicit()
    B0T = s.BlockMatrix([[D2 + 2 * q * I / 3, -q * I], [-D2 / 3, -D2]]).as_explicit()
    B, Bs = B0 + C * Db, B0T - Db * CT
    Y, Z = B.inv(), Bs.inv()
    a = s.diag(s.Rational(1, 10**5), s.Rational(3, 10**5))
    ap = s.diag(s.Rational(1, 10**6), -s.Rational(2, 10**6))
    gamma = s.diag(s.Rational(1, 10**7), s.Rational(3, 10**7))
    A = s.BlockMatrix([[18 * a, -6 * a], [10 * a, -4 * a]]).as_explicit()
    H = s.BlockMatrix([[24 * a, -18 * a], [-18 * a, 12 * a]]).as_explicit()
    J = s.BlockMatrix([[O, O], [-10 * ap, 4 * ap]]).as_explicit()
    V = s.BlockMatrix([[9 * gamma, -3 * gamma], [-3 * gamma, gamma]]).as_explicit()
    R = Bs * A + A.T * B - Db * H * Db - Db * J + J.T * Db + V
    omega = Z * R * Y
    decomposed = (
        A * Y
        + Z * A.T
        - Z * Db * H * Db * Y
        - Z * Db * J * Y
        + Z * J.T * Db * Y
        + Z * V * Y
    )
    assert omega == decomposed
    assert omega != A * Y + Z * A.T - Z * Db * H * Db * Y + Z * V * Y
    F = s.diag(-4 * I - D2 / 5, s.Rational(8, 3) * (-I / 30 - D2 / 7))
    K = F.inv()
    E = Y * (s.eye(4) + K * omega).inv() * K * Z
    full = Bs * F * B + R
    assert full * E == s.eye(4) and E * full == s.eye(4)
    assert E != Y * K * Z
    assert K * omega != omega * K
    assert E != Y * K * (s.eye(4) + K * omega).inv() * Z


@pytest.mark.parametrize(
    "bound", (0, s.Rational(1, 10**12), s.Rational(3, 2), 10, 1000)
)
def test_independent_actual_local_perturbation_and_total_inverse_constants(bound):
    bound = s.Rational(bound)
    mass = s.Integer(1000)
    omega = (
        s.Rational(6825, 4) / 100**2
        + s.Rational(3575, 100**3) / mass
        + s.Rational(1000, 100**4)
    )
    assert omega == s.Rational(6825543, 40000000) < s.Rational(1, 4)
    gap = (32 + 13 * bound) / 5
    assert (bound + s.Rational(1, 4)) / gap < s.Rational(5, 13) < s.Rational(1, 2)
    assert 1 / (gap - bound - s.Rational(1, 4)) == 20 / (123 + 32 * bound)
    assert s.Rational(5, 2) ** 2 * 20 == 125
    assert 384 * 125 == 48000
    # Elementary exponential lower polynomial controls the entire M>=0 class.
    assert 2 * mass * (1 + 32 + s.Rational(32**2, 2)) - 20 > 100 * mass


@pytest.mark.parametrize("module", (factorization, coefficients, resolvent))
def test_exact_module_checks_gates_and_serialization(module):
    data = module.data()
    for name, value in data["checks"].items():
        values = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.factor(v) == 0 for v in values), name
    assert all(bool(value) for value in data["gates"].values())
    serializer.serialize(
        {key: value for key, value in data.items() if key not in ("checks", "gates")}
    )
