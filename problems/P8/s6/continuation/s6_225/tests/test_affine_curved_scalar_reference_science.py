"""Independent curved geometry, local action and noncommuting inverse checks."""

from functools import cache

import pytest
import sympy as s
from p8_affine import verify as serializer
from p8_vacuum_affine_curved_scalar_reference import audit, local, resolvent
from p8_vacuum_affine_curved_scalar_reference import geometry as g


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


def test_scoped_stages_and_unchanged_frontier():
    for stage in (
        "actual_curvature_coordinates",
        "complete_finite_local_Hessian",
        "ordered_curved_reference_inverse",
    ):
        assert audit.require_stage(stage) == stage
    for t in (-s.Rational(1, 2), 0, s.Rational(1, 2)):
        assert audit.require_scope(t)[0] == t
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.validate_scope(audit.frontier(), audit.matching())


@cache
def independent_linear_four_geometry():
    t, x, y, z = s.symbols("independent_eta x y z", real=True)
    coords = (t, x, y, z)
    eps = s.Symbol("independent_variation")
    k = s.Symbol("independent_transfer", positive=True)
    phase = s.exp(s.I * k * x)
    a = s.Function("independent_scale", positive=True)(t)
    w = s.Function("independent_w")(t)
    c = s.Function("independent_c")(t)
    metric = s.diag(a * a, *[-a * a * (1 + 2 * eps * phase * v) for v in (w - c, w, w)])
    inverse = s.diag(
        1 / (a * a), *[-(1 - 2 * eps * phase * v) / (a * a) for v in (w - c, w, w)]
    )

    def lin(value):
        value = s.expand(value)
        return value.coeff(eps, 0) + eps * value.coeff(eps, 1)

    Gamma = {}
    for aa in range(4):
        for bb in range(4):
            for cc in range(4):
                Gamma[aa, bb, cc] = lin(
                    inverse[aa, aa]
                    * (
                        s.diff(metric[aa, cc], coords[bb])
                        + s.diff(metric[aa, bb], coords[cc])
                        - s.diff(metric[bb, cc], coords[aa])
                    )
                    / 2
                )
    Ric = s.Matrix(
        4,
        4,
        lambda bb, dd: lin(
            sum(
                s.diff(Gamma[aa, dd, bb], coords[aa])
                - s.diff(Gamma[aa, aa, bb], coords[dd])
                + sum(
                    Gamma[aa, aa, rr] * Gamma[rr, dd, bb]
                    - Gamma[aa, dd, rr] * Gamma[rr, aa, bb]
                    for rr in range(4)
                )
                for aa in range(4)
            )
        ),
    )
    R = lin(sum(inverse[aa, aa] * Ric[aa, aa] for aa in range(4)))
    Weyl = {}
    for aa in range(4):
        for bb in range(4):
            for cc in range(4):
                for dd in range(4):
                    riem = lin(
                        metric[aa, aa]
                        * (
                            s.diff(Gamma[aa, dd, bb], coords[cc])
                            - s.diff(Gamma[aa, cc, bb], coords[dd])
                            + sum(
                                Gamma[aa, cc, rr] * Gamma[rr, dd, bb]
                                - Gamma[aa, dd, rr] * Gamma[rr, cc, bb]
                                for rr in range(4)
                            )
                        )
                    )
                    Weyl[aa, bb, cc, dd] = lin(
                        riem
                        - (
                            metric[aa, cc] * Ric[bb, dd]
                            - metric[aa, dd] * Ric[bb, cc]
                            - metric[bb, cc] * Ric[aa, dd]
                            + metric[bb, dd] * Ric[aa, cc]
                        )
                        / 2
                        + R
                        * (
                            metric[aa, cc] * metric[bb, dd]
                            - metric[aa, dd] * metric[bb, cc]
                        )
                        / 6
                    )
    return t, eps, k, phase, a, w, c, R, Weyl


def test_independent_first_curvature_full_256_Weyl_components_and_density():
    t, eps, k, phase, a, w, c, R, Weyl = independent_linear_four_geometry()
    h = s.diff(a, t) / a
    S = (
        s.diff(w, t, 2)
        + 3 * h * s.diff(w, t)
        + 2 * k * k * w / 3
        - s.diff(c, t, 2) / 3
        - h * s.diff(c, t)
    )
    W = -k * k * w - s.diff(c, t, 2)
    assert s.simplify(-R.coeff(eps, 0) - 6 * s.diff(a, t, 2) / a**3) == 0
    assert s.simplify(-R.coeff(eps, 1) * a * a / (6 * phase) - S) == 0
    contraction = 0
    nonzero = 0
    signatures = (1, -1, -1, -1)
    for indices, value in Weyl.items():
        assert s.simplify(value.coeff(eps, 0)) == 0
        flat = s.simplify(value.coeff(eps, 1).subs(a, 1).doit())
        assert s.simplify(value.coeff(eps, 1) / a**2 - flat) == 0
        nonzero += flat != 0
        contraction += s.prod(signatures[i] for i in indices) * flat**2 / phase**2
    assert nonzero == 24
    assert s.simplify(Weyl[0, 1, 0, 1].coeff(eps, 1) / (a * a * phase) - W / 3) == 0
    assert s.simplify(contraction - s.Rational(4, 3) * W**2) == 0
    # The mixed Hessian has two opposite legs, so the original -1/30
    # coefficient yields -4/45 WD WG, not -2/45 or a fitted counterterm.
    assert -s.Rational(1, 30) * 2 * s.Rational(4, 3) == -s.Rational(4, 45)


def test_independent_ADM_warped_spatial_curvature_matches_literal_four_geometry():
    # Derive R from the ADM extrinsic-curvature and warped-three-metric formula,
    # independently of production's four-dimensional Christoffel contraction.
    w = g.ed * g.wd * g.pd + g.eg * g.wg * g.pg
    c = g.ed * g.cd * g.pd + g.eg * g.cg * g.pg
    D = lambda value: s.diff(value, g.eta)
    Dx = lambda value: s.diff(value, g.x)
    spatial = g.exponential(-2 * (w - c)) * (
        -4 * Dx(Dx(w)) - 2 * Dx(w) ** 2 - 4 * Dx(c) * Dx(w)
    )
    observed = g.jet(
        6 * g.U
        + 6 * D(D(w))
        + 18 * g.h * D(w)
        - 2 * D(D(c))
        - 6 * g.h * D(c)
        + 12 * D(w) ** 2
        - 8 * D(w) * D(c)
        + 2 * D(c) ** 2
        + spatial
    )
    assert s.simplify(observed - g.a**2 * g.literal()["Rold"]) == 0
    volume = g.a**4 * g.exponential(3 * w - c)
    R = observed / g.a**2
    m = s.Integer(1000)
    fixed = g.jet(volume * (5 * m**4 / 2 + 5 * m * m * R / 3 - R * R / 18))
    mixed = s.expand(fixed).coeff(g.ed, 1).coeff(g.eg, 1)
    WD, WG = -g.q * g.wd - D(D(g.cd)), -g.q * g.wg - D(D(g.cg))
    mixed -= s.Rational(4, 45) * WD * WG
    packet = local.data()
    assert (
        s.simplify(
            mixed
            - packet["curvature_square_factor"]
            - packet["complete_lower_order_local_remainder"]
        )
        == 0
    )


def test_independent_actual_clock_and_continuous_interval_bounds():
    t = s.Symbol("independent_proper_time", real=True)
    a = (1 + t * t) ** 2
    h = 4 * t * (1 + t * t)
    hp = s.factor(a * s.diff(h, t))
    U = s.factor(hp + h * h)
    endpoint = s.Rational(1, 2)
    # Positive coefficients and parity give whole-interval bounds, not samples.
    assert all(v >= 0 for v in s.Poly(h, t).all_coeffs())
    assert all(v >= 0 for v in s.Poly(hp, t).all_coeffs())
    assert all(v >= 0 for v in s.Poly(U, t).all_coeffs())
    assert h.subs(t, endpoint) == s.Rational(5, 2)
    assert hp.subs(t, endpoint) == s.Rational(175, 16)
    assert U.subs(t, endpoint) == s.Rational(275, 16)
    assert (a**4).subs(t, endpoint) == s.Rational(390625, 65536) < 6
    primitive = t / (2 * (1 + t * t)) + s.atan(t) / 2
    assert s.factor(s.diff(primitive, t) - 1 / a) == 0
    assert s.Poly(a - 1, t).all_coeffs() == [1, 0, 2, 0, 0]
    # a >= 1 implies the conformal interval length <= its proper length 1.
    assert 10 * h.subs(t, endpoint) ** 2 < 8**2
    assert 10 * hp.subs(t, endpoint) ** 2 < 35**2


@pytest.mark.parametrize("q", (s.Rational(0), s.Rational(7, 5), s.Rational(10000)))
def test_independent_complete_variable_coefficient_adjoint_and_deleted_term(q):
    t = s.Symbol("test_time", real=True)
    h = 1 + t + t * t / 3
    C = s.Matrix([[3 * h, -h], [0, 0]])
    Y = s.Matrix([t**3 * (1 - t) ** 3, t**4 * (1 - t) ** 4])
    Z = s.Matrix([2 * t**4 * (1 - t) ** 4, -(t**3) * (1 - t) ** 3 / 3])
    L = s.Matrix([[1, -s.Rational(1, 3)], [0, -1]])
    Q = q * s.Matrix([[s.Rational(2, 3), 0], [-1, 0]])
    primal = L * Y.diff(t, 2) + Q * Y + C * Y.diff(t)
    adjoint = L.T * Z.diff(t, 2) + Q.T * Z - (C.T * Z).diff(t)
    assert s.integrate(s.expand((Z.T * primal - adjoint.T * Y)[0]), (t, 0, 1)) == 0
    wrong = L.T * Z.diff(t, 2) + Q.T * Z - C.T * Z.diff(t)
    assert s.integrate(s.expand((Z.T * primal - wrong.T * Y)[0]), (t, 0, 1)) != 0


@pytest.mark.parametrize("q", (s.Rational(0), s.Rational(7, 5), s.Rational(10000)))
def test_independent_ordered_finite_causal_algebra_not_a_continuum_certificate(q):
    I, O = s.eye(2), s.zeros(2)
    D = s.Matrix([[1, 0], [-1, 1]])
    D2 = D * D
    H = s.diag(s.Rational(1, 5), s.Rational(2, 7))
    C = s.BlockMatrix([[3 * H, -H], [O, O]]).as_explicit()
    CT = s.BlockMatrix([[3 * H, O], [-H, O]]).as_explicit()
    Db = s.diag(D, D)
    B = s.BlockMatrix([[D2 + 2 * q * I / 3, -D2 / 3], [-q * I, -D2]]).as_explicit()
    BT = s.BlockMatrix([[D2 + 2 * q * I / 3, -q * I], [-D2 / 3, -D2]]).as_explicit()
    Bc, Bs = B + C * Db, BT - Db * CT
    F = s.diag(-4 * I - D2 / 5, s.Rational(8, 3) * (-I / 30 - D2 / 7))
    K = F.inv()
    J = Db.inv() * K
    Y, Z = Bc.inv(), Bs.inv()
    E = Y * J * Db * Z
    A = Bs * F * Bc
    assert A * E == s.eye(4) and E * A == s.eye(4)
    assert E != Z * J * Db * Y
    assert E != Y * J * Db * (BT - CT * Db).inv()
    density = s.diag(1, 2, 1, 2)
    force = density.inv() * A
    assert force * (E * density) == s.eye(4)
    assert (E * density) * force == s.eye(4)
    assert force * (density * E) != s.eye(4)


def test_independent_Volterra_and_ordered_kernel_constants():
    C0, J = s.Rational(5, 2), s.Rational(45, 2)
    assert C0 * 8 == 20
    assert C0 * (8 + 35) == s.Rational(215, 2) < 108
    u, v, T = s.symbols("u v window", positive=True)
    simplex = s.integrate(s.integrate(T - u, (v, 0, u)), (u, 0, T))
    assert simplex == T**3 / 6
    assert C0 * C0 * J / 6 == s.Rational(375, 16)
    assert C0 * C0 * J / 24 == s.Rational(375, 64) < 6
    assert 6 * 64 * C0 * C0 * J / 24 == 2250
    # No numerical exp(108T) evaluation or kappa cancellation enters this bound.


@pytest.mark.parametrize("module", (g, local, resolvent))
def test_exact_module_checks_gates_and_serialization(module):
    data = module.data()
    for key, value in data["checks"].items():
        values = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.factor(v) == 0 for v in values), key
    assert all(bool(value) for value in data["gates"].values())
    serializer.serialize(
        {key: value for key, value in data.items() if key not in ("checks", "gates")}
    )
