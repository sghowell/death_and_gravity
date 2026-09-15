"""Independent exact canonical, preservation, endpoint and scope diagnostics."""

import itertools

import pytest
import sympy as s
from p8_vacuum_affine_homogeneous_fold_dynamics import (
    audit,
    canonical,
    consistency,
    source,
)

ROWS = audit.residuals()
GATES = audit.gates()
BAD = audit.bad_cases()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_all_exact_residuals(name):
    value = ROWS[name]
    if isinstance(value, s.MatrixBase):
        assert all(x == 0 for x in value)
    else:
        assert value == 0
    assert not value.atoms(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_written_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize("name,call,args", BAD, ids=[row[0] for row in BAD])
def test_all_unsupported_inputs_rejected(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "kind,pair",
    [(kind, pair) for kind, pairs in source.ALLOWED.items() for pair in pairs],
)
def test_clock_jet_validates_before_warm_cache(kind, pair):
    expected = source.clock_jet(kind, *pair)
    assert source.clock_jet(kind, *map(s.Integer, pair)) == expected
    for index, value in enumerate(pair):
        for bad in (bool(value), float(value), s.Float(value)):
            modified = list(pair)
            modified[index] = bad
            with pytest.raises((TypeError, ValueError)):
                source.clock_jet(kind, *modified)


@pytest.mark.parametrize("order", range(3))
def test_full_homogeneous_jets_match_independent_S276_source_at_zero_trace(order):
    old = source.previous.constraint_jet(order)
    translated = old.subs(
        {
            source.previous.RHO: source.RHO,
            source.previous.PRESSURE: source.PRESSURE,
            source.previous.M: source.m,
            source.previous.CURVATURE: 0,
            source.previous.SHEAR: source.sh,
        },
        simultaneous=True,
    )
    current = source.clock_jet("C", 0, order).subs(source.ZERO_HEAVY).subs(source.p, 0)
    assert s.factor(current - translated) == 0


def test_solve_both_full_clock_constraints_for_fold_family():
    M2 = s.Symbol("independent_squared_matter")
    equations = [
        source.clock_jet("C", 0, k).subs(source.ZERO_HEAVY).subs(source.m**2, M2)
        for k in (0, 1)
    ]
    solution = s.solve(equations, (source.sh, M2), dict=True)
    assert len(solution) == 1
    assert s.factor(solution[0][source.sh] - source.FOLD_SHEAR) == 0
    assert s.factor(solution[0][M2] - source.FOLD_MATTER_SQUARE) == 0


@pytest.mark.parametrize(
    "pvalue",
    (
        s.Rational(-3, 7),
        s.Rational(-1, 10),
        s.Integer(0),
        s.Rational(1, 10),
        s.Rational(5, 4),
    ),
)
def test_complete_fold_and_both_M1_signs(pvalue):
    r, P = source.RHO, source.PRESSURE
    eps = source.PROFILE_BOUND
    bind = {r: eps / 7, P: -eps / 11, source.p: pvalue}
    shear = source.FOLD_SHEAR.subs(bind)
    mass2 = source.FOLD_MATTER_SQUARE.subs(bind)
    assert shear > s.Rational(49, 100) and mass2 > 6
    for sign in (-1, 1):
        rule = {
            **bind,
            source.sh: shear,
            source.m: sign * s.sqrt(mass2),
            **source.ZERO_HEAVY,
        }
        assert all(
            s.factor(source.clock_jet("C", 0, k).subs(rule)) == 0 for k in (0, 1)
        )


@pytest.mark.parametrize(
    "pvalue", (s.Rational(-1, 10), s.Rational(1, 10), s.Rational(3, 7))
)
def test_dropped_time_boundary_is_a_nonzero_error(pvalue):
    full = source.clock_jet("C", 1, 0).subs(source.ZERO_HEAVY)
    omitted = 3 * source.PRESSURE1 / 2 - source.RHO1
    defect = s.factor((full - omitted).subs(source.p, pvalue))
    assert defect == -9 * pvalue and defect != 0
    assert source.clock_jet("B", 1, 0) == 0
    assert source.clock_jet("B", 1, 1) == -6


def literal_bracket(A, G):
    alpha, h, Pa, PM, PB, PH = s.symbols(
        "test_alpha test_heavy test_Palpha test_PM test_Pbeta test_PH", real=True
    )
    kap = s.Integer(3)
    V = s.exp(3 * alpha)
    bind = {
        source.p: Pa / (3 * kap * V),
        source.m: PM / (kap * V),
        source.sh: PB**2 / (8 * kap**2 * V**2),
        source.eta: 10**100 * h,
        source.ph: PH / (kap * V),
    }
    a = A.subs(bind, simultaneous=True)
    H = kap * V * G.subs(bind, simultaneous=True)
    result = sum(
        s.diff(a, q) * s.diff(H, p) - s.diff(a, p) * s.diff(H, q)
        for q, p in ((alpha, Pa), (h, PH))
    )
    return s.factor(
        result.subs(
            {
                alpha: 0,
                Pa: 3 * kap * source.p,
                PM: kap * source.m,
                PB**2: 8 * kap**2 * source.sh,
                PH: kap * source.ph,
                h: source.eta / 10**100,
            },
            simultaneous=True,
        )
    )


@pytest.mark.parametrize(
    "A",
    (*canonical.FIELDS, source.p * source.m + source.eta * source.ph + source.sh**2),
)
def test_independent_literal_canonical_flow(A):
    H = source.clock_jet("H")
    assert s.factor(literal_bracket(A, H) - canonical.flow(A, H)) == 0


def test_second_preservation_from_independent_absolute_brackets():
    H, Hu = source.clock_jet("H"), source.clock_jet("H", 1, 0)
    C, Cu, Cuu = (source.clock_jet("C", k, 0) for k in range(3))
    K = Cu + literal_bracket(C, H)
    second = (
        Cuu + literal_bracket(Cu, H) + literal_bracket(C, Hu) + literal_bracket(K, H)
    )
    assert s.factor(source.at_fold(second) - consistency.coefficients()["c"]) == 0


def test_normalized_density_self_flow_is_not_antisymmetric():
    C = 1 + source.p**2 + source.m * source.sh + source.eta * source.ph
    actual = literal_bracket(C, C)
    assert s.factor(actual + C * s.diff(C, source.p)) == 0
    assert actual != 0


@pytest.mark.parametrize("signs", tuple(itertools.product((-1, 1), repeat=6)))
def test_profile_corner_diagnostics_not_a_substitute_for_enclosure(signs):
    rules = {
        q: sign * source.PROFILE_BOUND
        for q, sign in zip(source.PROFILE_JETS, signs, strict=True)
    }
    pc = s.factor(consistency.COMPATIBLE_P.subs(rules))
    assert abs(pc) < 21 * source.PROFILE_BOUND
    rules[source.p] = pc
    a, b, c = (
        s.factor(consistency.coefficients()[name].subs(rules))
        for name in ("a", "b", "c")
    )
    assert 28 < a < 29 and abs(b) < 1 and c > 60
    assert b * b - 4 * a * c < -6719
    assert s.factor(consistency.FIRST.subs(rules)) == 0


def test_polynomial_enclosure_uses_all_nonconstant_monomials():
    p = source.p
    r, P = source.RHO, source.PRESSURE
    expr = 7 - 5 * r + 3 * P * p - 2 * p**3 + r**2
    e = source.PROFILE_BOUND
    expected = 5 * e + 3 * e * (21 * e) + 2 * (21 * e) ** 3 + e**2
    assert consistency.remainder_bound(expr) == expected


def test_solved_independent_polynomial_fold_dynamics():
    tau, u, N, q, p = s.symbols(
        "test_tau test_time test_lapse test_q test_p", real=True
    )
    H = p + (N - 1) * q + (N - 1) ** 3 / 3
    C = s.diff(H, N)
    path = {u: -(tau**2), N: 1 - tau, q: -(tau**2), p: -2 * tau**3 / 3}
    assert s.factor(C.subs(path)) == 0
    ut = s.diff(path[u], tau)
    assert s.factor(s.diff(path[q], tau) - ut * s.diff(H, p).subs(path)) == 0
    assert s.factor(s.diff(path[p], tau) + ut * s.diff(H, q).subs(path)) == 0
    assert s.factor(s.diff(path[N], tau) / ut - 1 / (2 * tau)) == 0
    assert s.diff(C, N).subs(path) == -2 * tau


def test_first_compatibility_alone_does_not_give_real_time_passage():
    u, N, q, p, lam = s.symbols("test_u test_N test_q test_p test_lambda", real=True)
    H = N**3 / 3 + N * (q + u**2)
    C = s.diff(H, N)
    K = s.diff(C, u) + s.diff(C, q) * s.diff(H, p) - s.diff(C, p) * s.diff(H, q)
    assert C.subs({u: 0, N: 0, q: 0}) == 0 and K.subs({u: 0, N: 0, q: 0}) == 0
    quadratic = s.diff(C, N, 2) * lam**2 + s.diff(K, u)
    assert s.expand(quadratic - 2 * (lam**2 + 1)) == 0
    assert s.discriminant(quadratic, lam) == -16


@pytest.mark.parametrize(
    "D,J", ((s.Integer(5), s.Integer(2)), (s.Integer(29), s.Rational(3, 250)))
)
def test_endpoint_geometry_with_nonzero_higher_Taylor_contacts(D, J):
    tau = s.Symbol("test_endpoint_tau", positive=True)
    u = -D * J * tau**2 / 2 + tau**3 / 17
    N = 1 - J * tau + tau**2 / 19
    R = 1 + 2 * J * tau + tau**2 / 23 + tau**3 / 29
    alpha, beta = tau**2 / 7, tau**2 / 13
    scales = [
        s.exp(alpha + beta) * R ** (-s.Rational(1, 4)),
        s.exp(alpha - beta) * R ** (-s.Rational(1, 4)),
        s.exp(alpha) * R ** (-s.Rational(1, 4)),
    ]
    proper_tau = N * s.diff(u, tau)
    rates = [s.diff(A, tau) / (A * proper_tau) for A in scales]
    theta = s.factor(sum(rates))
    ricci = -2 * s.diff(theta, tau) / proper_tau - theta**2 - sum(x * x for x in rates)
    assert s.limit(tau * theta, tau, 0, dir="+") == 3 / (2 * D)
    assert s.limit(tau**3 * ricci, tau, 0, dir="+") == -3 / (D**2 * J)


def test_original_physical_lapse_in_volume_fixture():
    t = s.Symbol("test_volume_time", real=True)
    N = s.exp(t)
    R = s.exp(-2 * t)
    V = s.exp(3 * t * t) * R ** (-s.Rational(3, 4))
    correct = s.simplify(s.diff(V, t) / (N * V))
    wrong = s.simplify(s.diff(V, t) / (N * R ** (-s.Rational(1, 4)) * V))
    assert s.simplify(correct - (6 * t + s.Rational(3, 2)) * s.exp(-t)) == 0
    assert s.simplify(correct - wrong) != 0


@pytest.mark.parametrize("row", audit.qualifications())
def test_each_historical_physical_qualification_preserved(row):
    assert row in audit.previous.qualifications()


def test_complete_frontiers_and_scope_remain_open():
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 133
    assert len(audit.qualifications()) == 6
    assert audit.controls()["rejected_inputs"] == 70
    assert audit.require_parameters(audit.parameters()) == audit.parameters()
    assert "S275" in audit.observable()["not_established"]
