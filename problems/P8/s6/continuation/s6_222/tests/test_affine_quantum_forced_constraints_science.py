"""Independent full forced saddle, clock contact, ordered block and norm controls."""

import pytest
import sympy as s
from p8_affine import verify as serializer
from p8_vacuum_affine_quantum_forced_constraints import (
    audit,
    estimates,
    feedback,
    forces,
)

PACKETS = (forces.data, forces.contact_data, feedback.data, estimates.data)


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_integrated_exact_residual(name, value):
    rows = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(s.factor(v) == 0 for v in rows), name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_integrated_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_every_false_completion_or_unsupported_input_rejected(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


def test_current_scope_stages_and_all_previous_rows_unchanged():
    for t in (-s.Rational(1, 2), 0, s.Rational(1, 2)):
        assert audit.require_scope(t)[0] == t
    for stage in (
        "force_interface",
        "admissible_graph_identity",
        "first_prescribed_response",
    ):
        assert audit.require_stage(stage) == stage
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.frontier()) == 9
    assert audit.validate_scope(audit.frontier(), audit.matching())


@pytest.mark.parametrize("packet", PACKETS)
def test_exact_science_packet_and_serialization(packet):
    data = packet()
    for key, value in data["checks"].items():
        rows = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.factor(v) == 0 for v in rows), key
    assert all(bool(v) for v in data["gates"].values())
    serializer.serialize(
        {k: v for k, v in data.items() if k not in ("checks", "gates")}
    )


@pytest.mark.parametrize("theta", (-s.Rational(17, 10), 0, s.Rational(19, 10)))
@pytest.mark.parametrize("qval", (s.Rational(1, 10**12), 2, 10000))
def test_independent_forced_action_saddle_and_all_phase_equations(theta, qval):
    JJ = s.Rational(13, 10)
    ell = s.Rational(2, 25)
    w = s.Rational(-31, 1000)
    E = s.Rational(-23, 100)
    A = s.Rational(-11, 1000)
    T = s.Rational(19, 1000)
    delta = s.Rational(1, 3)
    H = s.Rational(-2, 7)
    v0, sg0, pv0, ps0 = map(s.Rational, ("1/7", "-2/9", "3/11", "-5/13"))
    gn, gz, gb = map(s.Rational, ("7/17", "-11/19", "13/23"))
    n, b, vd, sd = s.symbols("lapse shift velocity matter_velocity", real=True)
    # Separately entered scalar action, without the production L or its constraints.
    v, sg = s.symbols("field matter", real=True)
    L = -3 * vd**2 + (JJ + w * w / 2 - 3 * theta**2) * n * n + 6 * theta * n * vd
    L += (
        sd * sd / 2
        + w * n * sd
        - 3 * ell * vd * sg
        + 2 * b * (vd - theta * n)
        + ell * b * sg
    )
    L += (
        qval * v * v
        + 2 * E * qval * n * v
        - qval * sg * sg / 2
        + 3 * T * n * v
        + s.Rational(9, 2) * A * v * v
    )
    eq = s.Matrix(
        [
            s.diff(L, vd) - pv0,
            s.diff(L, sd) - ps0,
            s.diff(L, n) + gn + delta * gz,
            s.diff(L, b) + gb,
        ]
    ).subs({v: v0, sg: sg0})
    variables = s.Matrix([vd, sd, n, b])
    matrix = eq.jacobian(variables)
    constant = eq.subs(dict.fromkeys(variables, 0))
    values = matrix.inv() * (-constant)
    solution = dict(zip(variables, values))
    assignment = {
        forces.th: theta,
        forces.J: JJ,
        forces.l: ell,
        forces.w: w,
        forces.E: E,
        forces.A: A,
        forces.T: T,
        forces.delta: delta,
        forces.H: H,
        forces.q: qval,
        forces.v: v0,
        forces.sigma: sg0,
        forces.pv: pv0,
        forces.ps: ps0,
        forces.gn: gn,
        forces.gz: gz,
        forces.gb: gb,
    }
    predicted = forces.system()
    assert all(
        s.factor(v) == 0
        for v in values
        - s.Matrix(
            [
                predicted["solution"][key].subs(assignment)
                for key in (forces.vd, forces.sd, forces.n, forces.b)
            ]
        )
    )
    pd = (s.diff(L, v) + gz).subs({v: v0, sg: sg0}).subs(solution) - 3 * H * pv0
    psd = s.diff(L, sg).subs({v: v0, sg: sg0}).subs(solution) - 3 * H * ps0
    direct = s.Matrix([solution[vd], solution[sd], pd, psd])
    assert all(s.factor(v) == 0 for v in direct - predicted["flow"].subs(assignment))
    metric = s.Matrix([solution[n], v0 + delta * solution[n], solution[b]])
    assert all(s.factor(v) == 0 for v in metric - predicted["metric"].subs(assignment))
    # Both wrong classical-only constraints are actually detected.
    classical = predicted["metric"].subs(dict.fromkeys(forces.g, 0)).subs(assignment)
    assert s.factor(metric[0] - classical[0]) != 0
    assert s.factor(metric[2] - classical[2]) != 0


@pytest.mark.parametrize(
    "delta", (s.Rational(32, 125), s.Rational(1, 3), s.Rational(1, 2))
)
def test_independent_full_clock_composition_keeps_contact_once(delta):
    n, v, b = s.symbols("independent_n independent_v independent_b", real=True)
    zeta = v - s.log(1 + 2 * delta * ((1 + n) ** -2 - 1)) / 4
    metric = s.Matrix([n, zeta, b])
    M = s.Matrix([[2, 3, -1], [3, 5, 2], [-1, 2, 7]])
    current = s.Matrix([s.Rational(2, 7), s.Rational(-5, 11), s.Rational(3, 13)])
    functional = (current.T * metric)[0] + (metric.T * M * metric)[0] / 2
    observed = s.hessian(functional, (n, v, b)).subs({n: 0, v: 0, b: 0})
    first = s.Matrix([[1, 0, 0], [delta, 1, 0], [0, 0, 1]])
    contact = current[1] * (4 * delta * delta - 3 * delta)
    Rhat = M + s.diag(contact, 0, 0)
    assert (observed - first.T * Rhat * first).applyfunc(s.factor) == s.zeros(3)
    assert contact != 0
    assert observed != first.T * M * first
    assert observed != first.T * (Rhat + s.diag(contact, 0, 0)) * first


def test_independent_two_time_causal_block_solve_and_order_controls():
    # Two-time exact rational fixture, with nonsymmetric causal response and
    # distinct output density weights. This is not a discretization certificate
    # for the actual Proca continuum response.
    Jcan = s.BlockMatrix(
        [[s.zeros(2), s.eye(2)], [-s.eye(2), s.zeros(2)]]
    ).as_explicit()
    Cs = []
    Ds = []
    Fs = []
    for th, E, l, J, T, q, d in (
        (
            s.Rational(-2, 5),
            s.Rational(-1, 3),
            s.Rational(1, 10),
            2,
            s.Rational(1, 100),
            3,
            s.Rational(2, 5),
        ),
        (
            s.Rational(3, 5),
            s.Rational(-2, 7),
            s.Rational(1, 12),
            3,
            s.Rational(-1, 100),
            5,
            s.Rational(1, 3),
        ),
    ):
        w = -l * E
        n = s.Matrix(
            [
                [
                    -(2 * E * q + 3 * T) / (2 * J),
                    3 * th * l / (2 * J),
                    th / (2 * J),
                    -w / (2 * J),
                ]
            ]
        )
        C = n.col_join(s.Matrix([[1, 0, 0, 0]]) + d * n).col_join(
            s.Matrix([[0, 0, s.Rational(1, 2), 0]])
        )
        D = s.Matrix(
            [
                [-s.Rational(1, 2) / J, -d / (2 * J), 0],
                [-d / (2 * J), -d * d / (2 * J), 0],
                [0, 0, -s.Rational(3, 2)],
            ]
        )
        Cs.append(C)
        Ds.append(D)
        Fs.append(-Jcan * C.T)
    C = s.diag(*Cs)
    D = s.diag(*Ds)
    F = s.diag(*Fs)
    G = s.BlockMatrix(
        [
            [s.eye(4) / 7, s.zeros(4)],
            [s.Matrix(4, 4, lambda i, j: s.Rational(1 + i + 2 * j, 50)), s.eye(4) / 9],
        ]
    ).as_explicit()
    Q = s.BlockMatrix(
        [
            [s.Matrix([[2, 1, 0], [-1, 3, 1], [1, 0, 2]]) / 100, s.zeros(3)],
            [
                s.Matrix([[1, 3, -2], [4, -1, 1], [2, 1, 5]]) / 100,
                s.Matrix([[3, -1, 2], [1, 2, 0], [0, 2, 1]]) / 100,
            ],
        ]
    ).as_explicit()
    W = s.diag(*([s.Rational(2, 3)] * 3 + [s.Rational(3, 5)] * 3))
    Qbar = W * Q
    Zfree = s.Matrix([s.Rational(i + 1, 11) for i in range(8)])
    drive = s.Matrix([s.Rational(2 * i - 3, 31) for i in range(6)])
    Tmap = D + C * G * F
    S = s.eye(6) - Qbar * Tmap
    g = S.inv() * (drive + Qbar * C * Zfree)
    Z = Zfree + G * F * g
    full = s.BlockMatrix(
        [[s.eye(8), -G * F], [-Qbar * C, s.eye(6) - Qbar * D]]
    ).as_explicit()
    reference = full.inv() * Zfree.col_join(drive)
    assert all(s.factor(v) == 0 for v in reference - Z.col_join(g))
    assert all(s.factor(v) == 0 for v in g - drive - Qbar * (C * Z + D * g))
    without_D = (s.eye(6) - Qbar * C * G * F).inv() * (drive + Qbar * C * Zfree)
    wrong_order = (s.eye(6) - Q * W * Tmap).inv() * (drive + Q * W * C * Zfree)
    assert any(s.factor(v) != 0 for v in without_D - g)
    assert any(s.factor(v) != 0 for v in wrong_order - g)
    assert Q * W != W * Q and Q != Q.T
    g1 = Qbar * C * Zfree
    g2 = Qbar * Tmap * g1
    deleted_second = Qbar * C * G * F * g1
    assert any(s.factor(v) != 0 for v in g2 - deleted_second)


def test_derivative_losing_bound_alone_neither_gives_contraction_nor_excludes_inverse():
    lam = s.Symbol("lambda", positive=True)
    # A local spatial comparison multiplier (1-Delta)^5 satisfies the stated
    # H8 -> H-3 estimate, yet is unbounded on H8 -> H8. It is only a logical
    # norm counterexample, not the actual Proca response.
    assert s.simplify(lam ** s.Rational(-3, 2) * lam**5 / lam**4) == lam ** s.Rational(
        -1, 2
    )
    assert s.limit(lam**5, lam, s.oo) == s.oo
    # An unbounded Neumann factor can still have a bounded full inverse.
    inverse = 1 / (1 + s.Rational(3, 2) * lam**5)
    assert s.simplify((1 + s.Rational(3, 2) * lam**5) * inverse) == 1
    for value in (1, 2, 100):
        assert 0 < inverse.subs(lam, value) < 1
        assert (s.Rational(3, 2) * lam**5).subs(lam, value) > 1
    d = estimates.data()
    assert "not proof" in d["precise_estimate_gap"]
    assert "No numerical" in d["comparison_phase_input"]
