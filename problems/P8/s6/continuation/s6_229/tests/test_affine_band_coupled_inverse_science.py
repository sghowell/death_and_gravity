"""Independent coupled principal, primitive, finite-transfer and causal-order checks."""

import pytest
import sympy as s
from p8_vacuum_affine_band_coupled_inverse import audit, coordinates, inverse, matching
from p8_vacuum_affine_quantum_forced_constraints import forces


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


def test_original_frontier_and_actual_restricted_scope():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 85
    assert audit.frontier() == audit.previous.frontier()
    assert all(row["status"] != "COMPLETE" for row in audit.frontier())
    assert "ACTUAL_SMOOTH_BOUNDED_EXTERNAL_MOMENTUM" in audit.ITEM["status"]
    assert "NOT_UNRESTRICTED_GRAPH" in audit.ITEM["status"]
    assert inverse.require_ball(s.Rational(7, 3)) == s.Rational(7, 3)
    assert inverse.require_ball(10**12) == 10**12


def test_independent_classical_top_terms_and_required_matter_shift():
    t = s.Symbol("t", real=True)
    eta, w, c, r = [s.Function(name)(t) for name in ("eta", "w", "c", "r")]
    delta, ell, H = [s.Function(name)(t) for name in ("delta", "ell", "H")]
    W = lambda f: s.diff(f, t) + 3 * H * f
    # Independent highest-derivative portion of the complete sourced classical action.
    vdot = (
        s.diff(w, t)
        - delta * s.diff(eta, t, 2)
        - s.diff(delta, t) * s.diff(eta, t)
        + H * s.diff(eta, t)
    )
    sigmadot = s.diff(r + ell * eta, t)
    wb = -ell * (1 - 3 * delta)
    lag = (
        -3 * vdot**2
        + sigmadot**2 / 2
        + wb * s.diff(eta, t) * sigmadot
        - 3 * ell * vdot * (r + ell * eta)
        + 2 * s.diff(c, t) * vdot
    )
    eq = lambda f: s.expand(
        s.diff(lag, f)
        - W(s.diff(lag, s.diff(f, t)))
        + W(W(s.diff(lag, s.diff(f, t, 2))))
    )
    principal = s.Matrix(
        [
            [s.factor(eq(f).coeff(s.diff(g, t, n))) for g in (eta, w, c, r)]
            for f, n in zip((eta, w, c, r), (4, 4, 4, 2))
        ]
    )
    assert principal == s.diag(-6 * delta**2, 0, 0, -1)
    no_shift = (
        s.diff(r, t) ** 2 / 2 + wb * s.diff(eta, t) * s.diff(r, t) - 3 * ell * vdot * r
    )
    old_r = s.expand(s.diff(no_shift, r) - W(s.diff(no_shift, s.diff(r, t))))
    assert s.factor(old_r.coeff(s.diff(eta, t, 2)) - ell) == 0


@pytest.mark.parametrize("n,j", [(n, j) for n in (2, 4) for j in range(n + 1)])
def test_independent_variable_coefficient_row_primitive_including_principal(n, j):
    x, y = s.symbols("x y", nonnegative=True)
    source = y**7 * (1 + y + y * y)
    coef = 1 + 2 * y + 3 * y * y + y**4
    original = s.integrate(
        (x - y) ** (n - 1) * coef * s.diff(source, y, j) / s.factorial(n - 1), (y, 0, x)
    )
    start = 1 if j == n else 0
    kernel = sum(
        (-1) ** r
        * s.binomial(j, r)
        * (x - y) ** (n - 1 - j + r)
        * s.diff(coef, y, r)
        / s.factorial(n - 1 - j + r)
        for r in range(start, j + 1)
    )
    rebuilt = s.integrate(kernel * source, (y, 0, x))
    if j == n:
        rebuilt += coef.subs(y, x) * source.subs(y, x)
    assert s.expand(original - rebuilt) == 0
    if j:
        incorrect = coef.subs(y, x) * s.integrate(
            (x - y) ** (n - 1) * s.diff(source, y, j) / s.factorial(n - 1), (y, 0, x)
        )
        assert s.expand(original - incorrect) != 0


def test_independent_dimensional_invariant_contractions_and_finite_transfer_phase():
    d = s.Symbol("d", positive=True)
    C = 1 - 2 / d + 1 / (2 * (d + 2))
    trace = (d - 1) * (d * d - 5 * d + 8) / 4
    B = s.factor((trace - d * C) / d**2)
    phys = C * s.Matrix([[12, -4], [-4, 4]]) + B * s.Matrix([[36, -12], [-12, 4]])
    assert phys.subs(d, 3) == matching.data()["fixed_physical_scalar_matrix"]
    assert (
        phys.diff(d).subs(d, 3).applyfunc(s.cancel)
        == matching.data()["fixed_physical_first_dimension_jet"]
    )
    assert s.diff(B, d).subs(d, 3) == s.Rational(4, 225)
    x, P, a, m, tau = s.symbols("inverse_k P a m tau", positive=True)
    mu = s.Symbol("direction_cosine", real=True)
    # Sum of massive conformal frequencies at k and -k+P, not a massless state.
    scaled = s.sqrt(1 + a * a * m * m * x * x) + s.sqrt(
        1 - 2 * x * P * mu + x * x * (P * P + a * a * m * m)
    )
    expansion = s.series(scaled, x, 0, 3).removeO()
    assert (
        s.expand(
            expansion
            - (2 - P * mu * x + (a * a * m * m + P * P * (1 - mu * mu) / 2) * x * x)
        )
        == 0
    )
    assert s.diff(s.exp(-s.I * P * mu * tau), tau).subs(tau, 0) == -s.I * P * mu
    assert s.exp(-s.I * P * mu * tau).subs(tau, 0) == 1


def test_independent_weighted_adjoint_and_prepared_constraint_control():
    t = s.Symbol("t", real=True)
    eta, w, c, r = [s.Function(name)(t) for name in ("eta", "w", "c", "r")]
    H, delta, q, ell = [s.Function(name)(t) for name in ("H", "delta", "q", "ell")]
    En, Ev, Eb, Es = [s.Function(name)(t) for name in ("En", "Ev", "Eb", "Es")]
    W = lambda f: s.diff(f, t) + 3 * H * f
    original = (
        En * s.diff(eta, t)
        + Ev * (w + H * eta - delta * s.diff(eta, t))
        + Eb * (s.diff(c, t) + q * eta)
        + Es * (r + ell * eta)
    )
    boundary = (En - delta * Ev) * eta + Eb * c
    adj = (
        (-W(En) + H * Ev + W(delta * Ev) + q * Eb + ell * Es) * eta
        + Ev * w
        - W(Eb) * c
        + Es * r
    )
    assert s.expand(original - adj - W(boundary)) == 0
    a = (1 + t * t) ** 2
    bad = 1 / a**3
    assert s.factor(s.diff(bad, t) + 3 * s.diff(a, t) / a * bad) == 0
    assert bad.subs(t, -s.Rational(1, 2)) != 0


def test_independent_full_mixed_row_inverse_and_right_density_order():
    D = s.Matrix([[2, 0], [-1, 3]])
    I2 = s.eye(2)
    L = s.Matrix([[1, -s.Rational(1, 3)], [0, -1]])
    Ft = s.Matrix([[-4, 0], [s.Rational(1, 5), -4]])
    Fs = s.Matrix([[-s.Rational(1, 30), 0], [s.Rational(1, 7), -s.Rational(1, 30)]])
    gamma = s.Rational(1, 11)
    Lbig = s.kronecker_product(L, I2)
    Qref = gamma * Lbig.T * s.diag(Ft, s.Rational(8, 3) * Fs) * Lbig
    A = s.diag(-2, -3)
    Fref = s.diag(A, Qref, -I2)
    B0 = s.diag(
        A.inv(),
        Lbig.inv()
        * s.diag(Ft.inv(), s.Rational(3, 8) * Fs.inv())
        * Lbig.T.inv()
        / gamma,
        -I2,
    )
    assert Fref * B0 == s.eye(8) and B0 * Fref == s.eye(8)
    rowI = s.diag(D**-4, D**-4, D**-4, D**-2)
    V = s.zeros(8)
    for i in range(4):
        for j in range(4):
            V[2 * i : 2 * i + 2, 2 * j : 2 * j + 2] = s.Matrix(
                [
                    [(i + 1) / (s.Integer(19) * (j + 2)), 0],
                    [(i + j + 1) / s.Integer(37), (j + 1) / (s.Integer(23) * (i + 2))],
                ]
            )
    assert B0 * V != V * B0
    E = (s.eye(8) + B0 * V).inv() * B0 * rowI
    T = rowI.inv() * (Fref + V)
    assert E * T == s.eye(8) and T * E == s.eye(8)
    assert (s.eye(8) + V * B0).inv() * B0 * rowI != E
    assert rowI * (s.eye(8) + B0 * V).inv() * B0 != E
    W = s.kronecker_product(s.eye(4), s.diag(s.Rational(2, 5), s.Rational(3, 7)))
    raw = W.inv() * T
    assert raw * (E * W) == s.eye(8) and (E * W) * raw == s.eye(8)
    assert W * E != E * W


def test_independent_variable_pivot_commutator_and_actual_auxiliary_feedback():
    t = s.Symbol("time", real=True)
    delta = 1 / (2 * (1 + t * t) ** 3)
    pivot = -6 * delta**2
    eta = t**4
    assert (
        s.factor(
            s.diff(pivot * eta, t) - pivot * s.diff(eta, t) - s.diff(pivot, t) * eta
        )
        == 0
    )
    assert s.diff(pivot, t) * eta != 0
    data = forces.system()
    sets = [
        {
            forces.th: s.Rational(1 + i, 5),
            forces.E: s.Rational(2 - i, 7),
            forces.l: s.Rational(1, 10 + i),
            forces.w: -s.Rational(2 + i, 17),
            forces.J: s.Integer(2 + i),
            forces.A: s.Rational(1, 29 + i),
            forces.T: s.Rational(1, 31 + i),
            forces.delta: s.Rational(1, 2 + i),
            forces.H: s.Rational(i - 1, 13),
            forces.q: s.Integer(3 + 2 * i),
        }
        for i in range(2)
    ]

    def lift(name):
        matrix = data[name]
        out = s.zeros(2 * matrix.rows, 2 * matrix.cols)
        for i in range(matrix.rows):
            for j in range(matrix.cols):
                out[2 * i : 2 * i + 2, 2 * j : 2 * j + 2] = s.diag(
                    *(matrix[i, j].subs(values) for values in sets)
                )
        return out

    K, C, F, Da = [lift(name) for name in ("K", "C", "F", "D")]
    G = (s.kronecker_product(s.eye(4), s.Matrix([[7, 0], [-1, 8]])) - K).inv()
    assert Da.rank() == 4
    Q = s.zeros(6)
    for i in range(3):
        for j in range(3):
            Q[2 * i : 2 * i + 2, 2 * j : 2 * j + 2] = s.Matrix(
                [
                    [(i + 1) / (s.Integer(41) * (j + 2)), 0],
                    [(j + 2) / s.Integer(53), (j + 1) / (s.Integer(47) * (i + 2))],
                ]
            )
    Ta = Da + C * G * F
    B = s.eye(6) - Q * Ta
    h = s.Matrix([s.Rational(i + 1, 19) for i in range(8)])
    e = s.Matrix([s.Rational(7 - i, 23) for i in range(6)])
    Zh = G * h
    g = B.inv() * (e + Q * C * Zh)
    Z = Zh + G * F * g
    metric = C * Z + Da * g
    assert G.inv() * Z - F * g == h
    assert g == e + Q * metric
    block = G.inv().row_join(-F).col_join((-Q * C).row_join(s.eye(6) - Q * Da))
    assert block.inv() * h.col_join(e) == Z.col_join(g)
    assert (s.eye(6) - Q * (C * G * F)).inv() * (e + Q * C * Zh) != g
    assert (s.eye(6) - Ta * Q).inv() * (e + Q * C * Zh) != g


@pytest.mark.parametrize(
    "call",
    [
        coordinates.principal_data,
        coordinates.ward_data,
        coordinates.residual_data,
        matching.data,
        inverse.data,
    ],
)
def test_each_source_packet_uses_exact_zero_and_true_gates(call):
    data = call()
    for value in data["checks"].values():
        vv = (
            value.applyfunc(s.cancel)
            if isinstance(value, s.MatrixBase)
            else s.cancel(value)
        )
        assert vv == (s.zeros(*vv.shape) if isinstance(vv, s.MatrixBase) else 0)
    assert all(bool(value) for value in data["gates"].values())
