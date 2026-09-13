"""Independent full-profile, mixed-order, causal and current-force diagnostics."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_quantum_forced_constraints import forces
from p8_vacuum_affine_two_mass_coupled_inverse import assembly, audit, normal


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_every_complete_exact_identity(name):
    value = audit.residuals()[name]
    assert all(x == 0 for x in value) if isinstance(value, s.MatrixBase) else value == 0


@pytest.mark.parametrize("name", tuple(audit.gates()))
def test_every_complete_written_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_every_unsupported_scope_is_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_entire_current_parameters_frontier_and_counts():
    g = normal.reference
    assert audit.require_parameters(g.PROCA_MASS2, g.HEAVY_MASS2, g.KAPPA) == (
        g.PROCA_MASS2,
        g.HEAVY_MASS2,
        g.KAPPA,
    )
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 104
    assert audit.matching()[:-1] == audit.previous.matching()
    assert audit.validate_scope(audit.frontier(), audit.matching()) is True
    assert all(v is True for k, v in audit.controls().items() if k != "rejected_inputs")
    assert (
        len(audit.residuals()),
        audit.scalar_entry_count(),
        len(audit.gates()),
        audit.rejected_inputs(),
    ) == (67, 231, 38, 244)


@pytest.mark.parametrize("case", range(12))
def test_independent_unexpanded_whole_profile_density_and_duplicate_defect(case):
    with mp.workdps(65):
        dd = mp.mpf(1) / (2 * (1 + mp.mpf(case - 5) ** 2 / 100) ** 3)
        aa = (1 + mp.mpf(case - 5) ** 2 / 100) ** 2
        rho = mp.mpf(case + 2) / 7
        pressure = -mp.mpf(case + 1) / 11
        nD, nG = mp.mpf(case + 1) / 13, mp.mpf(2 - case) / 9
        vD, vG = mp.mpf(3 - case) / 8, mp.mpf(case + 2) / 17
        A, B = -pressure, -(rho + pressure) / 2

        def whole(e, f):
            N = 1 + e * nD + f * nG
            v = e * vD + f * vG
            return (
                aa**3
                * N
                * mp.exp(3 * v)
                * (1 + 2 * dd * (N**-2 - 1)) ** (-mp.mpf(3) / 4)
                * (A + B * (N**-2 - 1))
            )

        literal = mp.diff(whole, (0, 0), (1, 1))
        dJ = (21 * dd**2 - 3 * dd) * A / 2 + (1 - 6 * dd) * B
        T = (1 + 3 * dd) * A - 2 * B
        matched = aa**3 * (
            2 * dJ * nD * nG + 3 * T * (nD * vG + nG * vD) + 9 * A * vD * vG
        )
        assert abs(literal - matched) < mp.mpf("1e-58") * max(1, abs(literal))
        base, gaussianP, gaussianH = mp.mpf(3), mp.mpf(5), mp.mpf(7)
        current = base + matched
        fullH = gaussianH + literal
        once = current + gaussianP + fullH - literal
        assert abs(once - (base + gaussianP + fullH)) < mp.mpf("1e-58")
        assert abs((current + gaussianP + fullH) - once) > mp.mpf("1e-5")


@pytest.mark.parametrize("case", range(12))
def test_independent_full_forced_saddle_with_all_new_profile_coefficients(case):
    dd = s.Rational(1, 2) - s.Rational(case, 100)
    rho = s.Rational(case + 3, 71)
    pressure = -s.Rational(case + 1, 67)
    AH, BH = -pressure, -(rho + pressure) / 2
    dJH = (21 * dd**2 - 3 * dd) * AH / 2 + (1 - 6 * dd) * BH
    TH = (1 + 3 * dd) * AH - 2 * BH
    Jold, Aold, Told = s.Rational(5 + case, 3), -s.Rational(1, 59), s.Rational(1, 61)
    E = 1 - 3 * dd
    ell = s.Rational(1, 20)
    theta = s.Rational(case - 5, 13)
    q = s.Rational(case + 1, 7)
    N, V, Vd, S, Sd, B = (
        forces.n,
        forces.v,
        forces.vd,
        forces.sigma,
        forces.sd,
        forces.b,
    )
    wb = -ell * E
    oldfull = (
        -3 * Vd**2
        + (Jold + wb**2 / 2 - 3 * theta**2) * N**2
        + 6 * theta * N * Vd
        + Sd**2 / 2
        + wb * N * Sd
        - 3 * ell * Vd * S
        + 2 * B * (Vd - theta * N)
        + ell * B * S
        + q * V**2
        + 2 * E * q * N * V
        - q * S**2 / 2
        + 3 * Told * N * V
        + 9 * Aold * V**2 / 2
    )
    full = oldfull + dJH * N**2 + 3 * TH * N * V + 9 * AH * V**2 / 2
    equations = [
        s.diff(full, Vd) - forces.pv,
        s.diff(full, Sd) - forces.ps,
        s.diff(full, N) + forces.gn + dd * forces.gz,
        s.diff(full, B) + forces.gb,
    ]
    matrix, rhs = s.linear_eq_to_matrix(equations, (Vd, Sd, N, B))
    solution = matrix.inv() * rhs
    substitute = {
        forces.J: Jold + dJH,
        forces.A: Aold + AH,
        forces.T: Told + TH,
        forces.th: theta,
        forces.E: E,
        forces.l: ell,
        forces.w: wb,
        forces.delta: dd,
        forces.q: q,
        forces.H: s.Rational(case - 4, 19),
    }
    expected = s.Matrix(
        [
            forces.system()["solution"][v].subs(substitute, simultaneous=True)
            for v in (Vd, Sd, N, B)
        ]
    )
    assert (solution - expected).applyfunc(s.factor) == s.zeros(4, 1)
    assert matrix.det() == -8 * (Jold + dJH)
    wrong = expected.subs({forces.gn: 0, forces.gz: 0, forces.gb: 0})
    assert wrong != expected
    metric = s.Matrix([solution[2], V + dd * solution[2], solution[3]])
    C, D = metric.jacobian(forces.Z), metric.jacobian(forces.g)
    assert (C - forces.system()["C"].subs(substitute, simultaneous=True)).applyfunc(
        s.factor
    ) == s.zeros(3, 4)
    assert (D - forces.system()["D"].subs(substitute, simultaneous=True)).applyfunc(
        s.factor
    ) == s.zeros(3)
    assert D.rank() == 2
    direct_without_H = oldfull
    assert s.diff(full - direct_without_H, N, 2) == 2 * dJH


@pytest.mark.parametrize(
    "order,derivative", [(n, j) for n in (2, 4) for j in range(n + 1)]
)
def test_independent_full_local_output_primitive_and_upper_contact(order, derivative):
    t, u = normal.t, normal.u
    coefficient = 1 + 2 * u + 3 * u * u + u**4
    history = u ** (order + derivative + 2) * (1 + u)
    lhs = s.integrate(
        (t - u) ** (order - 1)
        * coefficient
        * s.diff(history, u, derivative)
        / s.factorial(order - 1),
        (u, 0, t),
    )
    kernel, local = normal.local_primitive_kernel(order, derivative, coefficient)
    rhs = s.integrate(kernel * history, (u, 0, t)) + local * history.subs(u, t)
    assert s.expand(lhs - rhs) == 0
    if derivative == order:
        assert s.expand(lhs - s.integrate(kernel * history, (u, 0, t))) != 0
    if derivative > 0:
        wrong = (
            coefficient.subs(u, t)
            * s.integrate(
                (t - u) ** (order - derivative - 1)
                * history
                / s.factorial(order - derivative - 1),
                (u, 0, t),
            )
            if derivative < order
            else coefficient.subs(u, t) * history.subs(u, t)
        )
        assert s.expand(lhs - wrong) != 0


@pytest.mark.parametrize("case", range(6))
def test_independent_unequal_time_clock_maps_and_output_density(case):
    t = s.Rational(case - 2, 9)
    u = -s.Rational(case + 1, 13)
    scale = lambda x: (1 + x * x) ** 2
    delta = lambda x: 1 / (2 * (1 + x * x) ** 3)
    A = lambda x: s.Matrix([[1, 0, 0], [delta(x), 1, 0], [0, 0, 1]])
    Qp = s.Matrix([[2, 1, 0], [-1, 3, 1], [2, 0, 1]])
    Qh = s.Matrix([[1, 2, 3], [0, 4, 1], [-2, 1, 5]]) / 7
    contact = s.Rational(case + 1, 17) * s.diag(1, 0, 0)
    profile = s.Matrix([[3, 2, 0], [2, -1, 0], [0, 0, 0]]) / 19
    fullclockH = A(t).T * (Qh + contact) * A(u) + profile
    recovered = A(t).T.inv() * (fullclockH - profile) * A(u).inv()
    assert recovered == Qh + contact
    output = (Qp + recovered) / scale(t) ** 3
    wrongdensity = (Qp + recovered) / scale(u) ** 3
    assert output != wrongdensity
    wrongleg = A(t).T.inv() * (fullclockH - profile)
    assert wrongleg != recovered
    assert (fullclockH - profile) / scale(t) ** 3 == A(t).T * (
        (Qh + contact) / scale(t) ** 3
    ) * A(u)


@pytest.mark.parametrize("case", range(4))
def test_independent_complete_mixed_order_noncommuting_inverse_both_products(case):
    # Finite matrix diagnostic only: no discretized continuum convergence is claimed.
    N = 2
    I = s.eye(N)
    C = s.Matrix([[1, 0], [s.Rational(1, case + 3), 1]])
    D4 = C**4
    D2 = C**2
    primitive = s.diag(D4, D4, D4, D2)
    density = s.diag(*[s.Rational(case + j + 2, case + j + 3) for j in range(8)])
    T0 = s.diag(s.Rational(-3, 2), -s.Rational(3, 2) / (1 + s.Rational(1, 16)) ** 6)
    P0 = -s.Matrix([[4, 0], [s.Rational(1, 3), 4]])
    H0 = -s.Matrix([[800, 0], [s.Rational(2, 5), 800]])
    P2 = -s.Matrix([[s.Rational(1, 30), 0], [s.Rational(1, 11), s.Rational(1, 30)]])
    H2 = -s.Matrix([[s.Rational(20, 3), 0], [s.Rational(1, 13), s.Rational(20, 3)]])
    L = s.kronecker_product(normal.L0, I)
    middle = L.T * s.diag(P0 + H0, s.Rational(8, 3) * (P2 + H2)) * L
    reference = s.diag(T0, middle, -I)
    remainder = s.zeros(8)
    for i in range(8):
        for j in range(8):
            if (i + j + case) % 3 == 0:
                remainder[i, j] = s.Rational((i + 1) * (j + 2), 10000 + case)
    assert reference * remainder != remainder * reference
    B0 = reference.inv()
    normal_inverse = (s.eye(8) + B0 * remainder).inv() * B0
    raw = density.inv() * primitive.inv() * (reference + remainder)
    inverse_raw = normal_inverse * primitive * density
    assert inverse_raw * raw == s.eye(8)
    assert raw * inverse_raw == s.eye(8)
    wrong_order = primitive * density * normal_inverse
    assert wrong_order * raw != s.eye(8)
    assert (P0 + H0).inv() != P0.inv() + H0.inv()
    assert (P2 + H2).inv() != P2.inv() + H2.inv()


@pytest.mark.parametrize("case", range(4))
def test_independent_full_current_force_schur_recovery_with_general_drive(case):
    values = {
        forces.J: s.Rational(case + 5, 3),
        forces.A: -s.Rational(case + 1, 71),
        forces.T: s.Rational(case + 2, 79),
        forces.th: s.Rational(case - 2, 7),
        forces.E: -s.Rational(1, 3),
        forces.l: s.Rational(1, 20),
        forces.w: s.Rational(1, 60),
        forces.delta: s.Rational(4, 9),
        forces.q: s.Rational(case + 1, 5),
        forces.H: s.Rational(case - 1, 11),
    }
    d = forces.system()
    C, F, D = [d[name].subs(values, simultaneous=True) for name in ("C", "F", "D")]
    G = s.Matrix([[1, 1, 0, 0], [0, 2, 1, 0], [1, 0, 2, 1], [0, 0, 1, 3]]) / 17
    Qp = s.Matrix([[1, 2, 0], [0, 1, 1], [1, 0, 2]]) / 101
    Qh = s.Matrix([[3, 0, 1], [1, 2, 0], [0, 1, 1]]) / 103
    Q = Qp + Qh
    h = s.Matrix([case + 1, 2 - case, 3, 4]) / 19
    e = s.Matrix([1, case + 2, 3 - case]) / 23
    full = s.BlockMatrix([[s.eye(4), -G * F], [-Q * C, s.eye(3) - Q * D]]).as_explicit()
    rhs = s.Matrix.vstack(G * h, e)
    solution = full.inv() * rhs
    Taux = D + C * G * F
    force = (s.eye(3) - Q * Taux).inv() * (e + Q * C * G * h)
    phase = G * h + G * F * force
    assert solution == s.Matrix.vstack(phase, force)
    metric = C * phase + D * force
    assert force == e + Q * metric
    assert (s.eye(3) - Q * Taux).inv() != (s.eye(3) - Taux * Q).inv()
    omitted = (s.eye(3) - Q * C * G * F).inv() * (e + Q * C * G * h)
    assert omitted != force


@pytest.mark.parametrize("value", ("1", "2", "5", "10", "100", "1e4", "1e8", "1e20"))
def test_independent_complete_weak_log_weight_inequality(value):
    with mp.workdps(65):
        lam = mp.mpf(value)
        # x=lambda*tau keeps the whole log singularity resolved even for huge lambda.
        integrand = lambda x: mp.exp(-x) * (1 + mp.log(lam) - mp.log(x)) / lam
        intervals = (
            [mp.mpf(0)]
            + [x for x in (mp.mpf(1), mp.mpf(5), mp.mpf(20), mp.mpf(100)) if x < lam]
            + [lam]
        )
        actual = mp.quad(integrand, intervals)
        bound = (2 + mp.log(lam)) / lam
        assert 0 < actual < bound
        assert bound <= 2 / mp.sqrt(lam)
        beta = mp.sqrt(lam) / 7
        chosen = (4 * beta + 1) ** 2
        assert 2 * beta / mp.sqrt(chosen) < mp.mpf("0.5")


@pytest.mark.parametrize("case", range(6))
def test_independent_actual_weighted_residual_and_nonzero_initial_kernel(case):
    t = assembly.t
    aa = (1 + t * t) ** 2
    H = s.diff(aa, t) / aa
    delta = 1 / (2 * (1 + t * t) ** 3)
    q = s.Rational(case + 1, 7) / aa**2
    ell = s.Rational(1, 20)
    fs = [(t + s.Rational(1, 2)) ** (j + 3) * (1 + (case + 1) * t) for j in range(4)]
    eta, w, c, r = fs
    En, Ev, Eb, Es = [1 + (j + 1) * t + (case + 1) * t * t for j in range(4)]
    W = lambda f: s.diff(f, t) + 3 * H * f
    n = s.diff(eta, t)
    v = w + H * eta - delta * n
    b = s.diff(c, t) + q * eta
    sigma = r + ell * eta
    adapted = [-W(En) + H * Ev + W(delta * Ev) + q * Eb + ell * Es, Ev, -W(Eb), Es]
    boundary = aa**3 * ((En - delta * Ev) * eta + Eb * c)
    raw = aa**3 * (En * n + Ev * v + Eb * b + Es * sigma)
    assert (
        s.cancel(
            raw - aa**3 * sum(x * y for x, y in zip(adapted, fs)) - s.diff(boundary, t)
        )
        == 0
    )
    assert W(aa**-3).simplify() == 0
    assert aa.subs(t, -s.Rational(1, 2)) ** -3 != 0
    assert s.diff((-6 * delta**2) * eta, t) != (-6 * delta**2) * s.diff(eta, t)


@pytest.mark.parametrize("case", range(6))
def test_independent_actual_geometric_leading_ratio_and_finite_transfer_phase(case):
    with mp.workdps(70):
        point = mp.mpf(case - 3) / 10
        P = mp.mpf(case + 1) / 3
        a = lambda t: (1 + t * t) ** 2
        derivative = 4 * point * (1 + point * point)
        lag = mp.mpf("1e-20")
        conformal = mp.quad(lambda v: 1 / a(v), [point, point + lag])
        ratio = lag**5 / (a(point + lag) ** 4 * a(point) * conformal**5)
        first = (ratio - 1) / lag
        expected = -3 * derivative / (2 * a(point))
        assert abs(first - expected) < mp.mpf("1e-17")
        phase = mp.exp(-1j * P * conformal)
        assert abs((phase - 1) / lag + 1j * P / a(point)) < mp.mpf("1e-17")
        assert abs((phase - 1) / lag) > mp.mpf("0.1")
