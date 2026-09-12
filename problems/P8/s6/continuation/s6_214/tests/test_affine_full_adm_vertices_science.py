"""Independent full ADM Legendre, matrix-exponential, Fourier and form-norm checks."""

import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_full_adm_vertices import adm, audit, bounds, chart, vertices
from p8_vacuum_affine_spatial_current import hamiltonian as old


def literal_spacetime_action(N, h, beta, A0, A, F0, Fij, m):
    v = s.sqrt(h.det())
    E = F0 - Fij.T * beta
    hi = h.inv()
    magnetic = sum(
        hi[i, k] * hi[j, l] * Fij[i, j] * Fij[k, l]
        for i in range(3)
        for j in range(3)
        for k in range(3)
        for l in range(3)
    )
    return (
        v * (E.T * hi * E)[0] / (2 * N)
        - N * v * magnetic / 4
        + m * m * v * (A0 - (beta.T * A)[0]) ** 2 / (2 * N)
        - N * m * m * v * (A.T * hi * A)[0] / 2
    )


@pytest.mark.parametrize("case", range(4))
def test_literal_nonzero_shift_Legendre_and_covariant_action(case):
    R = s.Matrix([[2, case, 1], [0, 3, 1], [1, -1, 2]])
    h = R.T * R + s.eye(3)
    N = s.Rational(5 + case, 4)
    beta = s.Matrix([1, -2, case + 1]) / 13
    A = s.Matrix([2, 1, -1]) / 7
    A0 = s.Symbol("A0", real=True)
    Fij = s.Matrix([[0, 1, -2], [-1, 0, 3], [2, -3, 0]]) / 11
    F0 = s.Matrix(s.symbols("f0:3", real=True))
    pi = s.Matrix(s.symbols("pi0:3", real=True))
    m = s.Integer(1000)
    v = s.sqrt(h.det())
    L = literal_spacetime_action(N, h, beta, A0, A, F0, Fij, m)
    expected_pi = v * h.inv() * (F0 - Fij.T * beta) / N
    actual_pi = s.Matrix([s.diff(L, f) for f in F0])
    assert all(s.simplify(x) == 0 for x in actual_pi - expected_pi)
    velocity = N * h * pi / v + Fij.T * beta
    H = (pi.T * F0)[0] - L
    H = s.simplify(H.subs(dict(zip(F0, velocity)), simultaneous=True))
    divpi = s.Symbol("divpi", real=True)
    integrated = H - A0 * divpi
    solution = adm.temporal_constraint(N, v, beta, A, divpi, m)
    assert s.simplify(s.diff(integrated, A0).subs(A0, solution)) == 0
    B = s.Matrix([Fij[1, 2], Fij[2, 0], Fij[0, 1]])
    expected = N * (
        (pi.T * h * pi)[0] / (2 * v)
        + adm.magnetic_energy(h, B)
        + m * m * v * (A.T * h.inv() * A)[0] / 2
        + divpi**2 / (2 * m * m * v)
    )
    expected += (beta.T * Fij * pi)[0] - (beta.T * A)[0] * divpi
    assert s.simplify(integrated.subs(A0, solution) - expected) == 0
    # Independently contract the original four-metric field strength and vector.
    g = s.zeros(4)
    g[0, 0] = N * N - (beta.T * h * beta)[0]
    g[0, 1:] = -beta.T * h
    g[1:, 0] = -h * beta
    g[1:, 1:] = -h
    gi = g.inv()
    F = s.zeros(4)
    F[0, 1:] = F0.T
    F[1:, 0] = -F0
    F[1:, 1:] = Fij
    avec = s.Matrix([A0, *A])
    cov = (
        N
        * v
        * (
            -sum(
                gi[i, k] * gi[j, l] * F[i, j] * F[k, l]
                for i in range(4)
                for j in range(4)
                for k in range(4)
                for l in range(4)
            )
            / 4
            + m * m * (avec.T * gi * avec)[0] / 2
        )
    )
    assert s.simplify(cov - L) == 0
    assert s.simplify(solution + N * divpi / (m * m * v)) == (beta.T * A)[0]


@pytest.mark.parametrize("case", range(4))
def test_general_SPD_magnetic_density_identity(case):
    R = s.Matrix([[1, case + 1, -1], [2, 1, 0], [1, -1, 3]])
    h = R.T * R + s.eye(3)
    B = s.Matrix(s.symbols("B0:3", real=True))
    # F_ij=epsilon_ijk B_k, independently built entry by entry.
    F = s.Matrix([[0, B[2], -B[1]], [-B[2], 0, B[0]], [B[1], -B[0], 0]])
    hi = h.inv()
    v = s.sqrt(h.det())
    direct = (
        v
        * sum(
            hi[i, k] * hi[j, l] * F[i, j] * F[k, l]
            for i in range(3)
            for j in range(3)
            for k in range(3)
            for l in range(3)
        )
        / 4
    )
    assert s.simplify(direct - adm.magnetic_energy(h, B)) == 0


@pytest.mark.parametrize("case", range(3))
def test_literal_full_symmetric_log_metric_second_vertices(case):
    k = s.Matrix(([0, 0, 0], [3, 4, 12], [1000, -2000, 500])[case])
    q = s.Matrix(([1, 0, -1], [7, -2, 3], [-500, 2000, 3000])[case])
    D = s.Matrix([[1, 2, -1], [2, 3, 1], [-1, 1, 2]]) / 7
    G = s.Matrix([[2, 1, 3], [1, 1, -2], [3, -2, -1]]) / 11
    assert D * G != G * D and s.trace(D) != 0 and s.trace(G) != 0
    e, z = s.symbols("e z", real=True)
    a = s.Rational(1) + s.Rational(case, 5)
    m = s.Integer(1000)
    nD = s.Rational(2, 5)
    nG = s.Rational(-3, 7)
    Q = e * D + z * G
    E = s.eye(3) + Q + Q * Q / 2
    Ei = s.eye(3) - Q + Q * Q / 2
    halftrace = s.trace(Q) / 2
    vol = 1 + halftrace + halftrace**2 / 2
    voli = 1 - halftrace + halftrace**2 / 2
    # Literal h/sqrt(h),sqrt(h)h^-1 and1/sqrt(h), not B_D feature definitions.
    kinetic = E * voli / a
    mass = a * m * m * Ei * vol
    M = (1 + e * nD + z * nG) * s.diag(
        mass + old.cross(k).T * kinetic * old.cross(q),
        kinetic + voli * k * q.T / (a**3 * m * m),
    )
    first = M.applyfunc(lambda v: s.expand(v).coeff(e, 1).subs(z, 0))
    second = M.applyfunc(lambda v: s.expand(v).coeff(e, 1).coeff(z, 1))
    assert first == vertices.vertex(k, q, nD, s.zeros(3, 1), D, a, m)
    assert second == vertices.contact(k, q, nD, D, nG, G, a, m)


@pytest.mark.parametrize("case", range(4))
def test_literal_Fourier_shift_flux_and_longitudinal_omission(case):
    k = s.Matrix(([0, 0, 0], [3, 4, 12], [1000, -2000, 500], [10**9, 3, 4])[case])
    q = s.Matrix(([1, 0, -1], [7, -2, 3], [-500, 2000, 3000], [3, 10**9, -4])[case])
    b = s.Matrix([1, -2, 3]) / 13
    M = vertices.vertex(k, q, 0, b, s.zeros(3), s.Rational(5, 4), 1000)
    A = s.Matrix(s.symbols("A0:3", real=True))
    pi = s.Matrix(s.symbols("P0:3", real=True))
    # Created/destroyed Fourier pairing: derivative on A(q) and Pi(-k).
    F = s.I * (q * A.T - A * q.T)
    direct = (b.T * F * pi)[0] + s.I * (b.T * A)[0] * (k.T * pi)[0]
    assert s.expand(direct - (pi.T * M[3:, :3] * A)[0]) == 0
    assert (
        M
        == vertices.vertex(q, k, 0, b, s.zeros(3), s.Rational(5, 4), 1000).conjugate().T
    )
    if k != s.zeros(3, 1):
        incomplete = (b.T * F * pi)[0]
        assert s.expand(direct - incomplete) != 0


@pytest.mark.parametrize("time", (-0.5, 0, 0.25, 0.5))
@pytest.mark.parametrize("case", range(5))
@pytest.mark.parametrize("complex_direction", (False, True))
def test_all_momentum_energy_relative_first_and_second_bounds(
    time, case, complex_direction
):
    a = (1 + time * time) ** 2
    m = 1000.0
    k = np.array(
        ([0, 0, 0], [3, 4, 12], [1000, -2000, 500], [1e9, 3, 4], [1e12, -2e12, 3e12])[
            case
        ],
        dtype=float,
    )
    q = np.array(
        ([1, 0, -1], [7, -2, 3], [-500, 2000, 3000], [3, 1e9, -4], [-1e12, 3e12, 2e12])[
            case
        ],
        dtype=float,
    )
    D = np.array([[1, 2, -1], [2, 3, 1], [-1, 1, 2]], dtype=float) / 7
    G = np.array([[2, 1, 3], [1, 1, -2], [3, -2, -1]], dtype=float) / 11
    b = np.array([1, -2, 3], dtype=float) / 13
    bG = np.array([-2, 1, 4], dtype=float) / 17
    n = 0.2
    nG = -0.3
    if complex_direction:
        D = D + 1j * G / 3
        G = G + 1j * D / 5
        b = b + 1j * bG / 7
        bG = bG - 1j * b / 11
        n = n + 0.1j
        nG = nG - 0.2j

    def polar(momentum):
        r = np.linalg.norm(momentum)
        P = np.outer(momentum / r, momentum / r) if r else np.diag([1.0, 0, 0])
        T = np.eye(3) - P
        V = P / (np.sqrt(a) * m) + T / np.sqrt(a * m * m + r * r / a)
        K = np.sqrt(a) * T + P / np.sqrt(1 / a + r * r / (a**3 * m * m))
        normal = np.zeros((6, 6))
        normal[:3, :3] = V
        normal[3:, 3:] = K
        F = np.array(
            old.features(s.Matrix(momentum), s.Float(a), s.Float(m)), dtype=float
        )
        # Stable polar feature construction avoids subtracting large6x6 blocks.
        return F @ normal

    Fk, Fq = polar(k), polar(q)
    assert np.linalg.norm(Fk.T @ Fk - np.eye(6), 2) < 1e-5
    feature = np.array(
        vertices.feature(s.sympify(n), s.Matrix(b), s.Matrix(D), s.Float(a)),
        dtype=complex,
    )
    contact = np.array(
        vertices.contact_feature(s.sympify(n), s.Matrix(D), s.sympify(nG), s.Matrix(G)),
        dtype=complex,
    )
    sigma = abs(n) + 2.5 * np.linalg.norm(D, 2) + a * np.linalg.norm(b)
    sigmaG = abs(nG) + 2.5 * np.linalg.norm(G, 2) + a * np.linalg.norm(bG)
    assert np.linalg.norm(Fk.T @ feature @ Fq, 2) <= sigma * (1 + 1e-5)
    assert np.linalg.norm(Fk.T @ contact @ Fq, 2) <= sigma * sigmaG * (1 + 1e-5)
    shift = np.array(vertices.shift_feature(s.Matrix(b), s.Float(a)), dtype=complex)
    assert abs(np.linalg.norm(shift, 2) - a * np.linalg.norm(b)) < 1e-12


def test_nonzero_trace_constraint_first_and_second_omission():
    D = 2 * s.eye(3)
    G = 4 * s.eye(3)
    first = vertices.spatial_feature(D)
    contact = vertices.spatial_contact_feature(D, G)
    assert first[9, 9] == -3
    assert contact[9, 9] == 18
    assert vertices.spatial_feature(s.diag(1, -1, 0))[9, 9] == 0


def test_zero_ADM_contacts_are_not_zero_four_metric_chart_contacts():
    a = s.Symbol("a", positive=True)
    b = s.Matrix([1, -2, 3])
    zero = s.zeros(3)
    assert vertices.contact_feature(1, zero, 1, zero) == s.zeros(10)
    assert chart.second(1, s.zeros(3, 1), zero, 1, s.zeros(3, 1), zero, a)[0, 0] == 2
    assert chart.second(0, b, zero, 0, b, zero, a)[0, 0] == -28 * a * a
    assert vertices.shift_feature(b, a) != s.zeros(10)


@pytest.mark.parametrize("module", (adm, vertices, chart, bounds))
def test_complete_full_ADM_core(module):
    d = module.data()
    for value in d["checks"].values():
        values = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.cancel(x) == 0 for x in values)
    assert all(d["gates"].values())


def test_complex_reverse_pair_reality_for_all_metric_vertices():
    k = s.Matrix([3, 4, 12])
    q = s.Matrix([-2, 5, 7])
    Q = (
        s.Matrix(
            [
                [1 + s.I, 2 - s.I, 3],
                [2 - s.I, 3, 1 + 2 * s.I],
                [3, 1 + 2 * s.I, 2 - s.I],
            ]
        )
        / 13
    )
    beta = s.Matrix([1 + s.I, -2 + s.I / 3, 3 - 2 * s.I]) / 11
    n = s.Rational(2, 7) + s.I / 5
    first = vertices.vertex(k, q, n, beta, Q)
    reverse = vertices.vertex(q, k, s.conjugate(n), beta.conjugate(), Q.conjugate())
    assert (first - reverse.conjugate().T).applyfunc(s.cancel) == s.zeros(6)
    assert first != first.conjugate().T
    second = vertices.contact(k, q, n, Q, 2 * n, 3 * Q)
    reverse_second = vertices.contact(
        q, k, s.conjugate(n), Q.conjugate(), 2 * s.conjugate(n), 3 * Q.conjugate()
    )
    assert (second - reverse_second.conjugate().T).applyfunc(s.cancel) == s.zeros(6)


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_all_full_ADM_exact_residuals(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(v == 0 for v in entries)


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_all_unsupported_full_ADM_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_exact_full_ADM_audit_counts_and_scope():
    assert len(audit.residuals()) == 39 and audit.scalar_entry_count() == 697
    assert len(audit.gates()) == 35 and all(audit.gates().values())
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 136
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 70
    assert audit.gates()["growing_finite_band_bound_not_UV_limit"]


def test_corrected_matching_history_and_primitive_statuses_unchanged():
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len({row["id"] for row in audit.matching()}) == len(audit.matching())
    assert audit.previous.metadata_errata()[0]["checkpoint"] == "S6.211"


def test_full_ADM_scope_mutation_does_not_mutate_previous_cache():
    old = audit.previous.matching()
    rows = audit.matching()
    rows[-1]["status"] = "COMPLETE"
    with pytest.raises(ValueError):
        audit.validate_scope(audit.frontier(), rows)
    assert audit.previous.matching() == old
