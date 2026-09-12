"""Independent canonical branch, full ADM, geometry, curvature and WKB checks."""

import importlib.util
import sys
from functools import cache
from pathlib import Path

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_full_adm_vertices import vertices
from p8_vacuum_affine_ordered_scalar_symbol import (
    density,
    geometry,
    jets,
    matching,
    pairs,
)


def helper(name):
    key = "p8_s216_test_helper_" + name
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(
            key, Path(__file__).parent / (name + ".py")
        )
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


@pytest.mark.parametrize("di", range(3))
@pytest.mark.parametrize("gi", range(3))
def test_literal_Fock_and_real_covariance_fix_complex_pair_orientation(di, gi):
    ann = s.zeros(3)
    ann[0, 1] = 1
    ann[1, 2] = s.sqrt(2)
    q = (ann + ann.T) / s.sqrt(2)
    p = -s.I * (ann - ann.T) / s.sqrt(2)
    U = s.Matrix(
        [[s.Rational(3, 5), s.Rational(4, 5)], [-s.Rational(4, 5), s.Rational(3, 5)]]
    )
    J = s.Matrix([[0, 1], [-1, 0]])
    basis = (s.diag(1, 0), s.diag(0, 1), s.Matrix([[0, 1], [1, 0]]))
    D, G = basis[di], basis[gi]

    def observable(M, z):
        return sum(
            (M[i, j] * z[i] * z[j] / 2 for i in range(2) for j in range(2)), s.zeros(3)
        )

    zd = [U[0, 0] * q + U[0, 1] * p, U[1, 0] * q + U[1, 1] * p]
    HD, HG = observable(D, zd), observable(G, [q, p])
    literal = s.expand(s.I * (HD * HG - HG * HD)[0, 0])
    cov = -s.trace(D * U * (J * G + G * J.T) * U.T) / 4
    u = s.Matrix([1, -s.I]) / s.sqrt(2)
    bd = (u.T * U.T * D * U * u)[0]
    bg = (u.T * G * u)[0]
    assert s.simplify(literal - cov) == 0
    assert s.simplify(literal - s.im(s.conjugate(bd) * bg)) == 0
    phase = (s.Rational(3, 5) + s.I * s.Rational(4, 5)) ** 2
    stripped = (u.T * D * u)[0]
    frozen = -s.im(s.conjugate(stripped) * bg * s.conjugate(phase))
    if (di, gi) == (0, 2):
        assert s.simplify(literal - frozen) != 0


@pytest.mark.parametrize("n", range(1, 7))
def test_positive_phase_integral_and_every_endpoint(n):
    # Direct polynomial integration uses a preparation function with n zero
    # left derivatives; its nth derivative need not vanish at readout.
    t, x = s.symbols("t x", real=True)
    b = x**n
    Omega = s.Integer(3)
    exact = s.exp(s.I * Omega * t) * s.integrate(s.exp(-s.I * Omega * x) * b, (x, 0, t))
    B = sum(s.I * (-s.I) ** j * s.diff(t**n, t, j) / Omega ** (j + 1) for j in range(n))
    rem = (
        (-s.I) ** n
        * s.exp(s.I * Omega * t)
        * s.integrate(s.exp(-s.I * Omega * x) * s.diff(b, x, n) / Omega**n, (x, 0, t))
    )
    assert s.simplify(s.expand_complex(exact - B - rem)) == 0


@cache
def canonical_pairs():
    a = s.Rational(3, 2)
    m = s.Integer(2)
    r = s.Integer(5)
    k = s.Matrix([3, 4, 0])
    ell = s.Matrix([0, 5, 12])
    left = (s.Matrix([4, -3, 0]) / 5, s.Matrix([0, 0, 1]), k / 5)
    right = (s.Matrix([1, 0, 0]), s.Matrix([0, 12, -5]) / 13, ell / 13)
    Wl = (s.Rational(3, 2), s.Rational(3, 2), s.Rational(7, 5))
    Wr = (s.Integer(2), s.Integer(2), s.Rational(9, 4))
    Pl = tuple(-s.I * w - s.Rational(i + 1, 7) for i, w in enumerate(Wl))
    Pr = tuple(-s.I * w - s.Rational(i + 2, 11) for i, w in enumerate(Wr))
    O1 = s.sqrt((k.T * k)[0] + a * a * m * m) / r
    O2 = s.sqrt((ell.T * ell)[0] + a * a * m * m) / r
    Ds = (s.eye(3), s.diag(-1, -1, 2), s.Matrix([[1, 2, -1], [2, 3, 4], [-1, 4, 5]]))
    out = {}
    for h, D in enumerate(Ds):
        M = vertices.vertex(k, -ell, 0, s.zeros(3, 1), D, a, m)
        for i, ek in enumerate(left):
            Kl = 1 / a if i < 2 else (1 + (k.T * k)[0] / (a * a * m * m)) / a
            fl = s.sqrt(Kl / (2 * r * Wl[i] / a))
            pil = (r / a) * Pl[i] * fl / Kl
            vl = s.Matrix([*(fl * ek), *(pil * ek)])
            for j, el in enumerate(right):
                Kr = 1 / a if j < 2 else (1 + (ell.T * ell)[0] / (a * a * m * m)) / a
                fr = s.sqrt(Kr / (2 * r * Wr[j] / a))
                pir = (r / a) * Pr[j] * fr / Kr
                vr = s.Matrix([*(fr * el), *(pir * el)])
                sector = ("T" if i < 2 else "L") + ("T" if j < 2 else "L")
                actual = -(vl.T * M * vr)[0]
                predicted = (
                    r
                    / a
                    * pairs.readout(
                        k / r,
                        ell / r,
                        ek,
                        el,
                        D,
                        O1,
                        O2,
                        Wl[i],
                        Wr[j],
                        Pl[i],
                        Pr[j],
                        m * a / r,
                        sector,
                    )
                )
                out[h, i, j] = s.simplify(actual - predicted)
    return out


@pytest.mark.parametrize("h", range(3))
@pytest.mark.parametrize("i", range(3))
@pytest.mark.parametrize("j", range(3))
def test_full_nonzero_trace_original_ADM_pair(h, i, j):
    assert canonical_pairs()[h, i, j] == 0


def scalar_directions(dim, ch):
    I = s.eye(dim)
    S = s.diag(*([-1] * (dim - 1) + [dim - 1]))
    return {"trace_trace": (I, I), "trace_scalar": (I, S), "scalar_trace": (S, I)}[ch]


@cache
def literal_geometry(dim, ch):
    u, y = s.symbols("u y", real=True)
    k = s.Matrix([s.sqrt(1 - u * u)] + [0] * (dim - 2) + [u])
    ell = -k + s.Matrix([0] * (dim - 1) + [y])
    ell2 = 1 - 2 * u * y + y * y
    P = s.eye(dim) - k * k.T
    Q = s.eye(dim) - ell * ell.T / ell2
    I = s.eye(dim)
    D, G = scalar_directions(dim, ch)
    td, tg = s.trace(D), s.trace(G)
    B, C = D - td * I / 2, G - tg * I / 2

    def field_strength_matrix(M):
        out = s.zeros(dim)
        for i in range(dim):
            F = k * I[:, i].T - I[:, i] * k.T
            for j in range(dim):
                L = ell * I[:, j].T - I[:, j] * ell.T
                out[i, j] = s.trace(M) * s.trace(F.T * L) / 4 - s.trace(F.T * M * L)
        return out

    MD, MG = field_strength_matrix(D), field_strength_matrix(G)
    actual = {
        "00": s.trace(P * B * Q * C),
        "01": s.trace(P * B * Q * MG.T),
        "10": s.trace(P * MD * Q * C),
        "11": s.trace(P * MD * Q * MG.T),
        "TL": (ell.T * B * P * C * ell)[0] / ell2,
        "LT": (k.T * B * Q * C * k)[0],
        "LL": (k.T * B * ell)[0] * (k.T * C * ell)[0] / ell2,
        "LC": tg * (k.T * B * ell)[0],
        "CL": td * (k.T * C * ell)[0],
        "CC": td * tg * ell2,
    }
    expected = geometry.contractions(ch)
    return {
        key: s.simplify(
            value - expected[key].subs({geometry.d: dim, geometry.u: u, geometry.y: y})
        )
        for key, value in actual.items()
    }


@pytest.mark.parametrize("dim", (3, 4, 5, 6))
@pytest.mark.parametrize("ch", density.CHANNELS)
@pytest.mark.parametrize(
    "key", ("00", "01", "10", "11", "TL", "LT", "LL", "LC", "CL", "CC")
)
def test_integer_dimension_literal_full_field_strength_geometry(dim, ch, key):
    assert literal_geometry(dim, ch)[key] == 0


@cache
def curvature_rows(dim, ch):
    literal = helper("_literal_curvature")
    D, G = scalar_directions(dim, ch)
    inv = (
        s.trace(D * G),
        (D * G)[dim - 1, dim - 1],
        D[dim - 1, dim - 1] * G[dim - 1, dim - 1],
        s.trace(D) * s.trace(G),
        s.trace(D) * G[dim - 1, dim - 1],
        s.trace(G) * D[dim - 1, dim - 1],
    )
    g, want = matching.hessians(*inv)
    mapping = {
        literal.A[0]: jets.a,
        literal.A[1]: jets.a * s.diff(jets.a, jets.t),
        literal.A[2]: jets.a * s.diff(jets.a * s.diff(jets.a, jets.t), jets.t),
        literal.p: jets.p,
        literal.Gtime[0]: g[0],
        literal.Gtime[1]: jets.a * g[1],
        literal.Gtime[2]: jets.a**2 * g[2] + jets.a * s.diff(jets.a, jets.t) * g[1],
    }
    return {
        key: s.factor(
            literal.spatial_operator(value).subs(mapping) / jets.a
            - want[key].subs(geometry.d, dim)
        )
        for key, value in literal.literal_curvatures(D, G).items()
    }


@pytest.mark.parametrize("dim", (3, 4, 5, 6))
@pytest.mark.parametrize("ch", density.CHANNELS)
@pytest.mark.parametrize(
    "key", ("R_old", "R_squared", "Ricci_squared", "Riemann_squared")
)
def test_full_exponential_volume_literal_curvature_Hessian(dim, ch, key):
    assert curvature_rows(dim, ch)[key] == 0


@pytest.mark.parametrize(
    "dimension,transfer", (("3", "13"), ("3.125", "29"), ("2.875", "7"))
)
def test_all_350_endpoint_coefficients_from_independent_four_order_WKB(
    dimension, transfer
):
    h = helper("_full_reference")
    with mp.workdps(135):
        d, timeval, P = mp.mpf(dimension), mp.mpf(".25"), mp.mpf(transfer)
        angle = mp.mpf(".37")
        radius = mp.mpf("1e-7") / (1 + P / 1000)
        points = 24
        phases = [mp.exp(2j * mp.pi * i / points) for i in range(points)]
        values = [
            h.mode_endpoint_rows(radius * z, timeval, angle, P, d) for z in phases
        ]
        keys, fn = h.complete_dimensional_endpoint_function()
        expected = fn(d, timeval, angle, P, 1000)
        assert len(keys) == 350
        for (key, degree), target in zip(keys, expected):
            actual = mp.im(
                mp.fsum(values[i][key] * phases[i] ** (-degree) for i in range(points))
                / points
                / radius**degree
            )
            assert abs(actual - target) < mp.mpf("1e-32") * (1 + abs(target)), (
                key,
                degree,
                actual,
                target,
            )


@pytest.mark.parametrize("ch", density.CHANNELS)
def test_legacy_orientation_fails_actual_scalar_covariant_pole(ch):
    inv, g, pole = matching.fixed_pole()
    row = density.physical_channel_row(ch)
    cov = pole.subs(dict(zip(inv, row))).subs(geometry.d, 3)
    entries = density.spatial_difference()
    legacy = sum(
        (-1) ** j * entries[ch, j, r, 4 - j].subs(geometry.d, 3) * g[r]
        for r in range(3)
        for j in range(r, 5)
    )
    wanted = {
        "trace_trace": 2 * s.diff(jets.a, jets.t, 2) * jets.p**2 * g[0],
        "trace_scalar": 2 * s.diff(jets.a, jets.t) * jets.p**2 * g[1],
        "scalar_trace": -2
        * jets.p**2
        * (s.diff(jets.a, jets.t, 2) * g[0] + s.diff(jets.a, jets.t) * g[1]),
    }[ch]
    assert s.factor(cov - 16 * legacy - wanted) == 0
    assert wanted != 0


@pytest.mark.parametrize("index", range(6))
def test_local_finite_proper_time_formal_transpose_by_invariant(index):
    inv, _g, _ell, rows = matching.finite_input()
    transpose = [v.xreplace({inv[4]: inv[5], inv[5]: inv[4]}) for v in rows]
    residuals = (
        rows[2] - transpose[2],
        rows[1] + transpose[1] - 2 * s.diff(transpose[2], jets.t),
        rows[0]
        - transpose[0]
        + s.diff(transpose[1], jets.t)
        - s.diff(transpose[2], jets.t, 2),
    )
    assert all(s.factor(s.diff(v, inv[index])) == 0 for v in residuals)


def test_corrected_finite_is_not_frozen_evanescent_value():
    from p8_vacuum_affine_dimensional_spatial_symbol import matching as old

    inv, _, _, rows = matching.finite_input()
    difference = s.factor(
        rows[0].subs(dict.fromkeys(inv[3:], 0)) - old.finite_input()[3][0]
    )
    assert (
        s.factor(
            difference
            - jets.p**2
            * (3 * jets.t**2 + 1)
            * (-31 * inv[0] + 30 * inv[1])
            / (420 * s.pi**2)
        )
        == 0
    )
    assert difference != 0


def test_deleting_full_volume_changes_scalar_curvature_pole():
    literal = helper("_literal_curvature")
    D = G = s.eye(3)
    # The spatial Einstein coefficient itself carries the trace-volume terms.
    rows = literal.literal_curvatures(D, G)
    full = literal.spatial_operator(rows["R_old"])
    assert full != 0
    # A tracefree-only formula substituted on I/I misses this physical value.
    expected = literal.A[0] ** 2 * literal.p**2 * literal.Gtime[0]
    assert s.factor(full - expected) == 0
    wrong = -(literal.A[0] ** 2) * literal.p**2 * literal.Gtime[0] / 2
    assert s.factor(full - wrong) != 0


@pytest.mark.parametrize(
    "fn", (geometry.data, jets.phase_data, jets.data, density.data, matching.data)
)
def test_complete_packet_residuals_and_gates(fn):
    packet = fn()
    for key, value in packet["checks"].items():
        values = list(value) if isinstance(value, s.MatrixBase) else [value]
        assert all(s.cancel(v) == 0 for v in values), key
    assert all(bool(value) for value in packet["gates"].values())


@pytest.mark.parametrize("fixture", range(2))
@pytest.mark.parametrize(
    "key", ("R_old", "R_squared", "Ricci_squared", "Riemann_squared")
)
def test_noncommuting_nonzero_trace_literal_curvatures(fixture, key):
    literal = helper("_literal_curvature")
    Ds = (
        s.Matrix([[1, 2, -1], [2, 3, 4], [-1, 4, 5]]),
        s.Matrix([[2, 1, 3], [1, -1, 2], [3, 2, 4]]),
    )
    Gs = (
        s.Matrix([[3, -1, 2], [-1, 2, 1], [2, 1, -2]]),
        s.Matrix([[1, -2, 1], [-2, 5, 3], [1, 3, 2]]),
    )
    D, G = Ds[fixture], Gs[fixture]
    assert D * G != G * D
    inv = (
        s.trace(D * G),
        (D * G)[2, 2],
        D[2, 2] * G[2, 2],
        s.trace(D) * s.trace(G),
        s.trace(D) * G[2, 2],
        s.trace(G) * D[2, 2],
    )
    g, want = matching.hessians(*inv)
    mapping = {
        literal.A[0]: jets.a,
        literal.A[1]: jets.a * s.diff(jets.a, jets.t),
        literal.A[2]: jets.a * s.diff(jets.a * s.diff(jets.a, jets.t), jets.t),
        literal.p: jets.p,
        literal.Gtime[0]: g[0],
        literal.Gtime[1]: jets.a * g[1],
        literal.Gtime[2]: jets.a**2 * g[2] + jets.a * s.diff(jets.a, jets.t) * g[1],
    }
    rows = noncommuting_curvature_fixture(fixture)
    assert (
        s.factor(
            literal.spatial_operator(rows[key]).subs(mapping) / jets.a
            - want[key].subs(geometry.d, 3)
        )
        == 0
    )


@cache
def noncommuting_curvature_fixture(fixture):
    Ds = (
        s.Matrix([[1, 2, -1], [2, 3, 4], [-1, 4, 5]]),
        s.Matrix([[2, 1, 3], [1, -1, 2], [3, 2, 4]]),
    )
    Gs = (
        s.Matrix([[3, -1, 2], [-1, 2, 1], [2, 1, -2]]),
        s.Matrix([[1, -2, 1], [-2, 5, 3], [1, 3, 2]]),
    )
    return helper("_literal_curvature").literal_curvatures(Ds[fixture], Gs[fixture])


from p8_vacuum_affine_ordered_scalar_symbol import audit


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_every_exact_corrected_scalar_identity(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, s.MatrixBase) else [value]
    assert all(v == 0 for v in entries), name


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=lambda v: v if isinstance(v, str) else None
)
def test_every_out_of_scope_or_false_closure_input_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_exact_four_historical_withdrawals_and_unchanged_primitive_frontier():
    before = audit.previous.matching()
    after = audit.matching()
    changed = {a["id"] for a, b in zip(before, after) if a != b}
    assert changed == set(audit.CORRECTIONS)
    assert len(after) == len(before) + 1
    assert audit.frontier() == audit.previous.frontier()
    for row in after:
        if row["id"] in audit.CORRECTIONS:
            assert row["status"] == audit.CORRECTIONS[row["id"]]


def test_all_gates_and_explicit_nonclosure():
    assert all(v is True for v in audit.gates().values())
    assert "withdrawn" in audit.observable()
    assert "NOT_FULL_CURRENT_REASSEMBLY" in audit.ITEM["status"]
    assert audit.rejected_inputs() == 152


def test_local_finite_norm_is_not_a_full_scalar_current_bound():
    norm = matching.finite_coefficient_norm()
    assert norm["sum"] == s.Rational(66786194445073, 464486400) < 10**6
    assert "not the complete" in matching.data()["norm_scope"]
