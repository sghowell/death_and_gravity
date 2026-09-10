"""Independent graph, overlap, radial, complex matrix and OS Cauchy tests."""

import copy
import itertools

import mpmath as mp
import pytest
import sympy as sp
from p8_vacuum_fermion_quadratic_forests import (
    audit,
    calibration,
    catalog,
    pole,
)


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_exact_identity(name, value):
    if isinstance(value, sp.MatrixBase):
        assert value == sp.zeros(*value.shape), name
    else:
        assert value == 0, name


@pytest.mark.parametrize("name,value", list(audit.gates().items()))
def test_proof_gate(name, value):
    assert value is True, name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_unsupported_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("row", catalog.word_rows())
def test_all_connected_bridgeless_proper_subgraphs(row):
    a, b = row["boson_vertices"]
    edges = [(i, (i + 1) % 4) for i in range(4)] + [(a, b)]

    def connected(ids, vertices):
        if not vertices:
            return False
        seen = {min(vertices)}
        for _ in range(4):
            for i in ids:
                u, v = edges[i]
                if u in seen or v in seen:
                    seen |= {u, v}
        return seen == vertices

    found = []
    for mask in range(1, 31):
        ids = tuple(i for i in range(5) if mask & (1 << i))
        verts = {v for i in ids for v in edges[i]}
        if not connected(ids, verts):
            continue
        if all(connected(tuple(j for j in ids if j != i), verts) for i in ids):
            found.append(ids)
    assert set(found) == {tuple(c["edges"]) for c in catalog.graph_cycles(row)}
    for first, second in itertools.combinations(found, 2):
        assert set(first) & set(second)


def test_labelled_determinant_and_boson_exchange_quotient():
    words = set()
    for tail in itertools.permutations(("P2", "B1", "B2")):
        words.add(tuple("B" if w.startswith("B") else w for w in ("P1",) + tail))
    assert words == {r["word"] for r in catalog.word_rows()}
    for w in words:
        changed = tuple("P2" if x == "P1" else "P1" if x == "P2" else x for x in w)
        j = changed.index("P1")
        assert changed[j:] + changed[:j] in words


@pytest.mark.parametrize("m", (1, 2, 720))
def test_exact_radial_self_and_log_moments(m):
    with mp.workdps(40):
        m = mp.mpf(m)
        I0 = mp.quad(lambda x: x / (m * m + x) ** 3, [0, m * m, mp.inf])
        Ilog = mp.quad(
            lambda x: x * mp.log(1 + x / (m * m)) / (m * m + x) ** 3, [0, m * m, mp.inf]
        )
        assert abs(I0 * m * m - mp.mpf(".5")) < mp.mpf("1e-32")
        assert abs(Ilog * m * m - mp.mpf(".75")) < mp.mpf("1e-32")


@pytest.mark.parametrize("m", (1, 3, 720))
def test_middle_and_high_moments_separately(m):
    with mp.workdps(40):
        m = mp.mpf(m)
        mm = m * m
        raw = mp.quad(
            lambda x: x / (mm + x) ** 3 * mp.log((mm + 4 * (mm + x)) / (mm + x)),
            [0, mm, mp.inf],
        )
        comp = mp.quad(
            lambda y: y / (mm + y) ** 3 * mp.log((mm + y) / max(mm, y / 4)),
            [0, 4 * mm, mp.inf],
        )
        highcomp = mp.quad(
            lambda y: y / (mm + y) ** 3 * mp.log(y / (4 * mm)), [4 * mm, mp.inf]
        )
        assert 0 < raw < 1 / mm
        assert 0 < comp < 1 / mm
        assert 0 < highcomp < mp.mpf(".75") / mm
        raw_compact = mp.quad(lambda u: (1 - u) * mp.log(4 + u), [0, 1])
        assert abs(raw * mm - raw_compact) < mp.mpf("1e-30")
        # The complementary middle/high terms partition the full log moment;
        # neither equals the raw middle integral separately.
        assert abs((comp + highcomp) * mm - mp.mpf(".75")) < mp.mpf("1e-30")


@pytest.mark.parametrize("y", (".001", ".1", "1", "4", "5", "100", "1000000"))
def test_complementary_middle_domain_and_log_cap(y):
    with mp.workdps(40):
        y = mp.mpf(y)
        lower = max(0, y / 4 - 1)
        integral = mp.quad(lambda x: 1 / (1 + x), [lower, y])
        assert abs(integral - mp.log((1 + y) / max(1, y / 4))) < mp.mpf("1e-30")
        assert integral <= mp.log(5) < 2


def gamma_matrices():
    sigma = (
        mp.matrix([[0, 1], [1, 0]]),
        mp.matrix([[0, -1j], [1j, 0]]),
        mp.matrix([[1, 0], [0, -1]]),
    )
    g = []
    for s in sigma:
        g.append(
            mp.matrix(
                [
                    [0, 0, -1j * s[0, 0], -1j * s[0, 1]],
                    [0, 0, -1j * s[1, 0], -1j * s[1, 1]],
                    [1j * s[0, 0], 1j * s[0, 1], 0, 0],
                    [1j * s[1, 0], 1j * s[1, 1], 0, 0],
                ]
            )
        )
    g.append(mp.matrix([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]))
    return g


@pytest.mark.parametrize("sector", ("scalar", "gauge"))
@pytest.mark.parametrize("regime", ("middle", "high"))
@pytest.mark.parametrize("angle", (".1", "1.7", "4.2"))
def test_full_paired_matrix_tail_with_complementary_subtraction(sector, regime, angle):
    with mp.workdps(45):
        m = mp.mpf(720)
        q = mp.matrix([10, 20, 30, 40])
        l = (
            mp.matrix([80, 50, 20, -10])
            if regime == "middle"
            else mp.matrix([2400, 100, -50, 0])
        )
        x = sum(a * a for a in q)
        y = sum(a * a for a in l)
        Sq = m * m + x
        Sl = m * m + y
        z = 1 + 2 * mp.exp(1j * mp.mpf(angle))
        p = mp.matrix([1j * mp.sqrt(z), 0, 0, 0])
        gamma = gamma_matrices()
        I = mp.eye(4)

        def S(v):
            return (
                m * I - 1j * sum((gamma[i] * v[i] for i in range(4)), mp.zeros(4))
            ) / (m * m + sum(a * a for a in v))

        trace = lambda A: sum(A[i, i] for i in range(4))
        cap = lambda A, B: (
            trace(A * B)
            if sector == "scalar"
            else sum(trace(g * A * g * B) for g in gamma)
        )
        b = 1 if sector == "scalar" else 0
        chord = sum(a * a for a in q - l) + b

        def paired(t):
            Dq = S(q + t * p) * S(q)
            Dl = S(l) * S(l + t * p)
            return (
                cap(Dq, Dl) / chord
                - cap(Dq, S(l) ** 2) / (y + b)
                - cap(S(q) ** 2, Dl) / (x + b)
            )

        tail = paired(1) - sum(
            mp.diff(paired, 0, n) / mp.factorial(n) for n in range(4)
        )
        base = 32 * 360**4
        if regime == "middle":
            bound = base * (
                4 / (Sq**3 * Sl * chord)
                + 1 / (Sq**3 * Sl * (y + b))
                + 1 / (Sl**3 * Sq * (x + b))
            )
        else:
            assert y >= 4 * Sq
            bound = base * (
                32 / (Sq ** mp.mpf("2.5") * y ** mp.mpf("2.5"))
                + 1 / (Sl**3 * Sq * (x + b))
            )
        if sector == "gauge":
            bound *= 4
        assert abs(tail) < bound
        # Trace reversal validates the symmetric half-domain for the complete word.
        raw = cap(S(q + p) * S(q), S(l) * S(l + p))
        swapped = cap(S(l + p) * S(l), S(q) * S(q + p))
        assert abs(raw - swapped) < mp.mpf("1e-40")


@pytest.mark.parametrize("power", range(5))
def test_arbitrary_regulator_local_polynomial_is_annihilated(power):
    s, e = sp.symbols("s epsilon")
    a, b, c, d = sp.symbols("a b c d")
    f = (a + b * s) / e**power + (c + d * s) * e
    assert sp.expand(f - f.subs(s, 1) - (s - 1) * sp.diff(f, s).subs(s, 1)) == 0


@pytest.mark.parametrize("location", (4, 8, 20, 4 + 3j))
@pytest.mark.parametrize("angle", (".1", "2.3", "4.7"))
def test_second_Cauchy_OS_bound_on_independent_analytic_function(location, angle):
    with mp.workdps(40):
        location = mp.mpc(location)
        s = 1 + mp.mpf(".5") * mp.exp(1j * mp.mpf(angle))
        f = lambda z: 1 / (location - z)
        remainder = (f(s) - f(1) - (s - 1) / (location - 1) ** 2) / (s - 1) ** 2
        E = 1 / (abs(location - 1) - 2)
        assert abs(remainder) <= E / 3


@pytest.mark.parametrize("i", range(9))
def test_each_frontier_row_is_required(i):
    rows = audit.frontier()
    del rows[i]
    with pytest.raises(ValueError):
        audit.validate_frontier(rows)


@pytest.mark.parametrize(
    "field", ("finite_MS_mass", "finite_MS_slope", "canonical_two_loop_pole")
)
def test_partial_scope_cannot_be_promoted(field):
    rows = copy.deepcopy(audit.frontier())
    row = next(r for r in rows if r["id"] == audit.TARGETS[0])
    row[field] = "COMPLETE"
    with pytest.raises(ValueError):
        audit.validate_frontier(rows)


def test_actual_nonlocal_bound_and_uncomputed_local_scope():
    d = calibration.data()["enclosure"]
    assert 0 < d["nonlocal_on_shell_divided_remainder_upper"] < sp.Rational(1, 10**799)
    assert "uncomputed" in d["scope"]
    assert (
        pole.enclosure(720, 0, 1, 144)["nonlocal_on_shell_divided_remainder_upper"] == 0
    )
    assert sum(r["status"] == "UNEVALUATED" for r in audit.frontier()) == 2
