"""Independent graph, routing, angular, Dirac and radial checks."""

import itertools

import mpmath as mp
import numpy as np
import pytest
import sympy as sp
from p8_vacuum_fermion_proper_references import dirac
from p8_vacuum_fermion_vertex_chord import audit, catalog, forest, regions, tail


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_all_exact_identities(name):
    v = audit.residuals()[name]
    assert all(x == 0 for x in (list(v) if isinstance(v, sp.MatrixBase) else [v]))


@pytest.mark.parametrize("name", tuple(audit.gates()))
def test_all_proof_gates(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=lambda v: v if isinstance(v, str) else None
)
def test_invalid_inputs(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def connected(edges):
    if not edges:
        return False
    vertices = {v for e in edges for v in e}
    reached = {next(iter(vertices))}
    while True:
        new = reached | {v for e in edges if reached.intersection(e) for v in e}
        if new == reached:
            return reached == vertices
        reached = new


@pytest.mark.parametrize("word", catalog.selected_words())
def test_every_proper_1PI_edge_subgraph_independently(word):
    ends = [i for i, x in enumerate(word) if x == "B"]
    edges = [(i, (i + 1) % 6) for i in range(6)] + [tuple(ends)]
    found = []
    for mask in range(1, 127):
        ids = [i for i in range(7) if mask & (1 << i)]
        sub = [edges[i] for i in ids]
        vertices = {v for e in sub for v in e}
        if not connected(sub) or len(sub) - len(vertices) + 1 <= 0:
            continue
        if all(
            connected(sub[:i] + sub[i + 1 :])
            and {v for e in sub[:i] + sub[i + 1 :] for v in e} == vertices
            for i in range(len(sub))
        ):
            nf = sum(i < 6 for i in ids)
            found.append((nf, len(ids) - nf, 4 - nf - 2 * (len(ids) - nf)))
    assert sorted(found) == [(2, 1, 0), (4, 1, -2), (6, 0, -2)]


def test_all_labelled_six_vertex_orders_independently():
    words = set()
    for w in itertools.permutations((1, 2, 3, 4, "B0", "B1")):
        i = w.index(1)
        w = w[i:] + w[:i]
        ends = [j for j, x in enumerate(w) if isinstance(x, str)]
        d = abs(ends[0] - ends[1])
        if min(d, 6 - d) == 2:
            words.add(tuple("B" if isinstance(x, str) else x for x in w))
    assert words == set(catalog.selected_words()) == set(catalog.marked_box_vertices())
    assert len(words) == 24


@pytest.mark.parametrize("word", catalog.selected_words())
def test_short_arc_routing_keeps_real_chord(word):
    w = catalog.short_arc_rotation(word)
    assert w[0] == w[2] == "B"
    assert len(w[3:]) == 3
    q, l, p = sp.symbols("q l p")
    assert (q + p) - (l + p) == q - l


@pytest.mark.parametrize(
    "qscale,lscale,direction,b",
    tuple(itertools.product((0, 1, 10), (30, 300), (-1, 0, 1), (0, 1))),
)
def test_high_region_paired_matrix_kernel(qscale, lscale, direction, b):
    gamma = [np.array(g.tolist(), dtype=complex) for g in dirac.gamma_matrices()]
    m = 720.0
    q = np.array([m * qscale, 0, 0, 0], dtype=complex)
    l = np.array(
        [m * lscale * direction, m * lscale * (1 - abs(direction)), 0, 0], dtype=complex
    )
    Sq = m * m + float(np.vdot(q, q).real)
    y = float(np.vdot(l, l).real)
    assert y >= 4 * Sq
    delta1 = np.array([1j, 0, 0, 0]) * np.sqrt(Sq) / 100
    delta2 = np.array([0, 1, 0, 0]) * np.sqrt(Sq) / 100

    def propagator(k):
        D = m * np.eye(4) + 1j * sum(
            (k[i] * gamma[i] for i in range(4)), np.zeros((4, 4), complex)
        )
        return np.linalg.inv(D)

    S0 = propagator(l)
    B1 = 1 / (float(np.vdot(q - l, q - l).real) + b)
    B0 = 1 / (y + b)
    pair = B1 * propagator(l + delta1) @ propagator(l + delta2) - B0 * S0 @ S0
    upper = 32 * np.sqrt(Sq) / y**2.5
    assert np.linalg.norm(pair, 2) < upper
    assert abs(B1 - B0) <= 12 * np.sqrt(Sq) / y**1.5


@pytest.mark.parametrize("x", (sp.Rational(1, 10), 1, 10, 10000))
def test_LOW_raw_inner_integral_with_independent_quadrature(x):
    with mp.workdps(40):
        xx = mp.mpf(str(x.evalf(35) if isinstance(x, sp.Basic) else x))
        S = 1 + xx

        def density(y):
            return xx * y / (max(xx, y) * S * S * (1 + y) * min(S, 1 + y) ** 2)

        actual = mp.quad(density, [0, xx, 4 * S])
        analytic = xx**2 / (2 * S**4) + xx * mp.log((1 + 4 * S) / S) / S**4
        assert mp.almosteq(actual, analytic)


def test_integrated_regional_majorants_independently():
    with mp.workdps(50):
        raw = mp.mpf(1) / 6 + mp.quad(lambda u: u * (1 - u) * mp.log(4 + u), [0, 1])
        sub = mp.quad(lambda u: u * (1 - u) * (mp.log(4 + u) - mp.log(u)), [0, 1])
        assert 0 < raw < mp.mpf(1) / 2
        assert 0 < sub < mp.mpf(17) / 36
        assert mp.almosteq(raw - mp.mpf(1) / 6, sub - mp.mpf(5) / 36)


@pytest.mark.parametrize("x", (0, 1, 10, 10000))
def test_HIGH_radial_moment_independently(x):
    with mp.workdps(45):
        actual = mp.quad(lambda y: y ** (-mp.mpf(3) / 2), [4 * (1 + x), mp.inf])
        assert abs(actual - 1 / mp.sqrt(1 + x)) < mp.mpf("1e-20")


@pytest.mark.parametrize("n", (1, 2, 3, 4, 8))
def test_positive_soft_degree_and_log_moments(n):
    with mp.workdps(40):
        k = mp.mpf(n) / 2
        moment = mp.quad(lambda u: u ** (k - 1) * (1 - u), [0, 1])
        logmoment = mp.quad(lambda u: -(u ** (k - 1)) * (1 - u) * mp.log(u), [0, 1])
        assert abs(moment - 1 / (k * (k + 1))) < mp.mpf("1e-18")
        assert abs(logmoment - (1 / k**2 - 1 / (k + 1) ** 2)) < mp.mpf("1e-17")


@pytest.mark.parametrize("seed", range(10))
def test_S4_degree_two_Gram_orbits(seed):
    pairs = [(i, i) for i in range(4)] + list(itertools.combinations(range(4), 2))
    symbols = sp.symbols("gram0:10")
    lookup = dict(zip(pairs, symbols))
    i, j = pairs[seed]
    total = sum(
        lookup[tuple(sorted((p[i], p[j])))] for p in itertools.permutations(range(4))
    )
    expected = sum(symbols[:4]) / 4 if i == j else sum(symbols[4:]) / 6
    assert sp.expand(total / 24 - expected) == 0


@pytest.mark.parametrize("i", range(9))
def test_no_primitive_row_can_be_dropped(i):
    rows = audit.frontier()
    del rows[i]
    with pytest.raises(ValueError):
        audit.validate_frontier(rows)


def test_direct_gauge_and_scalar_UV_leading_matrices():
    gamma = dirac.gamma_matrices()
    Sleading = -sp.I * gamma[0]
    scalar = Sleading * Sleading
    gauge = -sum((g * scalar * g for g in gamma), sp.zeros(4))
    assert scalar == -sp.eye(4)
    assert gauge == 4 * sp.eye(4)


def test_exact_counts_and_complete_selected_rows():
    assert len(audit.residuals()) == 66
    assert audit.scalar_entry_count() == 81
    assert len(audit.gates()) == 29
    assert audit.rejected_inputs() == 66
    assert len(audit.controls()) == 9
    rows = audit.frontier()
    bounded = {r["id"] for r in rows if r["status"] == "BOUNDED_PAIRED_PRIMITIVE_ROW"}
    assert bounded == {"scalar_Phi4_W0_F4", "gauge_Phi4"}
    assert sum(r["status"] == "UNEVALUATED" for r in rows) == 5
    assert "raw HIGH" in regions.data()["scope"]
    assert "canonical" in forest.data()["scope"]


def test_separate_scalar_gauge_limits_and_full_word_count():
    r = tail.enclosure(720, 1, 0, 144)
    assert r["gauge_vertex_chord_b2_absolute_upper"] == 0
    assert r["cumulatively_bounded_words_per_sector"] == 60
    assert r["remaining_words_in_these_two_primitive_rows"] == 0
    assert (
        tail.enclosure(720, 0, 1, 144)["combined_vertex_chord_b2_absolute_upper"] == 0
    )
