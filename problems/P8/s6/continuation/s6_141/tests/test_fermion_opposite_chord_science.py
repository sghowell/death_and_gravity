"""Independent graph, routing, angular, Dirac and radial checks."""

import itertools

import mpmath as mp
import pytest
import sympy as sp
from p8_vacuum_fermion_opposite_chord import angular, audit, catalog, domain, tail
from p8_vacuum_fermion_proper_references import dirac


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


def test_all_labelled_permutations_independent():
    words = set()
    for w in itertools.permutations((1, 2, 3, 4, "B0", "B1")):
        i = w.index(1)
        w = w[i:] + w[:i]
        ends = [j for j, x in enumerate(w) if isinstance(x, str)]
        if abs(ends[0] - ends[1]) == 3:
            words.add(tuple("B" if isinstance(x, str) else x for x in w))
    assert words == set(catalog.selected_words())
    assert len(words) == 12


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
def test_independent_all_proper_1PI_edge_subgraphs(word):
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
            found.append(tuple(ids))
    expected = [c["edges"] for c in catalog.simple_cycles(word)]
    assert set(found) == set(expected)
    assert len(found) == 3


@pytest.mark.parametrize("word", catalog.selected_words())
def test_every_routing_with_independent_integer_momentum_conservation(word):
    momenta = {
        1: sp.Matrix([1, 2, 3, 4]),
        2: sp.Matrix([2, 3, -1, 5]),
        3: sp.Matrix([-1, 4, 7, 2]),
    }
    momenta[4] = -sum(momenta.values(), sp.zeros(4, 1))
    r = domain.routing(word)
    A = sum((momenta[i] for i in r["arc_l_external"]), sp.zeros(4, 1))
    B = sum((momenta[i] for i in r["arc_q_external"]), sp.zeros(4, 1))
    q, l = sp.Matrix([3, 2, 1, -1]), sp.Matrix([7, 4, -2, 1])
    assert A + B == sp.zeros(4, 1)
    assert (q + A) - (l + A) == q - l
    assert len(r["arc_l_external"]) == len(r["arc_q_external"]) == 2


@pytest.mark.parametrize(
    "x,y,b", tuple(itertools.product((1, 4, 10000), (1, 3, 100), (0, 1)))
)
def test_angular_integral_independent_quadrature(x, y, b):
    with mp.workdps(40):
        xx, yy = mp.mpf(x), mp.mpf(y)

        def f(theta):
            if x == y and b == 0:
                return mp.cos(theta / 2) ** 2 / xx
            return mp.sin(theta) ** 2 / (
                xx + yy + b - 2 * mp.sqrt(xx * yy) * mp.cos(theta)
            )

        actual = 2 * mp.quad(f, [0, mp.pi / 4, mp.pi]) / mp.pi
        upper = 1 / max(xx, yy)
        if b == 0:
            assert mp.almosteq(actual, upper)
        else:
            assert 0 < actual < upper


@pytest.mark.parametrize("n", (1, 2, 3, 4, 8, 20))
def test_unbounded_radial_integral_independently(n):
    with mp.workdps(50):
        if n == 1:
            # Map the slowly decaying infinite tail to a smooth compact integral.
            actual = 8 * mp.quad(lambda u: (1 - u * u) * u ** (n - 1), [0, 1])
        else:
            actual = 4 * mp.quad(
                lambda t: t / (1 + t) ** (2 + mp.mpf(n) / 2), [0, 1, mp.inf]
            )
        assert mp.almosteq(actual, mp.mpf(16) / (n * (n + 2)))


@pytest.mark.parametrize(
    "qscale,lscale", ((0, 0), (0, 100), (100, 0), (1, 1), (10, 1000))
)
def test_joint_Dirac_resolvent_bound_for_asymmetric_loop_momenta(qscale, lscale):
    gamma = dirac.gamma_matrices()
    m = 720
    for scale in (qscale, lscale):
        q = sp.Matrix([m * scale, 0, 0, 0])
        S = m * m + (q.T * q)[0]
        delta = sp.sqrt(S) * sp.Matrix([sp.I, 1, 0, 0]) / 100
        D = m * sp.eye(4) + sp.I * sum(
            ((q[i] + delta[i]) * gamma[i] for i in range(4)), sp.zeros(4)
        )
        inv = D.inv()
        # Positive leading principal minors prove the squared operator-norm
        # bound directly, without substituting a Frobenius-norm estimate.
        gap = (4 * sp.eye(4) / S - inv.conjugate().T * inv).applyfunc(sp.simplify)
        assert all(bool(sp.simplify(gap[:j, :j].det()) > 0) for j in range(1, 5))
        assert sp.Abs(delta[0]) + sp.Abs(delta[1]) < sp.sqrt(S) / 20


@pytest.mark.parametrize("seed", range(10))
def test_S4_degree_two_only_diagonal_and_offdiagonal_Gram_orbits(seed):
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
def test_no_missing_primitive_row_allowed(i):
    rows = audit.frontier()
    del rows[i]
    with pytest.raises(ValueError):
        audit.validate_frontier(rows)


def test_exact_counts_and_partial_frontier():
    assert len(audit.residuals()) == 70
    assert audit.scalar_entry_count() == 88
    assert len(audit.gates()) == 29
    assert audit.rejected_inputs() == 66
    assert len(audit.controls()) == 9
    rows = audit.frontier()
    assert sum(r["status"] == "UNEVALUATED" for r in rows) == 5
    assert (
        sum(
            r["status"] == "PARTIAL_SELF_ENERGY_AND_OPPOSITE_CHORD_BOUNDS" for r in rows
        )
        == 2
    )
    assert "not" in domain.data()["scope"]
    assert "1/max" in angular.data()["normalized_angular_average"]


def test_scalar_and_gauge_pieces_are_separate():
    result = tail.enclosure(720, 1, 0, 144)
    assert result["gauge_opposite_chord_b2_absolute_upper"] == 0
    assert (
        result["combined_opposite_chord_b2_absolute_upper"]
        == result["scalar_opposite_chord_b2_absolute_upper"]
    )
    assert result["cumulatively_bounded_words_per_sector"] == 36
    assert result["remaining_vertex_words_per_sector"] == 24
