"""Independent word, complex-domain, kernel and convergent-tail tests."""

import itertools

import mpmath as mp
import pytest
import sympy as sp
from p8_vacuum_fermion_self_energy_chord import audit, catalog, domain, kernel, tail


@pytest.mark.parametrize("name", tuple(audit.residuals()))
def test_every_exact_scalar_and_matrix_entry(name):
    value = audit.residuals()[name]
    entries = list(value) if isinstance(value, sp.MatrixBase) else [value]
    assert all(v == 0 for v in entries)


@pytest.mark.parametrize("name", tuple(audit.gates()))
def test_every_proof_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=lambda v: v if isinstance(v, str) else None
)
def test_every_unsupported_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_independent_all_labelled_six_vertex_permutations():
    words = set()
    for perm in itertools.permutations((1, 2, 3, 4, "B0", "B1")):
        word = tuple("B" if str(v).startswith("B") else v for v in perm)
        i = word.index(1)
        words.add(word[i:] + word[:i])
    assert words == {r["word"] for r in catalog.catalog()}
    assert len(words) == 60


def test_independent_three_shortest_chord_distances():
    counts = {1: 0, 2: 0, 3: 0}
    for row in catalog.catalog():
        ends = [i for i, v in enumerate(row["word"]) if v == "B"]
        d = (ends[1] - ends[0]) % 6
        counts[min(d, 6 - d)] += 1
    assert counts == {1: 24, 2: 24, 3: 12}
    assert len(set(catalog.box_placements())) == 24


@pytest.mark.parametrize(
    "n,b,x",
    tuple(
        itertools.product(
            (2, 10, 360), (0, 1), (sp.Rational(1, 100), sp.Rational(1, 2), 1)
        )
    ),
)
def test_complex_parameter_gap_for_large_internal_momenta(n, b, x):
    m = 2 * n
    q = n * n - 1
    scale = n * n + 1
    assert scale * scale == m * m + q * q
    delta = sp.Rational(scale, 20) * (-1 + sp.I)
    p = q + delta
    Delta = x * m * m + (1 - x) * b + x * (1 - x) * p * p
    assert sp.re(Delta) >= x * m * m / 2
    assert sp.Abs(Delta) <= 2 * (m * m + q * q)


@pytest.mark.parametrize("m,t", tuple(itertools.product((2, 720), (0, 4, 100))))
def test_direct_whole_MS_kernel_operator_norm_in_a_complex_direction(m, t):
    with mp.workdps(40):
        mm = mp.mpf(m)
        q = mm * t
        scale = mp.sqrt(mm * mm + q * q)
        p = q + scale * (-1 + 1j) / 20
        Js0 = mp.quad(
            lambda x: mp.log((x * mm * mm + (1 - x) + x * (1 - x) * p * p) / (mm * mm)),
            [0, 1],
        )
        Js1 = mp.quad(
            lambda x: (
                (1 - x)
                * mp.log((x * mm * mm + (1 - x) + x * (1 - x) * p * p) / (mm * mm))
            ),
            [0, 1],
        )
        Jg0 = mp.quad(lambda x: mp.log(x + x * (1 - x) * p * p / (mm * mm)), [0, 1])
        Jg1 = mp.quad(
            lambda x: (1 - x) * mp.log(x + x * (1 - x) * p * p / (mm * mm)), [0, 1]
        )
        Cf = mp.mpf(4) / 3
        mass = mm * (Js0 + Cf * (-2 - 4 * Jg0)) / 144
        kinetic = (-Js1 + Cf * (-1 - 2 * Jg1)) / 144
        norm = max(abs(mass + 1j * p * kinetic), abs(mass - 1j * p * kinetic))
        upper = (3 + 12 * Cf) * scale * (mp.log(1 + t * t) + 6) / 144
        assert norm < upper


@pytest.mark.parametrize("R", (2, 3, 10, 100))
def test_degree_four_geometric_tail_bound_independently(R):
    exact = sp.Rational(1, R**3 * (R - 1))
    assert exact <= sp.Rational(2, R**4)
    partial = sum(sp.Rational(1, R**n) for n in range(4, 40))
    assert partial < exact


def test_radial_logarithmic_integral_by_independent_quadrature():
    with mp.workdps(50):
        value = mp.quad(
            lambda t: t * (mp.log(1 + t) + 6) / (1 + t) ** 4, [0, 1, mp.inf]
        )
        assert mp.almosteq(value, mp.mpf(41) / 36)


@pytest.mark.parametrize("seed", range(10))
def test_S4_degree_two_Gram_polynomial_has_only_two_orbits(seed):
    pairs = [(i, i) for i in range(4)] + list(itertools.combinations(range(4), 2))
    symbols = sp.symbols("gram0:10")
    lookup = dict(zip(pairs, symbols))
    i, j = pairs[seed]
    total = 0
    for perm in itertools.permutations(range(4)):
        total += lookup[tuple(sorted((perm[i], perm[j])))]
    expected = sum(symbols[:4]) / 4 if i == j else sum(symbols[4:]) / 6
    assert sp.expand(total / 24 - expected) == 0


@pytest.mark.parametrize("index", range(9))
def test_no_primitive_row_can_be_removed(index):
    rows = audit.frontier()
    del rows[index]
    with pytest.raises(ValueError):
        audit.validate_frontier(rows)


def test_exactly_two_new_partial_rows_not_complete_rows():
    rows = audit.frontier()
    partial = {
        r["id"] for r in rows if r["status"] == "PARTIAL_SELF_ENERGY_CHORD_BOUND"
    }
    assert partial == {"scalar_Phi4_W0_F4", "gauge_Phi4"}
    assert sum(r["status"] == "UNEVALUATED" for r in rows) == 5


def test_scalar_and_gauge_limits_and_no_double_count():
    d = tail.enclosure(720, 1, 0, 144)
    assert d["gauge_self_energy_chord_b2_absolute_upper"] == 0
    assert (
        d["combined_self_energy_chord_b2_absolute_upper"]
        == d["scalar_self_energy_chord_b2_absolute_upper"]
    )
    assert (
        tail.enclosure(720, 0, 1, 144)["combined_self_energy_chord_b2_absolute_upper"]
        == 0
    )


def test_scope_and_exact_counts():
    assert len(audit.residuals()) == 66
    assert audit.scalar_entry_count() == 84
    assert len(audit.gates()) == 30
    assert audit.rejected_inputs() == 66
    assert len(audit.controls()) == 9
    assert "not" in domain.data()["scope"]
    assert "whole" in kernel.__doc__.lower()
    assert "vertex-chord" in tail.data()["scope"]
