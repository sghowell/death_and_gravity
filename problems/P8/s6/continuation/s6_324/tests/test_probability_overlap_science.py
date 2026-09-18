"""Independent literal-residue, proper-face, overlap and scope checks."""

from itertools import permutations, product

import pytest
import sympy as s
from p8_vacuum_affine_three_real_probability_overlap import (
    audit,
    current,
    faces,
    measure,
    projection,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_each_exact_residual(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_each_written_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_false_scope_and_invalid_inputs_rejected(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("scale", (s.Rational(3, 7), s.Rational(11, 5)))
def test_literal_EH3_residue_at_independent_rotated_hard_null_direction(scale):
    e = current.e
    # A non-axis-aligned rational rotation with all10 independent soft entries
    # specialized nontrivially. The hard root remains conserved, not TT.
    O = s.ImmutableMatrix(
        [
            [1, 0, 0, 0],
            [0, s.Rational(3, 5), 0, s.Rational(4, 5)],
            [0, 0, 1, 0],
            [0, -s.Rational(4, 5), 0, s.Rational(3, 5)],
        ]
    )
    assert O.T * e.ETA * O == e.ETA
    q = scale * O * s.ImmutableMatrix([1, 0, 0, 1])
    B0 = s.ImmutableMatrix([[0, 0, 0, 0], [0, 2, 3, 0], [0, 3, -2, 0], [0, 0, 0, 0]])
    H0 = s.ImmutableMatrix([[2, 3, 5, 2], [3, 7, 11, 3], [5, 11, 13, 5], [2, 3, 5, 2]])
    A = s.ImmutableMatrix(
        [[2, 3, 5, 7], [3, 11, 13, 17], [5, 13, 19, 23], [7, 17, 23, 29]]
    )
    B, H = O * B0 * O.T, O * H0 * O.T
    assert H * e.ETA * q == e.VECTOR_ZERO and e.tr(H) != 0
    propagated = e.ETA * H * e.ETA - e.ETA * e.tr(H) / 2
    actual = -e.cubic_gravity((A, B, propagated), (e.VECTOR_ZERO, q, -q))
    expected = 2 * (q.T * A * q)[0] * s.trace(H * B)
    assert expected != 0 and s.factor(actual - expected) == 0
    assert s.factor(actual - expected / 2) != 0


def test_scalar_residue_retains_off_shell_trace_remainder():
    e = current.e
    p = s.ImmutableMatrix([2, 1, 0, 0])
    A = s.ImmutableMatrix(s.diag(1, 2, 3, 5))
    shell = current.dot(p, p) - 1
    assert shell != 0 and e.tr(A) != 0
    vertex = e.scalar_vertex(p, -p, (A,), 1)
    assert -vertex == 2 * (p.T * A * p)[0] - e.tr(A) * shell
    assert -vertex != 2 * (p.T * A * p)[0]


def test_both_remaining_null_emitters_required_by_exact_conservation():
    row = current.calibration()
    ps = row["whole_reduced_two_real_state"]
    qb, qc = row["whole_two_marked_null_legs"]
    ray = s.ImmutableMatrix([1, 1, 0, 0])
    tensor = lambda legs: sum((p * p.T / current.dot(p, ray) for p in legs), s.zeros(4))
    assert tensor((*ps, qb, qc)) * current.e.ETA * ray == current.e.VECTOR_ZERO
    assert tensor(ps) * current.e.ETA * ray == -qb - qc
    assert tensor((*ps, qb)) * current.e.ETA * ray == -qc
    assert qc != current.e.VECTOR_ZERO
    good, wrong = row["whole_correct_and_massive_only_soft_coefficients"]
    assert good != wrong and good != 0 and wrong != 0


def test_all_three_calibrations_use_original_source_and_exact_recoil():
    row = current.calibration()
    assert row["whole_original_parameters"] == source.original_parameters()
    samples = row["whole_three_exact_hierarchical_calibrations"]
    assert [r["h"] for r in samples] == [
        s.Rational(1, 1000),
        s.Rational(1, 10000),
        s.Rational(1, 100000),
    ]
    assert all(
        r["tree_count"] == 5116 and r["error_over_h_ceiling"] == 2 for r in samples
    )
    assert all(
        0 < sum(r["weights"]) < s.Rational(1, 8) and 0 < r["rho"] <= 1 for r in samples
    )
    assert (
        samples[2]["weights"][0] < samples[1]["weights"][0] < samples[0]["weights"][0]
    )
    assert "only at the three recorded samples" in row["whole_calibration_boundary"]


@pytest.mark.parametrize(
    "r", (s.Rational(1, 10**40), s.Rational(2, 3), s.Integer(2), s.Integer(10) ** 40)
)
def test_null_angular_current_bound_with_plus_cross_complex_tensor(r):
    # Rotate the null ray in its transverse plane; both tensor components enter.
    ct, st = s.Rational(3, 5), s.Rational(4, 5)
    ni = s.ImmutableMatrix([1, 0, 0, 1])
    nj = s.ImmutableMatrix(
        [
            1,
            2 * r * ct / (1 + r * r),
            2 * r * st / (1 + r * r),
            (1 - r * r) / (1 + r * r),
        ]
    )
    A = s.ImmutableMatrix(
        [
            [0, 0, 0, 0],
            [0, s.Rational(1, 2), s.I / 2, 0],
            [0, s.I / 2, -s.Rational(1, 2), 0],
            [0, 0, 0, 0],
        ]
    )
    coefficient = s.factor((nj.T * A * nj)[0] / current.dot(nj, ni))
    assert sum(s.expand(v * s.conjugate(v)) for v in A) == 1
    assert current.dot(nj, nj) == 0 and A * ni == current.e.VECTOR_ZERO
    assert s.expand(coefficient * s.conjugate(coefficient)) <= 4


def test_current_disc_and_derived_face_bounds_use_no_full_amplitude_disc():
    row = faces.data()
    assert row["whole_current_absolute_radius"] == s.Rational(1, 10**6)
    assert row["whole_current_derivative_caps"] == (129, 10**9, 10**15)
    assert (
        "No such absolute disc is asserted for G2 or G3"
        in row["whole_kinematic_domain"]
    )
    B, e, C2 = (
        faces.pair.LINEAR_SLOPE,
        faces.pair.RADIUS,
        faces.pair.REMAINDER_COEFFICIENT,
    )
    D2 = faces.tree.coefficient_envelope(2)
    L = (129 * B / e + 10**9 * D2) / s.sqrt(source.KAPPA)
    Q = (129 * C2 + 2 * 10**9 * B / e + 10**15 * D2) / s.sqrt(source.KAPPA)
    assert row["whole_derived_face_derivative_coefficients"] == (L, Q)
    assert 0 < L < faces.AXIS_CAP and 0 < Q < faces.PAIR_CAP
    assert B / e**2 < C2 and e / 8 < faces.CURRENT_RADIUS / 2


def rectangle(expr, variables):
    return s.expand(
        sum(
            (-1) ** (len(variables) - sum(bits))
            * expr.subs({v: 0 for v, bit in zip(variables, bits) if not bit})
            for bits in product((0, 1), repeat=len(variables))
        )
    )


def test_probability_transfer_against_independent_complex_polynomial_rectangle():
    a, b, c = variables = s.symbols("a b c", real=True)
    g = [
        1 + s.I,
        (2 - s.I) * a,
        (-1 + 3 * s.I) * b,
        (4 + s.I) * a * b,
        (3 + 2 * s.I) * c,
        (-2 + s.I) * a * c,
        (1 - 4 * s.I) * b * c,
        (2 + 5 * s.I) * a * b * c,
    ]
    G, H = sum(g), sum(g[:7])
    square = lambda z: s.expand(z * s.conjugate(z))
    signed = square(G) - square(H)
    transfer = sum(
        g[i] * s.conjugate(g[j]) for i, j in product(range(7), repeat=2) if i | j == 7
    )
    assert s.expand(rectangle(square(G), variables) - signed - transfer) == 0
    assert s.expand(transfer) != 0
    assert s.expand(rectangle(G, variables) - g[7]) == 0


@pytest.mark.parametrize("perm", tuple(permutations(range(3))))
def test_all_union_terms_and_face_symmetry(perm):
    remap = lambda mask: sum(((mask >> i) & 1) << perm[i] for i in range(3))
    proper = set(projection.data()["whole_proper_overlap_pairs"])
    assert len(proper) == 12
    assert {(remap(a), remap(b)) for a, b in proper} == proper
    all_pairs = [(a, b) for a, b in product(range(8), repeat=2) if a | b == 7]
    assert len(all_pairs) == 27 and sum(7 in pair for pair in all_pairs) == 15
    assert projection.data()["whole_term_inventory"] == {
        "all": 27,
        "contains_full": 15,
        "proper": 12,
    }


def test_pair_pair_overlap_is_nonzero_and_signed_measure_need_not_be_positive():
    a, b, c = variables = s.symbols("a b c", real=True)
    G = a * b + a * c
    assert rectangle(G, variables) == 0
    assert rectangle(G**2, variables) == 2 * a * a * b * c
    assert rectangle((1 + a * b - c) ** 2, variables) == -2 * a * b * c
    assert projection.data()["gates"]["dropping_pair_pair_terms_changes_result"] is True


def test_entropy_product_majorant_and_original_phase_space():
    a, b, c, x = s.symbols("a b c x", positive=True)
    product_majorant = 4 * a * s.sqrt(b * c) / (a * b * c)
    assert s.integrate(product_majorant, (a, 0, x), (b, 0, x), (c, 0, x)) == 16 * x * x
    angular = 8 * (4 * s.pi) ** 3 / (s.factorial(3) * 2**3 * (2 * s.pi) ** 9)
    assert s.factor(angular - 1 / (48 * s.pi**6)) == 0
    L, Q = faces.AXIS_CAP, faces.PAIR_CAP
    assert (
        measure.data()["whole_transfer_TV_coefficient"] == (L * Q + 2 * Q * Q) / s.pi**6
    )
    assert L * Q + 2 * Q * Q < s.Rational(1, 10**1424)


@pytest.mark.parametrize(
    "x", (s.Rational(1, 10**50), s.Rational(1, 128), s.Rational(1, 8))
)
def test_total_variation_bounds_and_optional_reference(x):
    bound = measure.total_variation_upper(x)
    assert bound == 2 * s.Rational(1, 10**1424) * x * x + s.Rational(1, 10**1340) * x**4
    assert 0 < bound < x
    assert measure.relative_reference_upper(x) == 2 * bound / s.sqrt(x)


@pytest.mark.parametrize(
    "bad",
    (
        True,
        False,
        None,
        "1",
        1.0,
        s.Float(1),
        s.I,
        0,
        -1,
        s.oo,
        s.Symbol("free"),
        s.Rational(1, 4),
    ),
)
def test_cutoff_guards_reject_inexact_or_out_of_domain_inputs(bad):
    with pytest.raises((TypeError, ValueError)):
        measure.total_variation_upper(bad)
    with pytest.raises((TypeError, ValueError)):
        measure.relative_reference_upper(bad)


def test_radiative_state_mismatch_keeps_other_energy_logarithm():
    a, b, eps, x = s.symbols("a b eps x", positive=True)
    density = (a + b) / (a * b)
    integral = s.integrate(density, (a, eps, x), (b, eps, x))
    assert s.simplify(integral - 2 * (x - eps) * s.log(x / eps)) == 0
    assert s.limit(integral, eps, 0, dir="+") == s.oo
    mismatch = measure.data()["whole_same_state_matching_remainder"]
    assert mismatch != 0 and len(mismatch.free_symbols) == 4
    assert (
        measure.data()["gates"]["no_virtual_or_integrated_counterterm_chosen"] is True
    )


def test_original_parameters_and_all_historical_frontiers_preserved():
    pars = source.original_parameters()
    assert pars["heavy"] == s.Integer(10) ** 200 / 512 + 2
    assert pars["kappa"] == s.Integer(10) ** 800
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 180
    assert len(audit.qualifications()) == 6 and audit.rejected_inputs() == 81
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.controls()) == 8
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    for word in ("real-virtual", "all-N", "evanescent", "Regge", "bounce", "state"):
        assert word in audit.observable()["not_established"]
