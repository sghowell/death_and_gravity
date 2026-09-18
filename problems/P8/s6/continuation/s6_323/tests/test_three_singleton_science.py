"""Independent generic, original-source, domain and scope checks."""

from functools import cache
from itertools import product

import pytest
import sympy as s
from p8_vacuum_affine_three_singleton_subtraction import (
    analytic,
    audit,
    forest,
    grouping,
    line,
    measure,
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


def test_generic_connected_numerator_against_independent_rational_recursion():
    g = line.generic_polynomials()
    values = {s.Symbol(name): s.Integer(i % 7 + 1) for i, name in enumerate(g.names)}
    values.update(
        dict(
            zip(
                s.symbols("a b c"),
                (s.Rational(1, 32), s.Rational(1, 64), s.Rational(1, 128)),
            )
        )
    )
    D = {mask: expr.as_expr().subs(values) for mask, expr in g.D.items()}
    assert all(value > 0 for value in D.values())

    @cache
    def current(mask):
        if not mask:
            return s.S.One
        out = s.S.Zero
        block = mask
        while block:
            rest = mask ^ block
            out += (
                g.vertex(block, rest).as_expr().subs(values) * current(rest) / D[mask]
            )
            block = (block - 1) & mask
        return out

    connected = current(7) - sum(
        current(7 ^ (1 << i)) * current(1 << i) for i in range(3)
    )
    connected += 2 * current(1) * current(2) * current(4)
    denominator = s.prod(D.values())
    assert connected * denominator == g.P.as_expr().subs(values)
    omitted_cumulant = connected - 2 * current(1) * current(2) * current(4)
    assert omitted_cumulant * denominator != g.P.as_expr().subs(values)


def test_generic_support_and_original_vertex_mapping():
    g = line.generic_polynomials()
    assert len(g.names) == 49 and len(g.P) == 62
    assert (min(map(sum, g.P)), max(map(sum, g.P))) == (6, 10)
    assert g.Q * s.prod(g.w) == g.P
    assert all(min(powers) >= 1 for powers in g.P)
    assert line.physical_mapping()["whole_exact_mapping_fixture_count"] == 16
    assert len(set(line.vector("v1"))) == len(line.vector("v1")) == 4
    assert all(line.generic_data()["gates"].values())


def test_original_forest_and_count_oracles_not_overwritten_by_loop_indices():
    row = forest.data()
    assert row["whole_central_graph_counts"] == (7, 19, 67, 307)
    assert row["whole_forest_graph_counts_by_central_leaves"] == {
        0: 1288,
        1: 1368,
        2: 804,
        3: 307,
    }
    assert len(row["whole_soft_canonical_power_stripped_parts"]) == 8
    assert all(
        value != 0
        for value in row["whole_soft_canonical_power_stripped_parts"].values()
    )
    x = s.Symbol("x")
    central_egf = s.exp(x) + 6 * s.exp(2 * x) / (2 - s.exp(x))
    for order, expected in enumerate((7, 19, 67, 307)):
        assert s.diff(central_egf, x, order).subs(x, 0) == expected
        assert row["checks"][f"central_EGF_coefficient{order}"] == 0
    assert row["checks"]["whole_forest_equals_original_root"] == 0
    assert row["checks"]["whole_eight_class_sum_equals_original_root"] == 0


@pytest.mark.parametrize(
    "number,sides,mask",
    (
        (2, (0, 1), 15),
        (3, (0, 0, 1), 15),
        (3, (1, 1, 1), 7),
    ),
)
def test_endpoint_owner_factors_reconstruct_whole_coefficients(number, sides, mask):
    row = grouping.group_case(number, sides, mask)
    assert row["squarefree_support"]
    for record in row["records"]:
        assert record["residual"] == record["cross_multiplication"] == 0
        assert (
            s.expand(
                record["quotient"] * record["untouched_factor"] - record["coefficient"]
            )
            == 0
        )
    data = grouping.data()
    assert data["whole_counts_by_external_multiplicity"] == {
        2: (64, 186),
        3: (128, 590),
    }
    assert all(data["whole_negative_controls"].values())


@pytest.mark.parametrize(
    "number,sides,mask",
    (
        (True, (0, 0), 0),
        (2.0, (0, 0), 0),
        (s.Integer(2), (0, 0), 0),
        (1, (0,), 0),
        (4, (0, 0, 0, 0), 0),
        (2, [0, 0], 0),
        (2, (0,), 0),
        (2, (0, True), 0),
        (2, (0, 2), 0),
        (2, (0, 0), True),
        (2, (0, 0), -1),
        (2, (0, 0), 16),
    ),
)
def test_grouping_guards_precede_cached_ring_construction(number, sides, mask):
    grouping.coefficient_ring(2)
    with pytest.raises((TypeError, ValueError)):
        grouping.group_case(number, sides, mask)


@pytest.mark.parametrize("number", (True, False, 2.0, s.Integer(2), None, "2", 0, 4))
def test_ring_entry_point_rejects_boolean_and_unsupported_aliases(number):
    grouping.coefficient_ring(2)
    with pytest.raises((TypeError, ValueError)):
        grouping.coefficient_ring(number)


def test_actual_scalar_pair_pole_is_not_declared_holomorphic():
    row = line.global_tube_control()
    assert row["connected_pair_pole_residue_nonzero"] is True
    assert row["pole_inside_naive_global_W_disc"] is True
    geometry = analytic.geometry()
    assert "not asserted for connected" in geometry["whole_positive_center_tube"]
    assert analytic.ETA == s.Rational(1, 10**13)
    assert analytic.SIGMA == s.Rational(1, 10**14)


def test_literal_origin_geometry_and_full_derived_budget():
    row = analytic.budgets()
    parts = row["whole_eight_core_rectangle_coefficients"]
    assert len(parts) == 8 and all(value > 0 for value in parts.values())
    assert sum(parts.values()) == row["whole_core_rectangle_coefficient"]
    assert (
        row["whole_full5116_rectangle_coefficient"]
        == sum(parts.values()) + row["whole_retained1349_rectangle_coefficient"]
    )
    assert 0 < row["whole_full5116_rectangle_coefficient"] < s.Rational(1, 10**670)
    assert row["whole_central_majorant_coefficients"]["matter"] == (
        s.Rational(5, 2),
        3088,
        18923696,
        174400047744,
    )
    group = row["whole_grouped_kernel_budgets"]
    assert group["common_C"] == 10**50 * (source.HEAVY_MASS2**2 + 1)
    assert group["second_remainder_B3"] == 10**4 * group["common_C"] / analytic.SIGMA**3
    assert group["first_remainder_B1"] == 10**4 * group["common_C"] / analytic.SIGMA**2


@pytest.mark.parametrize(
    "t,W",
    tuple(
        product(
            (s.Rational(1, 10**40), s.Rational(1, 64), s.S.One, s.Integer(3)),
            (s.Rational(1, 10**60), s.Rational(1, 10**16), s.Rational(1, 8)),
        )
    ),
)
def test_two_regime_envelopes_with_exact_hierarchical_parameters(t, W):
    sigma = analytic.SIGMA
    # Divide out the positive common C. These check the arithmetic, not
    # the analytic hypotheses proved in notes/grouping.md and analytic.md.
    if W <= sigma * t / 4:
        r3, r1 = 64 * W**2 / sigma**2, 16 * W / sigma
    else:
        assert t < 4 * W / sigma
        raw3 = 2 * t**2 * (t + W)
        raw1 = 2 * t**2
        assert raw3 <= 160 * W**3 / sigma**3
        assert raw1 <= 32 * W**2 / sigma**2
        r3, r1 = 256 * W**2 / sigma**3, 64 * W / sigma**2
    assert 4 * r3 < 10**4 * W**2 / sigma**3
    assert 4 * r1 < 10**4 * W / sigma**2


def rectangle(expression, variables):
    return s.expand(
        sum(
            (-1) ** (len(variables) - sum(bits))
            * expression.subs({v: 0 for v, bit in zip(variables, bits) if not bit})
            for bits in product((0, 1), repeat=len(variables))
        )
    )


def test_origin_polynomial_removed_but_probability_not_automatically_subtracted():
    a, b, c = variables = s.symbols("a b c", real=True)
    assert rectangle(1 + a + 2 * b + 3 * c, variables) == 0
    assert rectangle(a * b * c, variables) == a * b * c
    counterexample = 1 + a * b + c
    assert rectangle(counterexample, variables) == 0
    assert rectangle(counterexample**2, variables) == 2 * a * b * c


def test_entropy_compatible_faces_and_all_three_central_choices():
    a, b, c = s.symbols("a b c", positive=True)
    I = lambda u, v: (u + v) * s.log(u + v) - u * s.log(u) - v * s.log(v)
    J = c * I(a, b) + b * I(a, c) + a * I(b, c)
    assert s.diff(I(a, b), a, b) == 1 / (a + b)
    assert all(s.limit(J, variable, 0) == 0 for variable in (a, b, c))
    sample = (s.Rational(1, 128), s.Rational(1, 64), s.Rational(1, 32))
    assert (
        s.expand(
            analytic.rectangle_envelope(*sample) * 10**670
            - J.subs(dict(zip((a, b, c), sample)))
        )
        == 0
    )


@pytest.mark.parametrize(
    "bad", (True, False, None, "1", 1.0, s.Float(1), s.I, 0, -1, s.oo, s.Symbol("free"))
)
def test_energy_guards_precede_logarithms(bad):
    with pytest.raises((TypeError, ValueError)):
        analytic.rectangle_envelope(bad, s.Rational(1, 64), s.Rational(1, 64))


def test_original_parameters_total_energy_and_scope():
    assert source.original_parameters()["heavy"] == s.Integer(10) ** 200 / 512 + 2
    assert source.original_parameters()["kappa"] == s.Integer(10) ** 800
    with pytest.raises(ValueError):
        analytic.rectangle_envelope(*([s.Rational(1, 8)] * 3))
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 179
    assert len(audit.qualifications()) == 6 and audit.rejected_inputs() == 79
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.controls()) == 8
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    for word in (
        "probability",
        "real-virtual",
        "all-N",
        "evanescent",
        "Regge",
        "bounce",
    ):
        assert word in audit.observable()["not_established"]


def test_signed_measure_projectors_preserve_every_face_and_retain_overlap():
    faces, P, R = measure.projectors()
    assert P * P == P and R * R == R and P + R == s.eye(8)
    assert P * R == R * P == s.zeros(8)
    assert list(P[7, :]) == [1, -1, -1, 1, -1, 1, 1, 0]
    for E in faces:
        assert E * P == P * E == E and E * R == R * E == s.zeros(8)
    wrong = P - faces[0] * faces[1]
    assert any(E * wrong != E for E in faces)


def test_signed_density_retains_interference_and_is_not_probability_rectangle():
    B, R = 1 + s.I, 2 - s.I
    square = lambda z: s.expand(z * s.conjugate(z))
    assert square(B + R) - square(B) == 7
    assert square(R) == 5
    assert s.expand(B * s.conjugate(R) + s.conjugate(B) * R) == 2
    assert all(measure.identities()["gates"].values())


def test_complete_signed_measure_baseline_and_exact_phase_space():
    row = measure.integration()
    assert row[
        "whole_baseline_constant"
    ] == 7 * measure.finite_tree.coefficient_envelope(3)
    assert 0 < row["whole_baseline_constant"] < s.Rational(1, 10**754)
    angular = 8 * (4 * s.pi) ** 3 / (s.factorial(3) * 2**3 * (2 * s.pi) ** 9)
    assert s.factor(angular - row["whole_angular_polarization_symmetry_factor"]) == 0
    assert angular == 1 / (48 * s.pi**6)
    assert row["whole_amplitude_rectangle_constant"] == s.Rational(1, 10**670)


@pytest.mark.parametrize(
    "x", (s.Rational(1, 10**50), s.Rational(1, 128), s.Rational(1, 8))
)
def test_signed_measure_cutoff_bounds_and_reference_exponents(x):
    assert (
        measure.total_variation_upper(x)
        == s.Rational(1, 10**1424) * x**2 + s.Rational(1, 10**1340) * x**4
    )
    assert measure.relative_reference_upper(x) == 2 * s.Rational(
        1, 10**1424
    ) * x ** s.Rational(3, 2) + 2 * s.Rational(1, 10**1340) * x ** s.Rational(7, 2)
    assert 0 < measure.total_variation_upper(x) < 1
    assert not measure.total_variation_upper(x).has(s.Float)


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
def test_signed_measure_cutoff_guards(bad):
    for call in (measure.total_variation_upper, measure.relative_reference_upper):
        with pytest.raises((TypeError, ValueError)):
            call(bad)
