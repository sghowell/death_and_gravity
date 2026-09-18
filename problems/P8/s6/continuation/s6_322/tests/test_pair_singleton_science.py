"""Independent complete-pair, off-shell-core, anisotropic and scope tests."""

from itertools import product

import pytest
import sympy as s
from p8_vacuum_affine_pair_singleton_subtraction import (
    algebra,
    audit,
    hard,
    hierarchy,
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


@pytest.mark.parametrize("bits", tuple(product((0, 1), repeat=2)))
def test_complete_pair_hierarchy_all_basis_cases(bits):
    row = hierarchy.data()["whole_four_basis_pair_certificates"][
        "".join(map(str, bits))
    ]
    assert row["bounds"] == hierarchy.EXPECTED[bits]
    assert len(row["coefficient_certificates"]) == 24
    for certificate in row["coefficient_certificates"].values():
        if "zero" in certificate:
            assert certificate == {"zero": True}
        else:
            assert certificate["denominator_identity"] == 0
            assert (
                certificate["budget"]
                == certificate["coefficient_l1"] / certificate["positive_prefactor"]
            )


def test_pair_hierarchy_derivative_maxima_and_offdiagonal_norm():
    assert hierarchy.data()["whole_pair_hierarchy_caps"] == (530, 1566, 586, 3104)
    assert len(hierarchy.data()["whole_four_basis_pair_certificates"]) == 4
    # The frozen coefficients contain twice each off-diagonal entry.
    a, b, c, d, e, f = s.symbols("a b c d e f", real=True)
    T = s.Matrix([[a, b, c], [b, d, e], [c, e, f]])
    assert (
        s.expand(
            s.trace(T.T * T)
            - (a * a + 2 * b * b + 2 * c * c + d * d + 2 * e * e + f * f)
        )
        == 0
    )


def test_coefficient_majorant_and_unsupported_denominator_controls():
    z, r = s.symbols("z r", nonnegative=True)
    bound, record = hierarchy.coefficient_budget((3 + 2 * z * r) / (1 + r * r), z, r)
    assert bound == 5 and record["angular_power"] == 1
    for expression in (
        1 / r,
        1 / (1 + r),
        r**3 / (1 + r * r),
        1 / (1 - r * r),
        1 / (z * (1 + r * r)),
        s.Symbol("free") / (1 + r * r),
    ):
        with pytest.raises(ValueError):
            hierarchy.coefficient_budget(expression, z, r)
    assert hierarchy.coefficient_budget(0, z, r) == (0, {"zero": True})


def test_generic_offshell_vertex_and_shift_are_not_TT_shortcuts():
    p = s.ImmutableMatrix([s.Rational(5, 4), s.Rational(3, 4), 0, 0])
    Q = s.ImmutableMatrix([s.Rational(1, 16), s.Rational(1, 32), 0, 0])
    A = s.ImmutableMatrix(s.diag(0, 1, 2, 3))
    k = s.ImmutableMatrix([s.Rational(1, 64), 0, s.Rational(1, 64), 0])
    assert algebra.dot(p, p) == 1 and algebra.dot(Q, Q) > 0
    assert algebra.numerator(p, A, Q) != (p.T * A * p)[0]
    no_trace = 2 * (p.T * A * k)[0] + (k.T * A * k)[0] + (k.T * A * Q)[0]
    exact = algebra.numerator(p + k, A, Q) - algebra.numerator(p, A, Q)
    assert s.factor(exact - no_trace) != 0
    B = s.ImmutableMatrix(s.diag(0, 1, -1, 0))
    q = s.ImmutableMatrix([s.Rational(1, 16), 0, 0, s.Rational(1, 16)])
    assert algebra.numerator(p, B, q) == (p.T * B * p)[0]
    assert all(algebra.identities()["gates"].values())


def test_original387_and140_and_independent_temporal_class():
    row = algebra.calibration()
    assert row["whole_marked_calibration_count"] == 2
    assert row["whole_temporal_class_inventory"] == 387
    assert row["gates"]["actual_complete_pair_has_nonzero_trace"]
    assert row["gates"]["actual_complete_pair_not_transverse"]
    assert row["checks"]["collapsed_core_equals_independent_temporal_class"] == 0


def test_three_anisotropic_points_and_phase_retained():
    row = hard.complex_calibration()
    records = row["whole_three_exact_anisotropic_points"]
    assert len(records) == 3 and records[0]["shift"] == 0
    assert records[1]["shift"] == s.I / 10**18
    assert records[2]["shift"] == -s.I / 10**18
    assert records[1]["weights"] == tuple(s.conjugate(w) for w in records[2]["weights"])
    for number in (1, 2):
        assert row["gates"][f"{number}_outside_old_relative_c_disc"]
        assert all(
            row["gates"][f"{number}_inside_anisotropic_disc{i}"] for i in range(3)
        )
        assert row["gates"][f"{number}_phased_triangle_budget"]


def test_exact_cauchy_product_budget_and_unresolved_complement():
    row = hard.budgets()
    assert row["whole_anisotropic_radius_eta"] == s.Rational(1, 10**13)
    assert 0 < row["whole_pair_product_derivative_coefficient"] < s.Rational(1, 10**723)
    assert 0 < row["whole_retained188_derivative_coefficient"] < s.Rational(1, 10**738)
    assert row["whole_unsubtracted_three_singleton_terms"] == 3767
    assert row["whole_component_budgets"]["matter_final"] == 13632146432
    assert row["whole_component_budgets"]["gravity_final"] == s.Rational(
        1952093286852638820978919582900252803544196317184, 9
    )


def test_entropy_faces_and_permutation_sum_without_extra_three():
    a, b, _c = s.symbols("a b c", positive=True)
    I = (a + b) * s.log(a + b) - a * s.log(a) - b * s.log(b)
    assert s.diff(I, a, b) == 1 / (a + b)
    assert s.limit(I, a, 0) == 0 and s.limit(I, b, 0) == 0
    sample = (s.Rational(1, 128), s.Rational(1, 64), s.Rational(1, 32))
    pieces = sum(
        hard.class_rectangle(*args)
        for args in (
            sample,
            (sample[0], sample[2], sample[1]),
            (sample[1], sample[2], sample[0]),
        )
    )
    assert s.expand(hard.all_nonsingleton_rectangle(*sample) - 2 * pieces) == 0


@pytest.mark.parametrize(
    "bad", (True, False, None, "1", 1.0, s.Float(1), s.I, 0, -1, s.oo, s.Symbol("free"))
)
def test_energy_guards_precede_logarithms(bad):
    with pytest.raises((TypeError, ValueError)):
        hard.class_rectangle(bad, s.Rational(1, 64), s.Rational(1, 64))


def test_total_energy_guard():
    with pytest.raises(ValueError):
        hard.all_nonsingleton_rectangle(
            s.Rational(1, 8), s.Rational(1, 8), s.Rational(1, 8)
        )


def test_original_parameter_record_without_diagnostic_replacement():
    record = source.original_parameters()
    assert record["heavy"] == s.Integer(10) ** 200 / 512 + 2
    assert record["kappa"] == s.Integer(10) ** 800
    assert record["cubic"] == source.CUBIC and record["contact"] == source.CONTACT


def test_all_inherited_frontiers_remain_exactly_unchanged():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 178
    assert len(audit.qualifications()) == 6 and audit.rejected_inputs() == 75
    assert audit.matching()[:-1] == audit.previous.matching()
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    for word in ("3767", "5116", "real-virtual", "evanescent", "Regge", "bounce"):
        assert word in audit.observable()["not_established"]
    assert len(audit.controls()) == 8
