"""Independent reconstruction, exact hierarchy and original scope checks."""

from itertools import product

import pytest
import sympy as s
from p8_vacuum_affine_three_soft_current_subtraction import (
    audit,
    hard,
    source,
    subtraction,
)
from sympy.polys.rings import ring

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_all_exact_residuals(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_all_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_original_scope_rejections(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "sector,bits",
    [(sector, bits) for sector in (0, 1) for bits in product((0, 1), repeat=3)],
)
def test_each_complete_sector_TT_hierarchy_and_faces(sector, bits):
    row = subtraction.case_data(sector, bits)
    assert row["graph_count"] == 4 and len(row["checks"]) == 66
    assert len(row["coefficient_certificates"]) == 42
    assert len(row["compatible_axis_certificates"]) == 18
    assert row["third_mixed_matrix_l1_bound"] < 721
    for derivative, value in row["lower_matrix_l1_bounds"].items():
        assert value >= 0 and value <= 10**6
        assert derivative in ("100", "010", "001", "110", "101", "011")
    assert all(
        r["sequential_limits_agree"] and r["joint_axis_support"]
        for r in row["compatible_axis_certificates"]
    )


@pytest.mark.parametrize(
    "bad", (True, False, 1.0, s.Integer(1), s.Rational(1, 2), None, -1, 2, "0", [], {})
)
def test_sector_validation_precedes_cache(bad):
    subtraction.case_data(0, (0, 0, 0))
    with pytest.raises(ValueError):
        subtraction.case_data(bad, (0, 0, 0))


@pytest.mark.parametrize(
    "bad",
    (
        None,
        "000",
        [0, 0, 0],
        (),
        (0, 0),
        (0, 0, 0, 0),
        (False, 0, 0),
        (0, True, 0),
        (1.0, 0, 0),
        (s.Integer(0), 0, 0),
        (2, 0, 0),
        (-1, 0, 0),
    ),
)
def test_basis_validation_precedes_cache(bad):
    with pytest.raises(ValueError):
        subtraction.case_data(0, bad)


def test_positive_coefficient_AMGM_certificate_and_controls():
    R, x, y = ring("x,y", s.QQ)
    value, midpoints = subtraction.coefficient_bound(
        x * y, x * x + y * y, subtraction.shifts(2)
    )
    assert value == s.Rational(1, 2) and midpoints == 1
    direct, _ = subtraction.coefficient_bound(
        3 * x * x + 2 * y * y, x * x + y * y, subtraction.shifts(2)
    )
    assert direct == 3
    for numerator, denominator in ((x, R.one), (x, R.zero), (x, R.one - x)):
        with pytest.raises(ValueError):
            subtraction.coefficient_bound(numerator, denominator, subtraction.shifts(2))


@pytest.mark.parametrize(
    "numerator,denominator",
    (("a", "1"), ("1+a", "1"), ("a", "a+b"), ("a*a", "a*a*(b+1)"), ("a+b", "a+c")),
)
def test_valuation_faces_against_independent_symbolic_limits(numerator, denominator):
    R, *_generators = ring("a,b,c", s.QQ)
    variables = dict(zip(("a", "b", "c"), R.symbols))
    P = R.from_expr(s.sympify(numerator, locals=variables))
    Q = R.from_expr(s.sympify(denominator, locals=variables))
    A, B = subtraction.face(P, Q, 0)
    target = s.limit(P.as_expr() / Q.as_expr(), R.symbols[0], 0, dir="+")
    assert s.factor(A.as_expr() / B.as_expr() - target) == 0


def test_divergent_face_is_not_accepted():
    R, a, _b = ring("a,b", s.QQ)
    with pytest.raises(ValueError):
        subtraction.face(R.one, a, 0)
    assert subtraction.face(R.zero, a, 0) == (R.zero, R.one)


def test_squarefree_Taylor_inverse_has_no_unexamined_series_remainder():
    ea, eb, ec = source.VARIABLES
    polynomial = 2 + ea - 3 * eb + 5 * ec + ea * eb + eb * ec
    inverse = source.inverse(polynomial)
    assert source.jet(polynomial * inverse) == 1
    assert source.jet(ea**2 + eb**2 + ec**2) == 0
    with pytest.raises(ValueError):
        source.inverse(ea)


def test_TT_frame_symbolic_geometry_and_complex_basis_normalization():
    x, y = s.symbols("x y", real=True)
    U, V = source.frame(x, y)
    n = source.direction(x, y)
    assert s.factor(U.dot(U) - 1) == 0 and s.factor(V.dot(V) - 1) == 0
    assert s.factor(U.dot(V)) == s.factor(U.dot(n)) == s.factor(V.dot(n)) == 0
    plus = U * U.T - V * V.T
    cross = U * V.T + V * U.T
    assert s.factor(s.trace(plus.T * plus) - 2) == 0
    assert s.factor(s.trace(cross.T * cross) - 2) == 0
    assert s.factor(s.trace(plus.T * cross)) == 0
    made = (plus + s.I * cross) / 2
    assert s.simplify(s.trace(s.conjugate(made).T * made)) == 1


def test_literal_source_and_taylor_baseline_counts():
    assert len(source.taylor_calibration()) == 96
    assert len(source.angular_data()["checks"]) == 72
    assert all(hard.polynomial_calibration()["gates"].values())


def test_one_block_original47_and188_calibrations():
    row = hard.core_calibration()
    assert row["whole_original_parameter_calibration_points"] == 4
    assert row["whole_independent_spatial_root_calibrations"] == 24
    assert row["whole_temporal_class_inventory"] == 188
    assert row["whole_core_EGF_linear_coefficients"] == (
        44048,
        5566277620447838208000000,
    )


def test_exact_one_block_product_budget_and_unresolved_complement():
    row = hard.analytic_data()
    assert 0 < row["whole_product_rule_coefficient"] < s.Rational(1, 10**738)
    assert row["whole_remaining_temporal_terms"] == 4928
    assert row["gates"]["remaining4928_temporal_terms_not_declared_subtracted"]
    assert row["gates"]["no_full_all_N_real_virtual_Regge_or_original_P8_closure"]


def test_all_sixteen_cases_and_exact_scientific_row_inventory():
    row = subtraction.data()
    assert len(row["whole_reconstructed_sixteen_cases"]) == 16
    assert len(row["checks"]) == 1056
    assert row["whole_unit_complex_TT_hierarchy_caps"] == (
        4800000000000000000000,
        10**6,
        10**6,
        721,
    )
    assert all(row["gates"].values())


def test_independent_research_basis_bounds_reconstructed_exactly():
    first = subtraction.case_data(0, (0, 0, 0))
    second = subtraction.case_data(1, (0, 0, 0))
    assert first["third_mixed_matrix_l1_bound"] == s.Rational(987329876, 6089853)
    assert second["third_mixed_matrix_l1_bound"] == s.Rational(67469, 105)
    assert first["lower_matrix_l1_bounds"] == {
        "100": s.Rational(1587, 7),
        "010": s.Integer(264),
        "001": s.Integer(371),
        "110": s.Rational(11814044, 161805),
        "101": s.Rational(78703, 630),
        "011": s.Rational(16412020, 187131),
    }


def test_primitive_matching_and_historical_boundaries_retained():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 177
    assert len(audit.qualifications()) == 6 and audit.rejected_inputs() == 73
    assert audit.matching()[:-1] == audit.previous.matching()
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
    scope = audit.observable()
    for word in ("4928", "5116", "real-virtual", "evanescent", "Regge", "bounce"):
        assert word in scope["not_established"]
