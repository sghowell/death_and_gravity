"""Physical contours, full moments, literal insertions and strict scope tests."""

import pytest
import sympy as s
from p8_vacuum_affine_physical_loop_contours import (
    audit,
    calibration,
    contour,
    insertions,
    moments,
    parameters,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_residual(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_all_rejected_scopes(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("kind", parameters.KINDS)
@pytest.mark.parametrize("v", (-12, 0, s.Rational(4, 3), s.Rational(45, 8), 16))
def test_complete_invariant_domain(kind, v):
    virtualities = (0, 2) if kind == "triangle" else (0, 1, 2, 1)
    result = parameters.domain(kind, v, virtualities, 0)
    assert result[0] == kind and result[1] == v
    assert result[2] == virtualities and result[4] == source.HEAVY_MASS2


@pytest.mark.parametrize("kind", parameters.KINDS)
@pytest.mark.parametrize("degree", (0, 1, 2))
def test_all_order_Cauchy_rule(kind, degree):
    alpha = (degree,) + (0,) * (2 if kind == "triangle" else 5)
    assert moments.derivative_bound(kind, alpha) == s.factorial(
        degree
    ) * 4096**degree * moments.master_bound(kind)
    assert (
        moments.derivative_bound(kind, alpha, 1)
        == s.factorial(degree)
        * 4096**degree
        * moments.master_bound(kind)
        * 128
        / source.HEAVY_MASS2
    )


@pytest.mark.parametrize("p,r", tuple(moments.MOMENT_BUDGETS))
def test_all_needed_weighted_moments(p, r):
    coefficient, power = moments.MOMENT_BUDGETS[(p, r)]
    assert (
        moments.moment_bound(p, r) == s.Integer(coefficient) / source.HEAVY_MASS2**power
    )


@pytest.mark.parametrize("kind", parameters.KINDS)
@pytest.mark.parametrize("line_type", ("light", "heavy"))
def test_all_light_and_heavy_line_bounds(kind, line_type):
    value = insertions.line_bound(kind, line_type)
    assert value > 0 and not value.has(s.Float)
    assert insertions.line_bound(kind, line_type, numerator_norm2=0) == 0


def test_nonzero_literal_TT_numerator_not_Ward_guess():
    assert insertions.data()["whole_literal_TT_numerator"] != s.zeros(4)
    assert "xa*d_gamma" in insertions.data()["whole_full_TT_line_integral"]
    assert "Gamma" in insertions.data()["whole_full_TT_line_integral"]


def test_both_contours_have_fixed_endpoints_and_correct_signs():
    assert contour.light(s.S.Zero) == 0 and contour.light(s.S.One) == 1
    assert contour.heavy(s.S.Zero) == 0 and contour.heavy(s.S.One) == 1
    assert s.im(contour.light(s.Rational(1, 5))) > 0
    assert s.im(contour.heavy(s.Rational(1, 5))) < 0


def test_every_actual_calibration_count():
    assert calibration.data()["whole_exact_calibration_counts"] == {
        "homotopy": 350,
        "gap": 175,
        "external_virtualities": 48,
        "external_pairs": 72,
        "cyclic_routings": 99,
        "routing_domain_gates": 513,
        "actual_TT_insertions": 66,
    }


def test_no_principal_sheet_claim():
    assert (
        "NOT a single first-sheet disk" in contour.data()["whole_continuation_boundary"]
    )
    assert contour.RHO == s.Rational(1, 4096)
    assert contour.LIGHT_GAP == s.Rational(1, 32)


def test_original_parameters_and_frontiers_unchanged():
    n, g = source.HEAVY_MASS2, source.CUBIC
    assert source.CONTACT == -g * g * (3 / (n - 2) - 2 / (n - 2) ** 2)
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 197
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "curved" in audit.observable()["not_established"]
    assert "internal-graviton" in audit.observable()["not_established"]


def test_module_counts_and_uncompleted_amplitude_scope():
    assert len(parameters.data()["checks"]) == 7
    assert len(contour.data()["checks"]) == 19
    assert len(moments.data()["checks"]) == 7
    assert len(insertions.data()["checks"]) == 23
    assert len(calibration.data()["checks"]) == 78
    assert "not" in insertions.data()["whole_scope"].lower()
    assert "finite counterterms" in audit.observable()["not_established"]
