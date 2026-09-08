"""Independent contractions and controls for the two named curvature terms."""
from fractions import Fraction

import pytest
import sympy as sp
from p8_affine import connection
from p8_affine_kinetic import curvature


@pytest.mark.parametrize("name", tuple(curvature.checks()))
def test_exact_checks(name):
    value = curvature.checks()[name]
    assert all(entry == 0 for entry in value) if isinstance(value, sp.MatrixBase) else value == 0


@pytest.mark.parametrize("name,expected", (("axial", -3), ("stf", 2)))
def test_direct_trace_form(name, expected):
    vector = curvature.sectors()[name]
    matrices = [sp.Matrix([[vector[connection.index(a, b, i)] for b in range(4)]
                           for a in range(4)]) for i in range(1, 4)]
    assert sum(sp.trace(matrix*matrix) for matrix in matrices)/2 == expected
    assert all(sp.trace(matrix) == 0 for matrix in matrices)


@pytest.mark.parametrize("value", (1, -1, Fraction(3, 7), sp.Rational(-8, 3)))
def test_both_overall_signs(value):
    result = curvature.sign_control(value)
    assert result["one_strictly_negative"]
    assert result["coefficient_product"] < 0


@pytest.mark.parametrize("value", (0, True, False, 1.0, sp.Float(1), "1", sp.I,
                                   sp.oo, -sp.oo, sp.nan, sp.sqrt(2)))
def test_reject_wrong_rank_or_inexact_input(value):
    with pytest.raises((TypeError, ValueError)):
        curvature.sign_control(value)
