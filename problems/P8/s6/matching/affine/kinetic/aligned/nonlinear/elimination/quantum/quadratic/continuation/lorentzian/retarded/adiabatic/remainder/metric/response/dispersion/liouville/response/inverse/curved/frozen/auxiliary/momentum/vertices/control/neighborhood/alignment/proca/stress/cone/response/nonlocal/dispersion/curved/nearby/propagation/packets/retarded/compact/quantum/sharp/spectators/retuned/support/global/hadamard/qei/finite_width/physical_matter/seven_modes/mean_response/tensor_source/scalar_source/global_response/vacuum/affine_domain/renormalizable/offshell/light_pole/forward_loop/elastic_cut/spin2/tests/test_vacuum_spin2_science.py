"""Independent Ward, transfer, tensor, symmetry and exact calibration tests."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_vacuum_spin2 import audit, calibration, pole, stress, triangles


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_stress_triangle_Ward_tensor_or_calibration_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_continuous_bound_or_explicit_scientific_scope_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args",
    calibration.bad_cases(),
    ids=[r[0] for r in calibration.bad_cases()],
)
def test_reject_unsupported_exact_calibration(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_control_and_unclosed_obligation(name, value):
    assert bool(value), name


@pytest.mark.parametrize("radius", (1, 2, 3, Fraction(1, 2), sp.Rational(7, 2)))
def test_exact_transfer_disc_remainders(radius):
    d = calibration.point(10**400, radius)
    r = sp.Rational(radius)
    S = d["spin_two_form_factor_slope_upper"]
    assert d["transfer_disc_radius"] == r
    assert d["form_factor_quadratic_remainder_upper_coefficient"] == S / (4 - r) > 0
    assert d["negative_t_vertex_limit_error_per_abs_t_upper"] == 2 * S / (
        10**800 * (4 - r)
    )


def test_validation_not_bypassed_by_equal_cached_float():
    calibration.point(10**400, 1)
    for r in (True, 1.0, sp.Float(1)):
        with pytest.raises((TypeError, ValueError)):
            calibration.point(10**400, r)


@pytest.mark.parametrize("M", (2, 10, 100))
def test_two_triangle_integrand_at_actual_negative_transfer_sign(M):
    d = triangles.data()
    subs = {d["M"]: M, d["g"]: 1, d["x"]: sp.Rational(1, 3), d["z"]: sp.Rational(2, 5)}
    f = d["once_Ward_subtracted_unit_square_integrand"]
    assert sp.factor(f.subs(subs).subs(d["t"], 0)) == 0
    assert f.subs(subs).subs(d["t"], -1) < 0
    assert sp.diff(f, d["t"]).subs(subs).subs(d["t"], 0) > 0


def test_both_stress_placements_needed_for_Ward_weight():
    d = triangles.data()
    x = d["x"]
    d1 = d["on_shell_two_point_denominator"]
    L = d["light_stress_triangle"].subs(d["t"], 0)
    H = d["heavy_stress_triangle"].subs(d["t"], 0)
    assert sp.factor(L + H - x * (1 - x) / d1) == 0
    assert sp.factor(L - x * (1 - x) / d1) != 0
    assert sp.factor(H - x * (1 - x) / d1) != 0


def test_independent_pair_parameter_slope_weights_sum():
    d = triangles.data()
    x, z = d["x"], d["z"]
    weights = x * x * (1 - x) ** 3 * z * (1 - z) + x**3 * (1 - x) ** 2 * z * (1 - z)
    assert sp.factor(sp.integrate(weights, (z, 0, 1)) - x * x * (1 - x) ** 2 / 6) == 0


def test_local_quartic_projection_and_mass_independence():
    d = stress.data()
    assert d["projected_local_quartic_stress_bubble"] == 0
    names = {str(s) for s in d["projected_stress_vertex"].free_symbols}
    assert "line_mass_squared" not in names


def test_single_identical_field_symmetry_factor_and_c2_convention():
    d = pole.data()
    J = triangles.data()["equal_mass_two_distinct_field_integral"]
    expr = d["single_identical_field_equal_mass_diagnostic_c2"]
    g = next(s for s in expr.free_symbols if str(s) == "equal_mass_cubic_squared")
    Pl = next(s for s in expr.free_symbols if str(s) == "positive_Planck_mass")
    assert sp.simplify(expr + g * J / (48 * sp.pi**2 * Pl**2)) == 0


def test_graviton_finite_piece_scales_as_inverse_Planck_squared():
    d = calibration.point(10**400)
    doubled = calibration.point(2 * 10**400)
    key = "negative_t_vertex_finite_correction_absolute_upper"
    assert d[key] == 4 * doubled[key] > 0


def test_actual_slope_vertex_and_transfer_limit_have_distinct_scales():
    d = calibration.point()
    assert 0 < d["spin_two_form_factor_slope_upper"] < sp.Rational(1, 10**405)
    assert (
        0
        < d["negative_t_vertex_finite_correction_absolute_upper"]
        < sp.Rational(1, 10**1205)
    )
    assert d["negative_t_vertex_limit_error_per_abs_t_upper"] > 0


def test_exact_counts_and_original_scope():
    assert len(audit.residuals()) == 34
    assert len(audit.gates()) == 28
    assert len(audit.controls()) == 10
    assert audit.rejected_inputs() == 30
