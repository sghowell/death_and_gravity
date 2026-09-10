"""Independent exact kernel, subtraction and input-domain regressions."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_vacuum_light_pole import audit, calibration, kernel, normalization


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_actual_Hessian_or_analytic_identity(name, value):
    if isinstance(value, sp.MatrixBase):
        assert value == sp.zeros(*value.shape), name
    else:
        assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_exact_bound_and_written_argument_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args",
    calibration.bad_cases(),
    ids=[r[0] for r in calibration.bad_cases()],
)
def test_reject_unsupported_radius(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_control_and_scope(name, value):
    assert bool(value), name


@pytest.mark.parametrize("radius", (1, 2, Fraction(1, 2), sp.Rational(3, 2)))
def test_exact_rational_disc_calibration(radius):
    d = calibration.point(radius)
    k = kernel.data()
    r = sp.Rational(radius)
    M = k["actual_heavy_mass_squared"]
    g = k["actual_cubic_squared"]
    assert d["complex_disc_radius_about_mass_one"] == r
    assert d["quadratic_self_energy_upper_coefficient"] == g / (
        864 * M * M * (1 - r / M)
    )
    assert (
        d["factored_inverse_propagator_error_upper"]
        == r * d["quadratic_self_energy_upper_coefficient"]
    )
    assert d["no_extra_one_loop_pole_proved_in_disc"]


def test_validation_not_bypassed_by_equal_boolean_or_float_cache_key():
    calibration.point(1)
    for value in (True, 1.0, sp.Float(1)):
        with pytest.raises((TypeError, ValueError)):
            calibration.point(value)


def test_parameter_analyticity_alone_does_not_prove_zero_exclusion():
    d = kernel.data()
    M = d["actual_heavy_mass_squared"]
    g = d["actual_cubic_squared"]
    near = calibration.point(M - g / 10000)
    assert near["complex_disc_radius_about_mass_one"] < M
    assert near["factored_inverse_propagator_error_upper"] > 1
    assert not near["no_extra_one_loop_pole_proved_in_disc"]


def test_midpoint_weight_is_strictly_positive():
    d = kernel.data()
    assert (
        sp.factor(
            d["normalized_parameter_weight"].subs(
                d["Feynman_parameter"], sp.Rational(1, 2)
            )
            - 1 / (2 * d["heavy_mass_squared"] + 1)
        )
        == 0
    )


def test_full_two_by_two_mixed_insertion_contains_off_diagonal_terms():
    d = normalization.data()
    assert d["mixed_insertion"].shape == (2, 2)
    assert d["mixed_insertion"][0, 1] != 0
    assert d["mixed_insertion"][0, 1] == d["mixed_insertion"][1, 0]


def test_subtracted_analytic_kernel_has_a_double_zero():
    d = kernel.data()
    s = d["Minkowski_invariant"]
    f = d["local_analytic_on_shell_kernel"]
    assert f.subs(s, 1) == 0
    assert sp.diff(f, s).subs(s, 1) == 0
    assert (
        sp.factor(sp.diff(f, s, 2).subs(s, 1) - d["normalized_parameter_weight"] ** 2)
        == 0
    )


def test_positive_curvature_integrand_and_not_mass_one_curvature():
    d = kernel.data()
    h = next(iter(d["positive_curvature_conversion_integrand"].free_symbols))
    f = d["positive_curvature_conversion_integrand"]
    assert sp.diff(f, h).subs(h, 1) == sp.Rational(1, 2)
    phi = sp.Symbol("constant_canonical_light_field", real=True)
    p = sp.Symbol("positive_once_subtracted_curvature_shift", positive=True)
    assert sp.diff(d["on_shell_potential_quadratic"], phi, 2) == 1 - p


def test_bubble_threshold_has_an_interior_parameter_zero():
    x, mu = sp.symbols("x mu", positive=True)
    denominator = x * mu**2 + 1 - x - x * (1 - x) * (mu + 1) ** 2
    assert sp.factor(denominator.subs(x, 1 / (mu + 1))) == 0


def test_unit_disc_bound_matches_general_api():
    d = kernel.data()
    q = calibration.point()
    assert (
        q["quadratic_self_energy_upper_coefficient"]
        == d["actual_unit_disc_quadratic_self_energy_upper"]
    )
    assert (
        q["positive_curvature_conversion_upper"]
        == d["actual_rational_curvature_shift_upper"]
    )
    assert (
        q["finite_kinetic_counterterm_upper"]
        == d["actual_rational_finite_kinetic_counterterm_upper"]
    )


def test_strict_kinetic_bound_uses_three_not_one_as_rational_ceiling():
    d = kernel.data()
    assert (
        sp.Rational(1, 10**208)
        < d["actual_rational_finite_kinetic_counterterm_upper"]
        < sp.Rational(3, 10**208)
    )


def test_exact_counts_and_scope():
    rows = audit.residuals()
    assert len(rows) == 21
    assert (
        sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in rows.values()
        )
        == 24
    )
    assert len(audit.gates()) == 26
    assert len(audit.controls()) == 10
    assert audit.rejected_inputs() == 16
