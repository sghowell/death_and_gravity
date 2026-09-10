"""Full-vertex normalization, subtracted kernels and exact loop margins."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_vacuum_forward_loop import (
    audit,
    calibration,
    domain,
    leading,
    normalization,
    subtraction,
)


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_full_vertex_kernel_or_counterterm_identity(name, value):
    if isinstance(value, sp.MatrixBase):
        assert value == sp.zeros(*value.shape), name
    else:
        assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_exact_majorant_and_written_proof_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args",
    calibration.bad_cases(),
    ids=[r[0] for r in calibration.bad_cases()],
)
def test_reject_unsupported_relative_tolerance(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_negative_control_and_scope(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "eps", (Fraction(1, 10**6), sp.Rational(1, 10**5), sp.Rational(1, 2))
)
def test_supported_exact_tolerance(eps):
    d = calibration.point(eps)
    assert d["relative_tolerance"] == sp.Rational(eps)
    assert d["requested_tolerance_proved"]


def test_validation_not_bypassed_by_boolean_or_float_alias():
    calibration.point(sp.Rational(1, 2))
    for value in (True, 0.5, sp.Float(0.5)):
        with pytest.raises((ValueError, TypeError)):
            calibration.point(value)


def test_stricter_tolerance_not_inferred_from_positive_tree_margin():
    d = calibration.point(sp.Rational(1, 10**8))
    assert not d["requested_tolerance_proved"]
    assert d["tree_plus_one_loop_b2_lower"] > 0


def test_exact_actual_error_retains_both_distinct_parts():
    d = calibration.point()
    assert (
        d["total_one_loop_b2_error_upper"]
        == d["angular_b2_error_upper"] + d["angular_zero_b2_error_upper"]
    )
    assert 0 < d["angular_zero_b2_error_upper"] < d["angular_b2_error_upper"]
    assert d["relative_one_loop_error_upper"] < sp.Rational(1, 10**6)


def test_full_quartic_tensor_and_half_trace_factor():
    d = normalization.data()
    assert d["actual_quartic_vertex_tensor"].shape == (4, 4)
    assert d["actual_quartic_vertex_tensor"][0, 3] != 0
    assert d["three_channel_symmetry_factor"] == sp.Rational(1, 2)
    assert d["four_dimensional_radial_channel_prefactor"] == 1 / (32 * sp.pi**2)


def test_high_radial_rewrite_is_nontrivial_and_all_terms_retained():
    d = leading.data()
    assert len(d["exact_UV_cancellation_terms"]) == 3
    assert all(v != 0 for v in d["exact_UV_cancellation_terms"])
    assert d["high_region_constant"] == sp.Rational(6639, 16)


def test_heavy_external_shift_geometric_bound_has_correct_maximum():
    y, M = sp.symbols("y M", positive=True)
    f = y / (M + y) ** 2
    assert sp.diff(f, y).subs(y, M) == 0
    assert f.subs(y, M) == 1 / (4 * M)
    assert sp.factor(1 / (4 * M) - f - (y - M) ** 2 / (4 * M * (M + y) ** 2)) == 0


@pytest.mark.parametrize("y", (0, 1, 32, 64))
def test_independent_complex_radial_vertex_and_derivative_bounds(y):
    d = leading.data()
    subs = {
        d["s"]: 2 + sp.I / 2,
        d["a"]: sp.Rational(1, 4),
        d["radial_y"]: y,
        d["M"]: 32,
        d["small_c"]: 1,
        sp.Symbol("quartic_L", positive=True): 100,
    }
    v = d["radial_V0"]
    norm2 = lambda z: sp.cancel(sp.expand_complex(z * sp.conjugate(z)))
    assert norm2(v.subs(subs)) <= (20 * sp.Rational(32) * (y + 1) / (32 + y)) ** 2
    assert norm2(sp.diff(v, d["s"]).subs(subs)) <= 36
    assert norm2(sp.diff(v, d["s"], 2).subs(subs)) <= sp.Rational(18, 32) ** 2


def test_full_parameter_Delta_remains_positive_at_boundary_parameter_faces():
    d = domain.data()
    a0, a1, a2, a3 = sp.symbols(
        "alpha_light0 alpha_light1 alpha_heavy2 alpha_heavy3", nonnegative=True
    )
    M = sp.Symbol("heavy_mass_squared", positive=True)
    for H in (0, sp.Rational(1, 2), 1):
        value = d["full_box_parameter_Delta"].subs(
            {
                a0: (1 - H) / 2,
                a1: (1 - H) / 2,
                a2: H / 2,
                a3: H / 2,
                d["forward_s"]: 3,
                M: 32,
            }
        )
        assert value >= sp.Rational(1, 4) + (32 - 1) * H


def test_local_counterterm_contains_contact_single_and_double_heavy_poles():
    d = subtraction.data()
    h = sp.symbols("h_s h_t h_u", real=True)
    p = sp.Poly(d["full_UV_counterterm_amplitude"], *h)
    assert p.total_degree() == 2
    assert p.coeff_monomial(h[0]) != 0
    assert p.coeff_monomial(h[0] ** 2) != 0
    assert p.coeff_monomial(1) != 0


def test_exact_counts_and_original_scope():
    rows = audit.residuals()
    assert len(rows) == 54
    assert (
        sum(
            v.rows * v.cols if isinstance(v, sp.MatrixBase) else 1
            for v in rows.values()
        )
        == 69
    )
    assert len(audit.gates()) == 44
    assert len(audit.controls()) == 11
    assert audit.rejected_inputs() == 16
