"""Independent finite checks, exact invalid-input controls and scope regressions."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_vacuum_gauge_yukawa_screen import (
    audit,
    calibration,
    cuts,
    flow,
    model,
    threshold,
)


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_named_exact_identity(name):
    value = audit.residuals()[name]
    if isinstance(value, sp.MatrixBase):
        assert all(x == 0 for x in value)
    else:
        assert value == 0


@pytest.mark.parametrize("name", list(audit.gates()))
def test_every_explicit_gate(name):
    assert audit.gates()[name] is True


@pytest.mark.parametrize(
    "name,call,args", cuts.bad_cases(), ids=[x[0] for x in cuts.bad_cases()]
)
def test_invalid_channel_input(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_both_flavor_masses_and_generator_count():
    d = model.data()
    assert (
        d["all_Dirac_flavors"]
        == d["active_Dirac_flavors"] + d["spectator_Dirac_flavors"]
        == 14
    )
    assert len(d["gauge_generators"]) == d["gauge_dimension"] == 8
    assert d["flavor_mass_matrix"].shape == (2, 2)
    assert d["perturbative_gauge_boson_mass_squared"] == [0] * 8


def test_positive_ray_with_independent_integer_bracket():
    d = flow.data()
    root = d["quartic_to_gauge_squared"]
    assert 1 < root < 2
    assert sp.expand(675 * root**2 + 5310 * root - 11552) == 0
    assert d["Yukawa_squared_to_gauge_squared"] == sp.Rational(19, 45)


@pytest.mark.parametrize("n", range(1, 9))
def test_log_field_series_with_independent_fraction(n):
    d = threshold.data()
    assert d["log_field_coefficients"][str(n)] == sp.Rational(Fraction(-2, n))


def test_operator_vertex_from_two_field_strengths():
    eta = sp.diag(1, -1, -1, -1)
    k = sp.Matrix([1, 0, 0, 1])
    l = sp.Matrix([1, 0, 0, -1])
    eps = sp.Matrix([0, 1, 0, 0])
    F1 = k * eps.T - eps * k.T
    F2 = l * eps.T - eps * l.T
    contraction = sum(
        F1[i, j] * (eta * F2 * eta)[i, j] for i in range(4) for j in range(4)
    )
    C = cuts.data()["coefficient"]
    assert (
        sp.expand(
            4 * C * contraction
            - cuts.data()["transverse_polarization_amplitudes"][0].subs(
                sp.Symbol("positive_gauge_energy", positive=True), 1
            )
        )
        == 0
    )


@pytest.mark.parametrize(
    "s", (Fraction(1, 2), 1, Fraction(3, 2), 2, Fraction(5, 2), 3, Fraction(7, 2))
)
def test_exact_overlap_values(s):
    d = cuts.point(s)
    assert d["s"] + d["u"] == 4
    assert d["normalized_forward_boundary_jump"] == 16 * sp.I * sp.pi * (
        sp.Rational(s) - 2
    )


def test_center_and_second_derivative_cancellation_do_not_erase_neighboring_cut():
    d = cuts.data()
    jump = d["normalized_forward_boundary_jump"]
    assert cuts.point(2)["normalized_forward_boundary_jump"] == 0
    assert cuts.point(Fraction(3, 2))["normalized_forward_boundary_jump"] != 0
    assert sp.diff(jump, d["s"], 2) == 0
    assert sp.diff(jump, d["s"]) != 0


def test_zero_operator_control():
    d = cuts.data()
    assert d["s_channel_Im_amplitude"].subs(d["coefficient"], 0) == 0
    assert d["forward_log_coefficient"].subs(d["coefficient"], 0) == 0
    assert threshold.data()["leading_Phi_squared_F_squared_coefficient"] != 0


def test_prospective_coefficient_is_bounded_but_not_new_model_b2():
    d = calibration.data()
    assert 0 < d["leading_threshold_absolute_rational_upper"] < sp.Rational(1, 10**810)
    assert "NOT_new_model" in next(k for k in d if k.startswith("old_canonical"))
    assert "NOT computed" in d["scope"]
    assert (
        d["old_canonical_amplitude_budget_NOT_new_model_amplitude"][
            "highest_loop_order"
        ]
        == 2
    )


def test_exact_audit_counts_and_nonclosure():
    assert len(audit.residuals()) == 132
    assert audit.scalar_entries() == 248
    assert len(audit.gates()) == 38
    assert len(audit.controls()) == 9
    assert audit.rejected_inputs() == 16
    assert audit.controls()["full_V_G_B_and_original_P8_not_closed"] is True
