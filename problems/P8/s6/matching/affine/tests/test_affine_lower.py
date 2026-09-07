"""Check the fixed lower-order ODE, boundary and source-transcription control."""
import pytest
import sympy as sp
from p8_affine import lower


@pytest.mark.parametrize("name", ["integrating_factor_derivative", "ODE_from_general_p",
                                  "forcing_from_general_p", "cubic_solves_Q1", "ODE_solves_Q2",
                                  "scalar_after_boundary", "clock_forcing_zero", "clock_p_phi_zero"])
def test_matching_algebra(name):
    assert lower.identities()[name] == 0


@pytest.mark.parametrize("name", ["ten_Hessian_jet_vector_divergence", "integration_by_parts_scalar_sign"])
def test_boundary_derived_with_all_jets(name):
    assert lower.divergence_identities()[name] == 0


def test_palatini_rejects_printed_formula_away_from_unit_X():
    result = lower.palatini_control()
    assert result["corrected_cross_residual"] == 0
    assert result["printed_cross_residual"] != 0
    assert result["printed_at_nonunit"] == -sp.Rational(135, 2)
    assert result["constant_X_normalization_can_hide_typo"] == 0


def test_primitive_stays_inside_the_clock_tube():
    equation = lower.ode()
    assert equation["basepoint"] == -1
    integral, = equation["q"].atoms(sp.Integral)
    assert integral.limits == ((equation["integration_variable"], -1, lower.d.x),)
    assert not equation["integrand"].has(sp.Float)


def test_zero_reference_data_for_every_time():
    equation = lower.ode()
    assert equation["q"].subs(lower.d.x, -1).doit() == 0
