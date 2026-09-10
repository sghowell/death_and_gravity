"""Independent off-shell identities, exact jet inputs and both norm domains."""

import pytest
import sympy as sp
from p8_offshell_vacuum import audit, band, calibration, jets, matching, norms
from p8_polynomial_vacuum import calibration as parameters


@pytest.mark.parametrize(
    "name,value", list(audit.residuals().items()), ids=list(audit.residuals())
)
def test_exact_pointwise_or_parameter_identity(name, value):
    assert value == 0, name


@pytest.mark.parametrize(
    "name,value", list(audit.gates().items()), ids=list(audit.gates())
)
def test_norm_bound_and_written_proof_gate(name, value):
    assert bool(value), name


@pytest.mark.parametrize(
    "name,call,args",
    calibration.bad_cases(),
    ids=[row[0] for row in calibration.bad_cases()],
)
def test_exact_input_boundary(name, call, args):
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize(
    "name,value",
    [(k, v) for k, v in audit.controls().items() if k != "rejected_inputs"],
    ids=[k for k in audit.controls() if k != "rejected_inputs"],
)
def test_scope_and_negative_control(name, value):
    assert bool(value), name


def test_complete_symmetric_4D_jet_space():
    assert len(jets.JETS) == 210
    assert jets.SIGNS == (1, -1, -1, -1)
    assert max(map(sum, jets.INDICES)) == 6
    assert jets.JETS[(1, 1, 1, 1)] in jets.BY_SYMBOL


def test_coordinate_derivatives_commute_on_an_inhomogeneous_polynomial():
    f = jets.PHI**2 * jets.JETS[(1, 0, 1, 0)]
    assert jets.derivative(jets.derivative(f, 0), 3) == jets.derivative(
        jets.derivative(f, 3), 0
    )


def test_input_validation_precedes_derivative_cache_lookup():
    jets.derivative(jets.PHI, 0)
    for axis in (False, 0.0, sp.Integer(0)):
        with pytest.raises(ValueError):
            jets.derivative(jets.PHI, axis)


def test_actual_field_map_and_current_derivative_orders():
    d = jets.data()
    order = lambda expr: max(
        sum(jets.BY_SYMBOL[v]) for v in expr.free_symbols.intersection(jets.BY_SYMBOL)
    )
    assert order(d["cubic_field_redefinition"]) == 4
    assert order(d["literal_centered_resolvent_truncation"]) == 6
    assert max(order(v) for v in d["full_boundary_current"]) == 5


def test_independent_constant_off_shell_field_map_calibration():
    d = jets.data()
    p = parameters.point()
    l, g, c = sp.symbols("quartic_lambda quartic_gamma redundant_c", real=True)
    replace = {v: 0 for v in jets.JETS.values()}
    replace[jets.PHI] = 1
    R = d["cubic_field_redefinition"].subs(replace, simultaneous=True)
    assert sp.expand(R - (-c / 6 + l - 2 * g / 3)) == 0
    assert (
        R.subs({l: p["lambda"], g: p["gamma"], c: p["cubic_squared"] / p["D"] ** 2}) < 0
    )


def test_exact_coefficient_group_majorants():
    n = norms.data()
    assert n["field_map_coefficient_groups"] == {
        "quartic_lambda": 13,
        "quartic_gamma": sp.Rational(434, 3),
        "redundant_c": sp.Rational(1, 6),
    }
    assert n["quartic_action_coefficient_groups"] == {
        "quartic_lambda": 81,
        "quartic_gamma": 729,
        "redundant_c": sp.Rational(13, 6),
    }


@pytest.mark.parametrize("B", (0, sp.Rational(1, 2), 1, 2))
def test_exact_general_jet_radius(B):
    d = calibration.point(38, B)
    n = norms.data()
    r = sp.Symbol("canonical_jet_radius", nonnegative=True)
    assert d["pointwise_field_redefinition_error"] == n[
        "general_radius_pointwise_error"
    ].subs(r, B)
    assert (
        d["six_derivative_jet_change_bound"]
        == n["all_first_six_derivatives_majorant"] * B**3
    )
    assert d["pointwise_field_redefinition_error"] >= 0


def test_large_spectral_window_and_common_window_kept_separate():
    wide = calibration.point()
    small = calibration.point(38)
    assert wide["radius"] == 10**144 > small["radius"] == 38
    assert (
        0
        < small["spectral_quartic_action_error_coefficient"]
        < wide["spectral_quartic_action_error_coefficient"]
        < sp.Rational(1, 10**400)
    )


def test_common_integrated_bound_uses_the_mapped_source_norm():
    d = calibration.point(38)
    b = band.data()
    n = norms.data()
    assert b["mapped_source_L2_squared_factor"] == (1 + n["actual_C_R"]) ** 4
    assert (
        b["combined_tree_action_error_coefficient"]
        == d["pointwise_field_redefinition_error"]
        + d["spectral_quartic_action_error_coefficient"] * (1 + n["actual_C_R"]) ** 4
    )
    assert b["combined_tree_action_error_coefficient"] < sp.Rational(1, 10**800)


def test_actual_quartic_coefficient_from_affine_source():
    d = matching.actual_target()
    gamma = matching.data()["gamma"]
    assert d["actual_canonical_A3_quartic"] == -2 * gamma


def test_exact_counts_and_original_scope():
    assert len(audit.residuals()) == 24
    assert len(audit.gates()) == 34
    assert len(audit.controls()) == 13
    assert audit.rejected_inputs() == 35
