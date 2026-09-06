import pytest
import sympy as sp
from p8_m1_weyl_scalar import geometry, normalization
from p8_m1_weyl_scalar import reduction as r


@pytest.mark.parametrize("function", [geometry.scalar_contraction_checks, geometry.gauge_and_ADM_checks,
                                      r.auxiliary_checks, r.field_map_checks, r.reconstruction_checks,
                                      r.rank_one_checks, normalization.generic_checks, normalization.baseline_checks])
def test_exact_geometry_constraint_and_phase_identities(function):
    assert all(sp.simplify(value) == 0 for value in function().values())


def test_crossing_constraints_and_matter_mixing_are_regular():
    data = r.phase_data()
    bounce = {r.H: 0, r.THETA: 0, r.DELTA: sp.Rational(1, 2),
              r.LAMBDA: -sp.Rational(1, 2), r.L: sp.Rational(1, 10)}
    lensing = sp.factor(data["lensing0"].subs(bounce))
    assert lensing == -(-20*r.Q*r.V+80*r.J*r.V+r.PS)/(40*r.J)
    assert sp.diff(lensing, r.PS) == -1/(40*r.J)
    assert sp.diff(lensing, r.P) == 0
    for key in ("n0", "b0", "h0", "lensing0"):
        value = data[key].subs(bounce)
        assert not value.has(sp.zoo, sp.nan)
        assert not sp.denom(sp.cancel(value)).has(r.THETA, r.LAMBDA)


def test_full_map_not_branch_map_off_shell():
    controls = r.controls()
    assert controls["omitting_E_p_field_map_term"] == r.EP**2/3
    assert all(sp.simplify(value) != 0 for value in controls.values())


def test_physical_maps_keep_matter_boundary_and_corrected_shift_distinct():
    rec = r.reconstruction()
    assert rec["pre_mixed_metric_momentum"] == r.P+3*r.L*r.S
    assert rec["matter_density_first_order_map"] == 3*r.L*rec["v1"]
    assert rec["b1"] == 4*r.LDOT
    assert r.cd(rec["n1"]+sp.Rational(8, 3)*r.DELTA*r.Q**2*r.phase_data()["lensing0"]/r.J) == 0
    assert r.old_flow_derivative(r.L) == -3*r.H*r.L
    assert r.old_flow_derivative(r.Q) == -2*r.H*r.Q


def test_canonical_generator_retains_both_channels_and_connection():
    w = sp.Matrix([[11, 2], [2, 13]])
    omega = sp.Matrix([[0, 3], [-3, 0]])
    eps = sp.Symbol("epsilon")
    data = normalization.canonical_generator(w, omega, [2, 3], [5, 7], 17, eps)
    assert data["old"][:2, :2] == -omega
    assert data["old"][2:, 2:] == -omega
    assert data["first_order_correction"].rank() == 1
    assert data["first_order_correction"][0, 3] != 0
    assert data["first_order_correction"][3, 0] != 0
    j = sp.BlockMatrix([[sp.zeros(2), sp.eye(2)], [-sp.eye(2), sp.zeros(2)]]).as_explicit()
    assert data["representative"].T*j+j*data["representative"] == sp.zeros(4)


def test_scalar_and_normalization_omission_controls():
    controls = geometry.controls()
    assert controls.pop("conformal_scalar_metric_has_no_Weyl") == 0
    assert all(sp.simplify(value) != 0 for value in controls.values())
    assert all(sp.simplify(value) != 0 for value in normalization.controls().values())


@pytest.mark.parametrize("chart", ["bad", "principal", "both"])
def test_unknown_phase_chart_fails(chart):
    with pytest.raises(ValueError):
        r.chart_data(chart)
