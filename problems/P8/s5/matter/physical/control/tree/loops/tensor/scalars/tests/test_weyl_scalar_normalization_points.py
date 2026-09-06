import pytest
import sympy as sp
from p8_m1_control import model
from p8_m1_weyl_scalar import bounds as b
from p8_m1_weyl_scalar import normalization


@pytest.mark.parametrize("chart,x,scale", [("unitary", sp.Rational(3, 5), sp.Rational(3, 2)),
                                          ("gamma", sp.Integer(0), sp.Rational(2, 3))])
def test_exact_full_normalization_at_rational_time_and_nonunit_ell(chart, x, scale):
    q = sp.Integer(10)**22
    point = normalization.point(chart, x, q, scale)
    mapping = {model.x: x, model.z: 1/(scale**2*q), model.l: (1-x*x)**sp.Rational(11, 2)/10}
    alpha_bar = model.coefficients(chart)["alpha"].subs(mapping, simultaneous=True)
    dimensions = sp.diag(1/scale if chart == "gamma" else 1, 1)
    assert all(sp.simplify(value) == 0 for value in point["T_fixed"].T*point["T_fixed"]-dimensions*alpha_bar*dimensions)
    assert point["boundary_fixed"] == point["boundary_fixed"].T
    assert point["connection_fixed"] == -point["connection_fixed"].T
    ry, rp = point["lensing_Y_row"], point["lensing_P_row"]
    assert sp.simplify((rp.T*rp)[0]) < b.NORMALIZED_ROW_P**2
    assert sp.simplify((ry.T*ry)[0]) < (b.NORMALIZED_ROW_Y*q)**2
    dual_squared = sp.simplify((ry.T*point["potential_fixed"].inv()*ry+rp.T*rp)[0])
    assert dual_squared < b.DUAL_ENERGY_ROW_SQUARED*q


def test_exact_point_api_rejects_rounded_or_uncovered_input():
    with pytest.raises(ValueError):
        normalization.point("gamma", 0.0, 10**22)
    with pytest.raises(ValueError):
        normalization.point("unitary", 0, 10**22)
    with pytest.raises(ValueError):
        normalization.point("gamma", 0, 10**22, 3)
