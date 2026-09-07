"""Independent arithmetic and omission controls for the fixed-frame dictionary."""

from fractions import Fraction

import pytest
import sympy as sp
from p8_variable_reduction import dictionary


def test_all_dictionary_residuals_vanish():
    result = dictionary.checks()
    assert len(result) == 27
    assert set(result.values()) == {sp.S.Zero}


def test_actual_center_and_best_planck_fit_from_fraction_data():
    for c in (Fraction(2001, 1000), Fraction(5, 2), Fraction(4)):
        for x in (Fraction(91, 100), Fraction(1), Fraction(109, 100)):
            got = dictionary.center(c, x, 7, 3)
            b1 = Fraction(7*32, 9)/(c*(c-2))
            kap = Fraction(7*9, 16)*c*(c-2)
            assert got["beta1"] == b1
            assert got["kappa"] == kap
            assert got["M_star_squared"] == 35
            assert got["CD"]["GT"] == 35*x
            assert got["formal"]["GT"] == 35
            assert got["normalized_F2_defect"] == (x-1)/2
            assert got["I_defect"] == -1


def test_scalar_clock_weights_from_actual_hessian_transformation():
    j, jp, x, h2, box, vh = sp.symbols("J Jprime X H2 box vHv", real=True)
    # phi_;mu nu=J theta_;mu nu+Jprime theta_,mu theta_,nu.
    transformed_l1 = j*j*h2+2*j*jp*vh+jp*jp*x*x
    transformed_l2 = (j*box+jp*x)**2
    assert sp.expand(transformed_l1-transformed_l2
                     -j*j*(h2-box*box)-2*j*jp*(vh-x*box)) == 0
    data = dictionary.reparametrization()
    z = data["symbols"]
    assert data["A1_theta"] == z["J"]**2*z["A1"]
    assert data["A3_theta"] == z["J"]**4*z["A3"]
    assert sp.cancel(data["I_theta"]-data["I_parent"]) == 0
    # Xi alone has weight J², rather than being invariant without X/GT.
    assert sp.factor(data["Xi_theta"]-data["Xi_parent"]) != 0


def test_full_weighted_ibp_jet_fixture():
    data = dictionary.coefficients()
    z = data["symbols"]
    p = z["phi"]
    q = 2+p+3*p*p
    kap = 5+7*p+11*p*p+13*p**3
    replacements = {z["q"]: q, z["kappa"]: kap}
    evaluate = lambda key: sp.expand(data[key].subs(replacements).doit()).subs(p, 0)
    assert evaluate("F2_X") == -28
    assert evaluate("A1") == -56
    assert evaluate("A2") == 56
    assert evaluate("A3") == 0
    assert evaluate("K_X") == 474
    assert evaluate("F_X2_coefficient") == 940


def test_constant_kappa_conformal_control():
    data = dictionary.coefficients()
    p = data["symbols"]["phi"]
    replacements = {data["symbols"]["kappa"]: sp.Integer(7),
                    data["symbols"]["q"]: 1+p+p*p}
    for key in ("F2_X", "A1", "A2", "A3", "K_X", "F_X2_coefficient"):
        assert sp.simplify(data[key].subs(replacements).doit()) == 0
    assert data["curvature_squared_coefficient"] != 0


def test_missing_kappa_derivatives_is_not_the_dictionary():
    data = dictionary.coefficients()
    z = data["symbols"]
    fake_a1 = 4*z["kappa"]*z["q"]**2
    assert sp.factor(fake_a1-data["A1"]) != 0
    assert data["A1"] == 2*data["F2_X"]


def test_coefficient_error_floor_and_approximate_saturation():
    got = dictionary.remainder_floor(Fraction(1, 10), 5)
    assert got["weighted_C1_floor"] == Fraction(9, 2)
    # Saturation of the necessary triangle bound; it is not a constructed EFT.
    delta_a1, delta_f2x = 0, Fraction(-9, 4)
    assert abs(delta_a1)+2*abs(delta_f2x) == got["weighted_C1_floor"]
    inv = dictionary.invariant_remainder_floor(Fraction(1, 10), Fraction(1, 5), 5)
    assert inv["weighted_C1_floor"] == Fraction(18, 5)


def test_value_only_matching_misses_the_derivative_defect():
    point = dictionary.center()
    assert point["formal"]["F2"] == point["CD"]["F2"]
    assert point["formal"]["F2_X"]-point["CD"]["F2_X"] == Fraction(5, 2)
    assert point["formal"]["I"] != point["CD"]["I"]


def test_no_matter_or_hierarchy_promotion_flags():
    got = dictionary.calibration()
    assert got["curvature_square_retained"]
    assert got["order_six_can_contribute"]
    assert not got["controlled_heavy_EFT_claimed"]
    assert not got["UV_or_cutoff_claimed"]
    assert not got["center"]["matter_background_retuned"]


def test_inverse_X_prevents_naive_six_derivative_target_assignment():
    epsilon, x, l3 = sp.symbols("epsilon X L3", positive=True)
    # Formal derivative grading: dphi has degree1 and ddphi degree2.
    actual_scaling = epsilon**6*l3/(epsilon**2*x)
    assert sp.cancel(actual_scaling-epsilon**4*l3/x) == 0
    assert sp.cancel(actual_scaling-epsilon**6*l3/x) != 0


def test_inverse_metric_source_sign_by_direct_single_component_variation():
    eps, kap, q, velocity, rr = sp.symbols("eps kap Q velocity R00", positive=True)
    # Ricci only has R00=rr. delta g^00=5*kappa*rr/(3Q).
    # Its three spatial inverse entries change by kappa*rr/(3Q).
    g00 = 1+eps*5*kap*rr/(3*q)
    spatial = -1+eps*kap*rr/(3*q)
    volume = (-g00*spatial**3)**sp.Rational(-1, 2)
    matter = volume*g00*velocity**2/2
    variation = sp.diff(matter, eps).subs(eps, 0)
    assert sp.factor(variation-sp.Rational(2, 3)*kap*rr*velocity**2/q) == 0
    contact = dictionary.contacts()
    z = contact["symbols"]
    expression = contact["chi_curvature_contact"].subs(
        {z["kappa"]: kap, z["Q"]: q, z["Ric_chichi"]: rr*velocity**2,
         z["R"]: rr, z["Y"]: velocity**2})
    assert sp.factor(expression-variation) == 0
    assert sp.factor(expression+variation) != 0


def test_source_schur_preserves_massive_exchange_and_contact():
    data = dictionary.constant_r_schur()
    z = data["symbols"]
    g, f, nu, d = z["G"], z["F"], z["nu"], z["D"]
    matrix = sp.Matrix([[g*d+nu, -nu], [-nu, f*d+nu]])
    assert sp.cancel(matrix.inv()[0, 0]-data["response"]) == 0
    heavy_contact = sp.limit(data["response"]-1/((g+f)*d), d, 0)
    assert sp.factor(heavy_contact-f*f/(nu*(g+f)**2)) == 0
    assert sp.factor(data["source_tensor_contact_coefficient"]-heavy_contact/2) == 0
    assert sp.factor(data["pure_chi_Y2_coefficient"]
                     -sp.Rational(2, 3)*data["source_tensor_contact_coefficient"]) == 0


@pytest.mark.parametrize("bad", [True, False, 0.1, sp.Float("0.1"), sp.oo,
                                  -sp.oo, sp.nan, "1", sp.Symbol("x")])
def test_strict_exact_input_types(bad):
    with pytest.raises(TypeError):
        dictionary.center(x=bad)
    with pytest.raises(TypeError):
        dictionary.remainder_floor(bad)


@pytest.mark.parametrize("kwargs", [{"c": 2}, {"c": 1}, {"c": 5},
                                    {"x": Fraction(9, 10)}, {"x": Fraction(11, 10)},
                                    {"parent_planck_squared": 0}, {"time_scale": -1}])
def test_invalid_center_domains(kwargs):
    with pytest.raises(ValueError):
        dictionary.center(**kwargs)


@pytest.mark.parametrize("error", [-1, 1, 2])
def test_invalid_relative_errors(error):
    with pytest.raises(ValueError):
        dictionary.remainder_floor(error)
    with pytest.raises(ValueError):
        dictionary.invariant_remainder_floor(invariant_error=error)


def test_nonpositive_target_normalization_rejected():
    with pytest.raises(ValueError):
        dictionary.remainder_floor(target_planck_squared=0)
    with pytest.raises(ValueError):
        dictionary.invariant_remainder_floor(target_planck_squared=-1)
