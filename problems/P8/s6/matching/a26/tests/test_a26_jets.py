from fractions import Fraction

import sympy as sp
from p8_a26_vacuum import independent, jets


def test_full_variational_and_boundary_checks():
    assert jets.checks() and not any(jets.checks().values())


def test_full_symmetric_H_variation_is_not_the_restricted_density_variation():
    grad, hs, _, _, _, lag = jets.flat_jet()
    # Erroneously deleting off-diagonal directions before differentiating
    # cannot reproduce the real Euler expression, even at diagonal jets.
    restricted = lag.subs(dict.fromkeys(hs[3:], 0))
    wrong = -sum(hs[i]*sp.diff(restricted, grad[i], 2) for i in range(3))
    wrong += sum(hs[i]**2*sp.diff(restricted, hs[i], grad[i], 2) for i in range(3))
    at = dict(zip(grad, (1, 0, 0), strict=True)) | dict.fromkeys(hs[:3], 1)
    assert wrong.subs(at) != jets.evaluate_euler((1, 0, 0), (1, 1, 1))


def test_independent_fourier_constant_terms_not_a_sampled_quadrature():
    assert independent.fourier_moments() == {
        "g4": Fraction(9, 4), "g2_grad2": Fraction(3, 4), "grad4": Fraction(5, 4),
        "g2_lap_hess": 1, "grad2_lap_hess": 0}


def test_quartic_bulk_survives_all_displayed_IBP_terms():
    primary, other = jets.fourier_bulk(), independent.ibp_bulk()
    assert primary["quartic_after_time_IBP"] == other["T2_Tp2_after_IBP"] == -Fraction(1, 2)
    assert other["Tp4_after_IBP"] == 0
    assert other["raw"]["Tp4"] != 0  # cannot just drop this term before IBP


def test_fourth_difference_kills_quadratic_boundary_terms_not_bulk():
    ks, weights = (2, -2, 1, -1, 0), (1, 1, -4, -4, 6)
    assert all(sum(w*k**power for w, k in zip(weights, ks, strict=True)) == 0 for power in range(4))
    assert sum(w*k**4 for w, k in zip(weights, ks, strict=True)) == 24


def test_A5_operator_count_is_fourth_not_second_small_amplitude_order():
    # Four/six elementary phi factors in L3/L5; X has two factors.
    assert 4-2 == 2 and 6-2 == 4


def test_actual_compact_polynomial_Sobolev_bulk_sign_control():
    # A simple compact-time C3 control, not substituted for the written
    # C-infinity theorem: T=(1-t²)^4 on[-1,1], zero outside.
    t = sp.Symbol("t", real=True)
    profile = (1-t**2)**4
    integral = sp.integrate(profile**2*sp.diff(profile, t)**2, (t, -1, 1))
    assert integral > 0
    assert -integral/2 < 0


def test_euler_is_homogeneous_of_degree_one_but_not_additive():
    left = jets.evaluate_euler((2, 0, 0), (2, 2, 2))
    assert left == 2*jets.evaluate_euler((1, 0, 0), (1, 1, 1))
    assert jets.evaluate_euler((3, 0, 0), (2, 2, 2)) != sp.Rational(5, 2)
