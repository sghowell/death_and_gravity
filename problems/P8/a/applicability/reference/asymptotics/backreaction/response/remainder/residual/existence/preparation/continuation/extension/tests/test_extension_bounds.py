"""Complete map arithmetic, pointwise barriers and unchanged preparation."""

import json

import sympy as sp
from p8a_continuation import verify as a13
from p8a_extension import bounds, independent, regularity
from p8a_preparation import verify as a11


def test_every_constant_has_an_independent_fraction_replay():
    replay = independent.replay(json.loads(a11.REPORT.read_text()), json.loads(a13.REPORT.read_text()))
    assert independent.serialize(bounds.calibration()) == replay


def test_same_source_is_not_scaled_to_the_new_time_interval():
    data = bounds.calibration()
    assert data["source_window"] == sp.Rational(1, 10**10)
    assert data["length"] == sp.Rational(1, 10**6)
    assert data["extension_factor"] == 10000
    assert data["source_off_from"] == data["source_window"]/2
    assert data["source_off_from"] != data["length"]/2
    assert data["source_free_span"] == sp.Rational(19999, 20000000000)


def test_complete_rhs_includes_all_six_pair_contributions():
    pairs = bounds.calibration()["pairs"]
    total = sum(pairs[name] for name in ["Einstein_remainder", "curvature", "source",
                                        "auxiliary_product", "local_Wick"])
    assert total+pairs["actual_mode_response"]["total"] == pairs["full_G"]
    assert pairs["Einstein_full_P"] > 3
    assert pairs["Einstein_remainder"] < sp.Rational(4, 10**6)
    assert all(pairs[name] > 0 for name in ["curvature", "source", "local_Wick"])


def test_full_contraction_selfmap_and_sharper_actual_distance():
    data = bounds.calibration()
    gate = data["gate"]
    assert data["pairs"]["full_G"] < sp.Rational(9, 10**6)
    assert gate["contraction"] < sp.Rational(3, 10**6)
    assert gate["center_image"] < sp.Rational(9, 10**8)
    assert gate["self_map"] < data["weighted_radius"]
    assert gate["fixed_point_weighted"] < sp.Rational(9, 10**8)
    assert gate["fixed_point_pointwise"] < sp.Rational(81, 10**8)
    assert all(value > 0 for value in data["strict_margins"].values())


def test_P_then_q_bounds_are_verified_before_their_use():
    data = bounds.calibration()
    assert data["auxiliary_P_derived"] < data["auxiliary_P_cap"]
    assert data["auxiliary_q_derived"] < data["auxiliary_q_cap"]
    assert data["auxiliary_P_cap"] == sp.Rational(1, 10**5)
    assert data["auxiliary_q_cap"] == sp.Rational(2, 10**8)


def test_weighted_ball_is_not_mistaken_for_a_pointwise_ball():
    data = bounds.calibration()
    assert data["sigma"]*data["length"] == 2
    assert data["geometry"]["pointwise_ball_X"] == 9*data["weighted_radius"]
    assert data["geometry"]["pointwise_ball_X"] > sp.Rational(1, 10**6)
    # The fixed point itself, rather than every trial in the new ball,
    # satisfies the old unweighted uniqueness radius globally.
    assert data["gate"]["fixed_point_pointwise"] < sp.Rational(1, 10**6)
    assert data["geometry"]["pointwise_ball_W"] == sp.Rational(9, 10**13)


def test_same_flat_start_highest_jet_bound_is_strict():
    data = regularity.calibration()
    assert data["highest_jet_contraction"] < sp.Rational(1, 10**4)
    assert data["same_zero_neighborhood_end"] == sp.Rational(1, 4*10**10)
    assert data["new_inverse_is_a_generic_C1_endomorphism"] is False
    assert all(value == 0 for value in regularity.identities().values())


def test_elementary_weight_majorant_has_an_exact_series_witness():
    partial = sum(sp.Rational(1, sp.factorial(n)) for n in range(4))
    tail = sp.Rational(1, 24)/(1-sp.Rational(1, 5))
    assert partial+tail == sp.Rational(87, 32)
    assert sp.Rational(87, 32)**2 < 9
