"""Algebraic and continuous audit gates for the finite prepared response."""
from functools import cache

import sympy as sp
from p8_vector_state import wkb

from . import envelopes, evolution, source, tail, tangent


@cache
def residuals():
    out = {**tangent.low_checks(), **tail.algebra_checks(), **evolution.algebra_checks()}
    for sector in ("T", "L"):
        reference = tangent.reference(sector)
        out.update({sector+"_varied_residual_low_"+str(j): value
                    for j, value in reference["varied_residual_low_coefficients"].items()})
        out.update({sector+"_subtracted_reference_low_"+str(j): value
                    for j, value in tail.reference_tail(sector)["low_tail_numerator_coefficients"].items()})
        out[sector+"_tenth_derivative_nonzero_control"] = reference["tenth_source_derivative_fixture"]+sp.Rational(7, 20736)
        for order in range(1, 5):
            out[sector+"_source_linear_coefficient_"+str(order)] = source.clean(
                sum(field*sp.diff(tangent.coefficient(sector, order), field) for field in source.n)
                -tangent.coefficient(sector, order))
    for order in range(1, 5):
        out["isotropic_zero_momentum_varied_reference_"+str(order)] = source.clean(
            (tangent.varied_P("T", order)-tangent.varied_P("L", order)).subs(source.z, 0))
    return out


@cache
def checks():
    m = wkb.MASS_TIME_MIN
    out = {}
    for sector in ("T", "L"):
        data = envelopes.reference(sector)
        out[sector+"_all_source_box_reconstructions"] = all(value == 0 for value in data["box_reconstructions"])
        out[sector+"_positive_base_reference_and_transport_bounds"] = all(data["base_reference_checks"].values())
        out[sector+"_source_A_variation_at_most_two"] = data["source_A_variation_upper"] <= 2
        out[sector+"_source_B_variation_at_most_two"] = data["source_B_variation_upper"] <= 2
        out[sector+"_source_log_frequency_variation_at_most_one"] = data["source_log_frequency_variation_upper"] <= 1
        out[sector+"_frequency_variation_upper_below_one_half"] = data["delta_W_over_frequency_upper"] < sp.Rational(1, 2)
        out[sector+"_reference_c_variation_upper_below_seven"] = data["delta_c_upper"] < 7
        for key in ("A", "B"):
            out[sector+"_unvaried_weight_"+key+"_at_most_one"] = wkb.box_bound(source.data(sector)[key])["absolute_upper"] <= 1
        out[sector+"_canonical_rate_at_most_three"] = wkb.box_bound(source.data(sector)["rate"])["absolute_upper"] <= 3
        out[sector+"_reference_p_below_two_frequency_f"] = 5+sp.Rational(3, 2)*m <= 2*m
        out[sector+"_transport_exponential_below_two"] = 2*data["base_residual_over_inverse_frequency_eighth_upper"]/m**9 < sp.Rational(1, 4)
        out[sector+"_reference_tail_majorants_nonnegative"] = tail.reference_tail(sector)["all_majorants_nonnegative"]
    old = evolution.constants()["frozen_initial_and_evolved_mixing"]
    B6, B8, B10 = [old["initial_mixing_envelopes"][j] for j in (6, 8, 10)]
    out["middle_initial_band_covered_by_global_sixth_power"] = B8/(4*m)**2 <= B6
    out["high_initial_band_covered_by_global_sixth_power"] = B10/(8*m)**4 <= B6
    out["evolved_global_mixing_below_one"] = (B6+old["evolution_mixing_envelope"]/m**4)/m**6 < 1
    out["Gronwall_geometric_exponential_bound"] = 1/(1-sp.Rational(1, 4)) < 2
    example = evolution.bound(10**24, 1000)
    out["nonlocal_prepared_response_below_1e_minus_43"] = example["nonlocal_subtracted_C10_to_C0_upper_bound"] < sp.Rational(1, 10**43)
    out["complete_prepared_mass_response_below_1e_minus_39"] = example["complete_mass_response_C10_to_C0_upper_bound"] < sp.Rational(1, 10**39)
    return {name: bool(value) for name, value in out.items()}
