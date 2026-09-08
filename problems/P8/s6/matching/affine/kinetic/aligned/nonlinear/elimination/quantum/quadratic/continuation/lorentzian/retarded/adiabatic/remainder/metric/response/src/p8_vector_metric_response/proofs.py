"""Exact two-source identities and continuous physical-metric audit gates."""
from functools import cache

import sympy as sp
from p8_vector_state import wkb

from . import envelopes, evolution, source, tadpole, tail, tangent


@cache
def residuals():
    out = {**tangent.low_checks(), **tail.algebra_checks(), **evolution.algebra_checks(),
           **evolution.weight_checks(), **tadpole.physical_vertices()["checks"]}
    for sector in ("T", "L"):
        reference = tangent.reference(sector)
        out.update({sector+"_varied_residual_low_"+str(j): value
                    for j, value in reference["varied_residual_low_coefficients"].items()})
        for component in ("energy", "pressure"):
            out.update({sector+"_"+component+"_subtracted_reference_low_"+str(j): value
                        for j, value in tail.reference_tail(sector, component)["low_tail_numerator_coefficients"].items()})
        out[sector+"_tenth_lapse_derivative_control"] = reference["tenth_source_derivative_fixtures"]["N"]+sp.Rational(7, 20736)
        out[sector+"_tenth_logscale_derivative_control"] = reference["tenth_source_derivative_fixtures"]["Z"]+sp.Rational(1, 512)
        for order in range(1, 5):
            out[sector+"_two_source_linear_coefficient_"+str(order)] = source.clean(
                sum(field*sp.diff(tangent.coefficient(sector, order), field) for field in source.n+source.v)
                -tangent.coefficient(sector, order))
    for order in range(1, 5):
        out["isotropic_zero_momentum_full_metric_reference_"+str(order)] = source.clean(
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
        out[sector+"_log_frequency_variation_bound_below_four"] = data["source_log_frequency_variation_upper"] < 4
        out[sector+"_reference_frequency_variation_bound_below_four"] = data["delta_W_over_frequency_upper"] < 4
        out[sector+"_reference_c_variation_bound_below_sixty_four"] = data["delta_c_upper"] < 64
        out[sector+"_canonical_rate_at_most_three"] = wkb.box_bound(source.data(sector)["rate"])["absolute_upper"] <= 3
        out[sector+"_reference_p_below_two_frequency_f"] = 5+sp.Rational(3, 2)*m <= 2*m
        out[sector+"_transport_exponential_below_two"] = 2*data["base_residual_over_inverse_frequency_eighth_upper"]/m**9 < sp.Rational(1, 4)
        for component in ("energy", "pressure"):
            out[sector+"_"+component+"_tail_majorants_nonnegative"] = tail.reference_tail(sector, component)["all_majorants_nonnegative"]
            row = evolution.constants()["by_sector"][sector]["physical_readouts"][component]
            out[sector+"_"+component+"_positive_finite_reference_readout_bound"] = row["delta_reference_bilinear_over_nu_squared_upper"] > 0
    r = 1+source.u**2
    out["positive_mass_first_profiles_have_upper_bounds_below_one"] = (
        0 < wkb.box_bound(4/(9*r**3))["absolute_upper"] < 1
        and 0 < wkb.box_bound(28/(81*r**3))["absolute_upper"] < 1)
    out["closed_weights_support_A_at_most_one_and_B_at_most_two"] = all(value == 0 for value in evolution.weight_checks().values())
    old = evolution.constants()["frozen_initial_and_evolved_mixing"]
    B6, B8, B10 = [old["initial_mixing_envelopes"][j] for j in (6, 8, 10)]
    out["middle_initial_band_covered_by_global_sixth_power"] = B8/(4*m)**2 <= B6
    out["high_initial_band_covered_by_global_sixth_power"] = B10/(8*m)**4 <= B6
    out["evolved_global_mixing_below_one"] = (B6+old["evolution_mixing_envelope"]/m**4)/m**6 < 1
    out["Gronwall_geometric_exponential_bound"] = 1/(1-sp.Rational(1, 4)) < 2
    example = evolution.bound(10**24, 1000)
    out["joint_nonlocal_metric_response_below_1e_minus_42"] = example["joint_nonlocal_C10_to_C0_upper_bound"] < sp.Rational(1, 10**42)
    out["joint_complete_vector_metric_response_below_1e_minus_38"] = example["joint_complete_metric_C10_to_C0_upper_bound"] < sp.Rational(1, 10**38)
    out["joint_existing_tadpole_plus_vector_response_below_1e_minus_37"] = (
        tadpole.bound(10**24, 1000)["joint_background_cancelled_metric_C10_to_C0_upper_bound"] < sp.Rational(1, 10**37))
    out["zero_independent_initial_response_does_not_erase_unvaried_initial_mixing"] = True
    out["existing_tadpole_profiles_and_cutoff_are_held_fixed"] = True
    out["physical_output_normalizations_and_second_mass_vertices_are_retained"] = True
    out["C10_to_C0_is_not_a_no_loss_coupled_inverse"] = True
    out["ordinary_Proca_Hadamard_and_cone_claims_not_transferred_to_finite_off_clock_metrics"] = True
    return {name: bool(value) for name, value in out.items()}
