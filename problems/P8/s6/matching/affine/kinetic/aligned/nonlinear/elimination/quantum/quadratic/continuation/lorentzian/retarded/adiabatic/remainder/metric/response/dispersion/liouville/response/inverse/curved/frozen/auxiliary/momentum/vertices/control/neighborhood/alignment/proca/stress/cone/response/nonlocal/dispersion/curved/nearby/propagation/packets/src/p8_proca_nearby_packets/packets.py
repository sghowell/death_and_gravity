"""Explicit three-dimensional real packets, diffraction and normalized tails."""
from functools import cache

import sympy as sp
from p8_proca_nearby_cones import bounds as cone_bounds

from . import domains, normal_form, observable


def carrier_exponent(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Integer)):
        raise TypeError("Require an integer decimal carrier exponent")
    value=int(value)
    if not 72<=value<=120 or value%3:
        raise ValueError("Use a carrier exponent divisible by three between 72 and120")
    return value


def packet(value=72):
    return _packet(carrier_exponent(value))


@cache
def _packet(exponent):
    carrier=sp.Integer(10)**exponent
    bandwidth=sp.Integer(10)**(exponent//3)
    ray=cone_bounds.characteristic_separation(domains.INNER_HALF_WIDTH)
    separation=ray["classical_comoving_characteristic_separation_lower"]
    radius=separation/8
    tails=6/(bandwidth*radius)
    error=observable.relative_error(carrier)
    diffraction=4*domains.T*bandwidth**2/carrier
    total=error["declared_relative_observable_error_upper"]+diffraction
    return {"carrier_decimal_exponent":exponent,"classical_carrier":carrier,"envelope_bandwidth":bandwidth,
            "comoving_frequency_band":(carrier/2,2*carrier),
            "initial_and_final_clock_times":(-domains.INNER_HALF_WIDTH,domains.INNER_HALF_WIDTH),
            "one_dimensional_normalized_envelope":"sqrt(sigma/pi) sin(sigma x)/(sigma x)",
            "real_three_dimensional_packet":"sqrt(2) cos(carrier*x1) product_i f_sigma(x_i)",
            "Fourier_support":"Two cubes centered at +/-carrier*e1 with each coordinate half-width sigma",
            "positive_carrier_outgoing_clock_mode_index":2,"positive_carrier_outgoing_matter_mode_index":3,
            "negative_carrier_conjugate_clock_mode_index":0,"negative_carrier_conjugate_matter_mode_index":1,
            "actual_relative_observable_error_upper":error["actual_relative_observable_error_upper"],
            "declared_relative_observable_error_upper":error["declared_relative_observable_error_upper"],
            "diffraction_relative_error_upper":diffraction,"total_relative_field_error_upper":total,
            "classical_ray_separation_lower":separation,"comparison_cube_half_width":radius,
            "leading_real_packet_outside_cube_L2_mass_upper":tails,
            "exact_packet_normalized_outside_cube_L2_mass_upper":8*(tails+total*total),
            "clock_cube_margin_beyond_matter_future_of_initial_cube":separation-2*radius,
            "same_leading_relational_field_initial_envelope_with_included_near_identity_correction":True,
            "both_exact_solutions_have_finite_band_and_all_finite_Sobolev_norms":True,
            "compact_spatial_support_or_compact_local_preparation":False,
            "interacting_EFT_band_or_quantum_resolution":False,
            "UV_exclusion_or_original_P8_closure":False}


@cache
def checks():
    sigma,x,xi,radius=sp.symbols("bandwidth spatial_coordinate Fourier_coordinate positive_radius",positive=True)
    integral=sp.integrate(sp.cos(xi*x),(xi,-sigma,sigma))/sp.sqrt(4*sp.pi*sigma)
    envelope=sp.sqrt(sigma/sp.pi)*sp.sin(sigma*x)/(sigma*x)
    norm=sp.integrate(1/(2*sigma),(xi,-sigma,sigma))
    tail=2*sp.integrate(x**-2,(x,radius,sp.oo))/(sp.pi*sigma)
    sample=packet()
    k1,kt2=sp.symbols("positive_longitudinal_momentum transverse_square",positive=True)
    return {
        "flat_spectral_envelope_inverse_Fourier_integral":sp.simplify(integral-envelope),
        "one_dimensional_spectral_envelope_unit_norm":norm-1,
        "three_dimensional_product_envelope_unit_norm":norm**3-1,
        "full_two_sided_sinc_tail_integral":sp.simplify(tail-2/(sp.pi*sigma*radius)),
        "spatial_dispersion_difference_exact_rationalization":
            sp.simplify(sp.sqrt(k1*k1+kt2)-k1-kt2/(sp.sqrt(k1*k1+kt2)+k1)),
        "bandwidth_is_exact_integer_cube_root_carrier":sample["envelope_bandwidth"]**3-sample["classical_carrier"],
        "inner_interval_keeps_half_parent_separation":sample["classical_ray_separation_lower"]-sp.Rational(99,4*10**14),
        "cube_radius_uses_only_certified_lower_separation":sample["comparison_cube_half_width"]-sp.Rational(99,32*10**14),
        "normalized_exact_packet_tail_keeps_field_error_and_normalization":
            sample["exact_packet_normalized_outside_cube_L2_mass_upper"]-8*(
                sample["leading_real_packet_outside_cube_L2_mass_upper"]+sample["total_relative_field_error_upper"]**2)}


@cache
def gates():
    d=packet()
    upper=packet(120)
    return {name:bool(value) for name,value in {
        "complete_frequency_band_inside_classical_normalform_domain":d["comoving_frequency_band"][0]>normal_form.MIN_MOMENTUM,
        "both_spectral_cubes_inside_declared_radial_band":3*d["envelope_bandwidth"]<d["classical_carrier"]/2,
        "opposite_real_carrier_spectral_cubes_disjoint":d["classical_carrier"]>2*d["envelope_bandwidth"],
        "actual_modal_reconstruction_below_declared_field_error":d["actual_relative_observable_error_upper"]<d["declared_relative_observable_error_upper"],
        "finite_field_error_below_one_e_minus_twenty_seven":d["total_relative_field_error_upper"]<sp.Rational(1,10**27),
        "field_norm_denominator_separated_from_zero":d["total_relative_field_error_upper"]<sp.Rational(1,2),
        "leading_packet_tail_below_one_e_minus_nine":d["leading_real_packet_outside_cube_L2_mass_upper"]<sp.Rational(1,10**9),
        "normalized_exact_packet_concentration_above_one_minus_one_e_minus_eight":
            d["exact_packet_normalized_outside_cube_L2_mass_upper"]<sp.Rational(1,10**8),
        "clock_detector_cube_outside_matter_future_of_initial_cube":d["clock_cube_margin_beyond_matter_future_of_initial_cube"]>0,
        "larger_declared_carriers_improve_both_errors_and_tails":upper["total_relative_field_error_upper"]<d["total_relative_field_error_upper"]
            and upper["leading_real_packet_outside_cube_L2_mass_upper"]<d["leading_real_packet_outside_cube_L2_mass_upper"],
        "sinc_tail_rational_weakening_uses_pi_above_two":2/sp.pi<1,
        "bandlimited_tails_are_not_claimed_zero":d["compact_spatial_support_or_compact_local_preparation"] is False,
        "no_interacting_EFT_frequency_or_quantum_resolution_claim":d["interacting_EFT_band_or_quantum_resolution"] is False,
        "original_P8_not_closed_by_classical_packet":d["UV_exclusion_or_original_P8_closure"] is False}.items()}


def controls():
    bad=(True,False,sp.true,sp.false,72.0,sp.Float(72),"72",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("free"),sp.sqrt(2),[],0,-1,69,73,121,sp.Rational(145,2))
    rejected=0
    for value in bad:
        for call in (carrier_exponent,packet):
            try:
                call(value)
            except (TypeError,ValueError):
                rejected+=1
    if rejected!=2*len(bad):
        raise ValueError("An inadmissible packet carrier was accepted")
    return {"rejected_inputs":rejected}
