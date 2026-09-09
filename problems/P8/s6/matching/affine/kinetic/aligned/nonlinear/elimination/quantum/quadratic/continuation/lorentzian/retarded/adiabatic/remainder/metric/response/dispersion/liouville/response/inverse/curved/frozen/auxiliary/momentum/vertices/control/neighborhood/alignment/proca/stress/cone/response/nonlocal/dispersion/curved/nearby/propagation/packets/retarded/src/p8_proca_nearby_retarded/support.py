"""Compact source and detector neighborhoods strictly outside the matter cone."""
from functools import cache

import sympy as sp
from p8_proca_nearby_cones import bounds as cone_bounds
from p8_proca_nearby_packets import domains


def ordered_times(start,end):
    for value in (start,end):
        if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
            raise TypeError("Retarded times must be exact real rational data")
    start,end=sp.Rational(start),sp.Rational(end)
    if not -domains.INNER_HALF_WIDTH<start<end<domains.INNER_HALF_WIDTH:
        raise ValueError("Require strictly ordered interior retarded times")
    return start,end


def neighborhood(start=-domains.T/4,end=domains.T/4):
    start,end=ordered_times(start,end)
    gap=(end-start)*sp.Rational(99,100)*cone_bounds.CLOCK_SPEED_EXCESS_LOWER
    radius=min((start+domains.INNER_HALF_WIDTH)/4,
               (domains.INNER_HALF_WIDTH-end)/4,(end-start)/8,gap/100)
    return {"source_time":start,"detector_time":end,
            "source_spatial_center":sp.zeros(3,1),
            "detector_spatial_center":"(integral_source^detector omega_clock du,0,0)",
            "clock_minus_matter_ray_radius_lower":gap,
            "common_time_and_Euclidean_space_neighborhood_radius":radius,
            "matter_ray_change_upper_from_time_neighborhoods":4*radius,
            "source_and_detector_spatial_distance_loss_upper":2*radius,
            "every_pair_matter_spacelike_margin_lower":gap-6*radius,
            "every_pair_positive_clock_time_margin_lower":end-start-2*radius,
            "source_lower_time_margin":start-radius+domains.INNER_HALF_WIDTH,
            "detector_upper_time_margin":domains.INNER_HALF_WIDTH-end-radius,
            "nonzero_retarded_pairing_for_some_real_smooth_compact_source_and_detector":True,
            "no_explicit_compact_probe_frequency_band_or_UV_matching_error":True}


@cache
def data():
    return {"default_neighborhood":neighborhood(),
            "other_ordered_neighborhoods":[neighborhood(-domains.T/8,domains.T/8),
                                         neighborhood(0,domains.T/4)],
            "front_non_cancellation":"positive clock delta shell plus a locally L2 remainder; matter shell is disjoint",
            "compact_probe_existence":"nonzero distribution on a product of source and detector neighborhoods; product mollifiers separate distributions",
            "retarded_zero_data_before_source_in_original_canonical_phase":True,
            "classical_local_relational_response_not_supported_in_physical_matter_cone":True,
            "unbounded_classical_frequency_front_is_not_finite_EFT_causality":True}


@cache
def checks():
    d=neighborhood()
    return {"interior_pair_retains_half_of_parent_packet_time_span":
                d["detector_time"]-d["source_time"]-domains.T/2,
            "actual_clock_front_has_positive_declared_matter_cone_gap":
                d["clock_minus_matter_ray_radius_lower"]-sp.Rational(99,8*10**14),
            "compact_neighborhood_keeps_time_and_space_margin_costs":
                d["every_pair_matter_spacelike_margin_lower"]-
                (d["clock_minus_matter_ray_radius_lower"]
                 -d["matter_ray_change_upper_from_time_neighborhoods"]
                 -d["source_and_detector_spatial_distance_loss_upper"]),
            "default_neighborhood_uses_one_hundredth_of_front_gap":
                d["common_time_and_Euclidean_space_neighborhood_radius"]
                -d["clock_minus_matter_ray_radius_lower"]/100}


@cache
def gates():
    cases=[data()["default_neighborhood"],*data()["other_ordered_neighborhoods"]]
    return {name:bool(value) for name,value in {
        "clock_speed_excess_is_actual_whole_interval_parent_lower":
            cone_bounds.CLOCK_SPEED_EXCESS_LOWER==sp.Rational(1,400000),
        "all_declared_source_detector_pairs_are_strictly_matter_spacelike":
            all(row["every_pair_matter_spacelike_margin_lower"]>0 for row in cases),
        "all_declared_pairs_have_positive_clock_time_order":
            all(row["every_pair_positive_clock_time_margin_lower"]>0 for row in cases),
        "source_and_detector_neighborhoods_stay_inside_actual_interval":
            all(row["source_lower_time_margin"]>0 and row["detector_upper_time_margin"]>0 for row in cases),
        "source_and_detector_neighborhoods_have_positive_radius":
            all(row["common_time_and_Euclidean_space_neighborhood_radius"]>0 for row in cases),
        "separate_matter_delta_front_cannot_cancel_the_clock_front":True,
        "compact_probe_existence_does_not_claim_band_limited_compact_data":True,
        "time_retarded_is_not_synonymous_with_support_in_physical_matter_cone":True,
        "nearby_classical_response_not_transferred_to_old_quantum_profile":True}.items()}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("free"),sp.sqrt(2),[])
    cases=[(value,domains.T/4) for value in bad]+[(-domains.T/4,value) for value in bad]
    cases.extend(((0,0),(domains.T/4,0),(-domains.T,0),(0,domains.T),
                  (-domains.INNER_HALF_WIDTH,0),(0,domains.INNER_HALF_WIDTH)))
    rejected=0
    for start,end in cases:
        for call in (ordered_times,neighborhood):
            try:
                call(start,end)
            except (TypeError,ValueError):
                rejected+=1
    if rejected!=2*len(cases):
        raise ValueError("An inadmissible local-source time pair was accepted")
    return {"rejected_inputs":rejected}
