"""Fresh real and complex enclosures; the old coordinate box is not an old solution."""
from functools import cache

import sympy as sp
from p8_proca_nearby_cones import bounds as box_arithmetic
from p8_proca_nearby_packets.domains import complex_annulus

from . import model

T=sp.Rational(1,10**7)
N0=sp.Rational(1000001,1000000)
RN=sp.Rational(4,10**7)
RZ=sp.Rational(3,10**5)


def pair_box(numerator,denominator):
    if numerator.ring!=model.FIELD.ring or denominator.ring!=model.FIELD.ring:
        raise ValueError("Require the specified exact physical-coordinate polynomial ring")
    pn,pd=(box_arithmetic.polynomial_box(p) for p in (numerator,denominator))
    if pd["lower"]*pd["upper"]<=0:
        raise ValueError("The whole coordinate-box denominator reaches zero")
    endpoints=[a/b for a in (pn["lower"],pn["upper"]) for b in (pd["lower"],pd["upper"])]
    return {"lower":min(endpoints),"upper":max(endpoints),"numerator":pn,"denominator":pd}


def rational_box(value):
    if value.field!=model.FIELD:
        raise ValueError("Require the specified exact physical-coordinate rational field")
    return pair_box(value.numer,value.denom)


def interval(row):
    return row["lower"],row["upper"]


def multiply(a,b):
    endpoints=[x*y for x in a for y in b]
    return min(endpoints),max(endpoints)


def center(value):
    # Keep one generator during native FracElement evaluation; evaluate the
    # last rational expression explicitly, avoiding its all-generator API edge.
    return sp.factor(value.evaluate([(model.uf,0),(model.zf,0)]).as_expr().subs(model.N,N0))


@cache
def enclosures():
    d=model.system()
    names=("clock_physical_speed_squared","clock_physical_subluminal_margin",
        "on_constraint_rescaled_lapse_Hessian","constraint_matter_square","M",
        "lapse_flow","rescaled_trace_flow","hat_Hubble","conformal_log_flow")
    return {name:rational_box(d[name]) for name in names}


@cache
def physical_hubble():
    d=model.system()
    return (d["hat_Hubble"]+d["conformal_log_flow"])/model.Nf


@cache
def physical_hubble_partials():
    H=physical_hubble()
    p,q=H.numer,H.denom
    # Exact uncancelled quotient rule. No scientific GCD patch, and no need
    # to form a very large combined acceleration fraction before enclosure.
    return {name:pair_box(p.diff(j)*q-p*q.diff(j),q*q)
        for j,name in enumerate(("u","N","z"))}


@cache
def acceleration():
    partials=physical_hubble_partials()
    e=enclosures()
    terms=(interval(partials["u"]),
        multiply(interval(partials["N"]),interval(e["lapse_flow"])),
        multiply(interval(partials["z"]),interval(e["rescaled_trace_flow"])))
    total=tuple(sum(term[j] for term in terms) for j in (0,1))
    answer=multiply(total,(1/(N0+RN),1/(N0-RN)))
    return {"exact_chain_rule":"dH_physical/dtau=(H_u+H_N*N_dot+H_z*z_dot)/N",
        "three_coordinate_time_chain_term_enclosures":terms,
        "physical_proper_time_acceleration_lower":answer[0],
        "physical_proper_time_acceleration_upper":answer[1],
        "all_real_points_in_full_coordinate_box_not_only_center_samples":True}


@cache
def checks():
    u,N,z=model.FIELD.ring.gens
    toy=3+2*u+5*(N-box_arithmetic._qq(N0))**2+7*u*z
    row=pair_box(toy,model.FIELD.ring.one)
    old=model.old.system()
    new=model.system()
    new_w=center(new["constraint_matter_square"])
    old_w=center(old["constraint_matter_square"])
    return {"explicit_time_box_matches_arithmetic_only_not_old_solution":T-box_arithmetic.TIME_WINDOW,
        "explicit_new_lapse_center_matches_arithmetic_center":N0-box_arithmetic.LAPSE_CENTER,
        "explicit_lapse_box_radius":RN-box_arithmetic.LAPSE_RADIUS,
        "explicit_trace_box_radius":RZ-box_arithmetic.TRACE_RADIUS,
        "uncancelled_pair_box_keeps_exact_center":row["numerator"]["center"]-3,
        "uncancelled_pair_box_keeps_every_deviation_monomial":
            row["numerator"]["deviation_upper"]-(2*T+5*RN**2+7*T*RZ),
        "new_physical_Hubble_zero_at_exact_central_datum":center(physical_hubble()),
        "new_lapse_velocity_zero_at_exact_central_datum":center(new["lapse_flow"]),
        "using_old_matter_charge_leaves_exact_nonzero_new_constraint_residual":sp.factor(
            center(new["C"]-old["C"])-center(new["M"])*(old_w-new_w))}


@cache
def gates():
    e=enclosures()
    ann={name:complex_annulus(value) for name,value in e.items()}
    acc=acceleration()
    d=model.system()
    old=model.old.system()
    rows={"new_clock_squared_speed_strictly_above_493_over_500":
            e["clock_physical_speed_squared"]["lower"]>sp.Rational(493,500),
        "new_clock_squared_speed_strictly_below_247_over_250":
            e["clock_physical_speed_squared"]["upper"]<sp.Rational(247,250),
        "new_clock_subluminal_squared_margin_above_three_over_250":
            e["clock_physical_subluminal_margin"]["lower"]>sp.Rational(3,250),
        "new_clock_subluminal_squared_margin_below_seven_over_500":
            e["clock_physical_subluminal_margin"]["upper"]<sp.Rational(7,500),
        "new_real_lapse_pivot_negative_between_minus_four_and_minus_two":
            -4<e["on_constraint_rescaled_lapse_Hessian"]["lower"]
            and e["on_constraint_rescaled_lapse_Hessian"]["upper"]<-2,
        "new_complex_lapse_pivot_nonzero":ann["on_constraint_rescaled_lapse_Hessian"]["lower"]>2,
        "new_matter_square_positive_and_complex_nonzero":
            e["constraint_matter_square"]["lower"]>sp.Rational(9,1000)
            and ann["constraint_matter_square"]["lower"]>sp.Rational(9,1000)
            and ann["constraint_matter_square"]["upper"]<sp.Rational(11,1000),
        "unchanged_curvature_lapse_derivative_has_nonzero_complex_modulus":
            ann["M"]["lower"]>sp.Rational(24,100) and ann["M"]["upper"]<sp.Rational(26,100),
        "new_complex_lapse_velocity_below_one_e_minus_four":ann["lapse_flow"]["upper"]<sp.Rational(1,10000),
        "new_complex_trace_velocity_below_nine":ann["rescaled_trace_flow"]["upper"]<9,
        "new_complex_hat_Hubble_below_one_over_60000":ann["hat_Hubble"]["upper"]<sp.Rational(1,60000),
        "new_physical_bounce_acceleration_above_39997_over_10000_on_whole_box":
            acc["physical_proper_time_acceleration_lower"]>sp.Rational(39997,10000),
        "new_physical_bounce_acceleration_below_4001_over_1000_on_whole_box":
            acc["physical_proper_time_acceleration_upper"]<sp.Rational(4001,1000),
        "new_central_charge_is_not_old_central_charge":
            center(d["constraint_matter_square"])!=center(old["constraint_matter_square"]),
        "old_charge_does_not_solve_new_constraint_at_the_same_central_lapse":
            center(d["C"]-old["C"])!=0}
    return {name:bool(value) for name,value in rows.items()}
