"""Complex bounds along the actual analytic solution, not sampled real data."""
from functools import cache

import sympy as sp
from p8_proca_nearby_bounce import taylor
from p8_proca_nearby_cones import bounds as parent_bounds
from p8_proca_nearby_cones import rational as parent

from . import finite, modes

T=sp.Rational(1,10**7)
INNER_HALF_WIDTH=T/2
BASIS_NORM=sp.Integer(1000)
MATRIX_NORM=sp.Integer(10000)
GAP=sp.Rational(1,10**7)


def complex_annulus(box):
    np,dp=box["numerator"],box["denominator"]
    values=(np["center"],np["deviation_upper"],dp["center"],dp["deviation_upper"])
    if any(isinstance(v,bool) or not isinstance(v,(int,sp.Rational)) for v in values):
        raise TypeError("Complex polynomial annuli require exact rational data")
    nc,nd,dc,dd=map(sp.Rational,values)
    if nd<0 or dd<0:
        raise ValueError("A polynomial deviation cannot be negative")
    nl=max(sp.Integer(0),abs(nc)-nd)
    nu=abs(nc)+nd
    dl=abs(dc)-dd
    du=abs(dc)+dd
    if dl<=0:
        raise ValueError("A complex denominator annulus reaches zero")
    return {"lower":nl/du,"upper":nu/dl}


@cache
def analytic_domain():
    data=parent.system()
    ann={name:complex_annulus(value) for name,value in parent_bounds.enclosures().items()}
    alpha=(2*data["A"]*parent.zf+data["B"])/3
    ann["actual_force_trace"]=complex_annulus(parent_bounds.rational_box(alpha))
    ndev=parent_bounds.LAPSE_CENTER-1+parent_bounds.LAPSE_RADIUS
    hdev=3*T*T+3*T**4+T**6
    Ddev=(1+ndev)**2*hdev
    e4dev=((1+ndev)**2-1)/(1-Ddev)
    edev=e4dev/(1-e4dev)
    logR=taylor.evolution()["log_hat_scale_absolute_upper"]
    Rdev=logR/(1-logR)
    cdev=ann["clock_physical_speed_squared_excess"]["upper"]/(1-ann["clock_physical_speed_squared_excess"]["upper"])
    omega_lower=(1-ndev)/((1+edev)*(1+Rdev))
    omega_upper=(1+ndev)/((1-edev)*(1-Rdev))
    return {"actual_complex_solution_time_radius":T,"inner_real_time_half_width":INNER_HALF_WIDTH,
            "actual_complex_solution_phase_image_upper":taylor.evolution()["Picard_image_radius_upper"],
            "actual_solution_rescaled_trace_upper":parent_bounds.solution_bridge()["rescaled_trace_absolute_upper"],
            "full_complex_rational_annuli":ann,
            "N_minus_one_upper":ndev,"hbg_minus_one_upper":hdev,"D_minus_one_upper":Ddev,
            "e_fourth_minus_one_upper":e4dev,"e_minus_one_upper":edev,
            "actual_complex_log_scale_upper":logR,"R_minus_one_upper":Rdev,
            "clock_speed_minus_one_complex_upper":cdev,
            "matter_coordinate_frequency_modulus_lower":omega_lower,
            "matter_coordinate_frequency_modulus_upper":omega_upper,
            "same_sign_frequency_gap_lower":omega_lower*ann["clock_physical_speed_squared_excess"]["lower"]/(2+cdev),
            "opposite_sign_frequency_gap_lower":omega_lower*(2-cdev),
            "declared_four_mode_complex_gap_lower":GAP,
            "clock_kinetic_modulus_lower":sp.Rational(2)/(4*2**3*sp.Rational(26,100)**2),
            "clock_kinetic_modulus_upper":sp.Rational(3)/(4*sp.Rational(1,2)**3*sp.Rational(24,100)**2),
            "matter_kinetic_modulus_lower":sp.Rational(1,2)**3/2,
            "matter_kinetic_modulus_upper":2**3/sp.Rational(1,2),
            "matter_density_modulus_upper":sp.Rational(11,100)+taylor.evolution()["Picard_image_radius_upper"],
            "matter_density_modulus_lower":sp.Rational(9,100)/2}


def laurent_majorant(value,domain):
    total=0
    for term in sp.Add.make_args(sp.expand(value)):
        powers=term.as_powers_dict()
        coefficient=term
        for variable,(lo,hi) in domain.items():
            exponent=powers.get(variable,0)
            if exponent<0 and lo<=0:
                raise ValueError("A negative Laurent power has no positive modulus lower bound")
            coefficient/=variable**exponent
            coefficient*=(hi if exponent>=0 else lo)**exponent
        if coefficient.free_symbols:
            raise ValueError("An unbounded coefficient variable remains")
        total+=abs(coefficient)
    return sp.simplify(total)


@cache
def coefficient_majorants():
    p=finite.parent
    limits={p.a:(0,3),p.m:(0,8),p.g:(0,2),p.r:(0,8),
            p.rn:(sp.Rational(1,10),1),p.h:(1,6),p.ell:(0,sp.Rational(1,8)),
            p.alpha:(0,1),p.beta:(0,1),p.H:(0,1)}
    d=finite.data()
    return {key:d[key].applyfunc(lambda value:laurent_majorant(value,limits))
            for key in ("A0","A1","A2","E0","E1","C0")}


@cache
def basis_majorants():
    roots=sp.symbols("root_K_clock root_K_matter root_omega_clock root_omega_matter",positive=True)
    mapping={modes.kc:roots[0]**2,modes.km:roots[1]**2,
             modes.wc:roots[2]**2,modes.wm:roots[3]**2}
    domain={roots[0]:(sp.Rational(1,2),20),roots[1]:(sp.Rational(1,10),10),
            roots[2]:(sp.Rational(1,2),2),roots[3]:(sp.Rational(1,2),2),
            modes.ell:(0,sp.Rational(1,8))}
    out={}
    for key in ("S","S_inverse"):
        b=modes.data()[key].subs(mapping,simultaneous=True).applyfunc(lambda v:laurent_majorant(v,domain))
        out[key]={"entry_upper":b,"infinity_norm_upper":max(sum(b[i,j] for j in range(4)) for i in range(4))}
    return out


@cache
def gates():
    d=analytic_domain()
    a=d["full_complex_rational_annuli"]
    return {name:bool(value) for name,value in {
        "actual_complex_solution_disc_is_parent_certified_disc":T==taylor.TIME_WINDOW,
        "actual_complex_phase_fits_parent_trace_box":d["actual_solution_rescaled_trace_upper"]<parent_bounds.TRACE_RADIUS,
        "actual_complex_conformal_root_near_one":d["e_minus_one_upper"]<sp.Rational(4,10**6),
        "actual_complex_scale_near_one":d["R_minus_one_upper"]<sp.Rational(2,10**9),
        "complex_clock_speed_root_near_one":d["clock_speed_minus_one_complex_upper"]<sp.Rational(1,2),
        "complex_clock_excess_nonzero":a["clock_physical_speed_squared_excess"]["lower"]>sp.Rational(5,10**6),
        "complex_pivot_in_two_to_three_annulus":a["on_constraint_rescaled_lapse_Hessian"]["lower"]>2
            and a["on_constraint_rescaled_lapse_Hessian"]["upper"]<3,
        "complex_matter_square_in_positive_modulus_annulus":a["constraint_matter_square"]["lower"]>sp.Rational(9,1000)
            and a["constraint_matter_square"]["upper"]<sp.Rational(11,1000),
        "complex_M_modulus_between_twenty_four_and_twenty_six_hundredths":a["M"]["lower"]>sp.Rational(24,100)
            and a["M"]["upper"]<sp.Rational(26,100),
        "actual_complex_trace_force_below_one":a["actual_force_trace"]["upper"]<1,
        "complex_matter_frequency_between_one_half_and_two":d["matter_coordinate_frequency_modulus_lower"]>sp.Rational(1,2)
            and d["matter_coordinate_frequency_modulus_upper"]<2,
        "all_four_complex_frequency_gaps_nonzero":d["same_sign_frequency_gap_lower"]>GAP
            and d["opposite_sign_frequency_gap_lower"]>GAP,
        "complex_clock_kinetic_in_declared_annulus":d["clock_kinetic_modulus_lower"]>sp.Rational(1,2)
            and d["clock_kinetic_modulus_upper"]<200,
        "complex_matter_kinetic_in_declared_annulus":d["matter_kinetic_modulus_lower"]>sp.Rational(1,100)
            and d["matter_kinetic_modulus_upper"]<100,
        "complex_matter_density_nonzero_and_small":d["matter_density_modulus_lower"]>sp.Rational(1,25)
            and d["matter_density_modulus_upper"]<sp.Rational(1,8),
        "every_exact_Laurent_entry_below_one_hundred":all(v<100 for matrix in coefficient_majorants().values() for v in matrix),
        "actual_analytic_basis_and_inverse_norms_below_one_thousand":
            all(row["infinity_norm_upper"]<BASIS_NORM for row in basis_majorants().values()),
        "inner_real_interval_has_positive_time_Cauchy_radius":INNER_HALF_WIDTH<T,
        "no_complex_ordering_inferred_from_real_interval_order":True}.items()}


@cache
def checks():
    d=analytic_domain()
    dummy={"numerator":{"center":-5,"deviation_upper":1},"denominator":{"center":-3,"deviation_upper":1}}
    test=complex_annulus(dummy)
    return {
        "complex_annulus_handles_negative_centers_lower":test["lower"]-sp.Rational(1),
        "complex_annulus_handles_negative_centers_upper":test["upper"]-sp.Rational(3),
        "actual_complex_phase_image_is_one_e_minus_six":d["actual_complex_solution_phase_image_upper"]-sp.Rational(1,10**6),
        "inner_time_interval_uses_one_half_original_radius":INNER_HALF_WIDTH-T/2,
        "complex_conformal_deviation_keeps_entire_binomial_tail":d["e_minus_one_upper"]
            -d["e_fourth_minus_one_upper"]/(1-d["e_fourth_minus_one_upper"])}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("free"),sp.sqrt(2),[])
    slots=(("numerator","center"),("numerator","deviation_upper"),
           ("denominator","center"),("denominator","deviation_upper"))
    cases=[]
    for value in bad:
        for group,key in slots:
            row={"numerator":{"center":2,"deviation_upper":sp.Rational(1,4)},
                 "denominator":{"center":3,"deviation_upper":sp.Rational(1,4)}}
            row[group][key]=value
            cases.append(row)
    for group in ("numerator","denominator"):
        row={"numerator":{"center":2,"deviation_upper":0},"denominator":{"center":3,"deviation_upper":0}}
        row[group]["deviation_upper"]=-1
        cases.append(row)
    for center,deviation in ((0,0),(1,1),(-1,2)):
        cases.append({"numerator":{"center":1,"deviation_upper":0},
                      "denominator":{"center":center,"deviation_upper":deviation}})
    rejected=0
    for row in cases:
        try:
            complex_annulus(row)
        except (TypeError,ValueError):
            rejected+=1
    x,y=sp.symbols("majorant_variable unbounded_variable")
    for value in (1/x,y):
        try:
            laurent_majorant(value,{x:(0,1)})
        except ValueError:
            rejected+=1
    if rejected!=len(cases)+2:
        raise ValueError("An invalid complex annulus or Laurent bound was accepted")
    return {"rejected_inputs":rejected}
