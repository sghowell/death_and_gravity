"""Direct complex-time logarithmic jets on the actual nearby solution box."""
from functools import cache

import sympy as sp
from p8_proca_nearby_cones import bounds as box
from p8_proca_nearby_cones import rational as parent
from p8_proca_nearby_packets import domains

N_DOT=sp.Rational(7,50000)
Z_DOT=sp.Integer(9)
H_UPPER=sp.Rational(1,62500)
E_LOG_DOT=sp.Rational(1,12500)
PARTIAL_CAPS={
    "D":(sp.Rational(7,10**7),sp.Rational(7,10**14),sp.Integer(0)),
    "h":(sp.Rational(7,10**7),sp.Integer(0),sp.Integer(0)),
    "M":(sp.Rational(2,10**6),sp.Integer(3),sp.Integer(0)),
    "on_constraint_rescaled_lapse_Hessian":(sp.Rational(6,10**5),sp.Integer(10),sp.Rational(5,10**5)),
    "clock_physical_speed_squared":(sp.Rational(5,10**4),sp.Integer(14),sp.Rational(13,10**5))}


def require_field(value):
    if getattr(value,"field",None)!=parent.FIELD:
        raise ValueError("Require the exact actual-background rational field")
    return value


def modulus_upper(value):
    return domains.complex_annulus(box.rational_box(require_field(value)))["upper"]


def logarithmic_partial_upper(value):
    value=require_field(value)
    result=[]
    for variable in parent.FIELD.ring.gens:
        total=sp.Integer(0)
        for poly in (value.numer,value.denom):
            center=box.polynomial_box(poly)
            lower=abs(center["center"])-center["deviation_upper"]
            if lower<=0:
                raise ValueError("A logarithmic jet numerator or denominator reaches zero")
            derivative=box.polynomial_box(poly.diff(variable))
            total+=max(abs(derivative["lower"]),abs(derivative["upper"]))/lower
        result.append(total)
    return tuple(result)


@cache
def native_bounds():
    d=parent.system()
    n=modulus_upper(d["constraint_time_forcing"])/2
    z=modulus_upper(d["trace_flow_at_fixed_lapse"])+modulus_upper(d["trace_lapse_chain_coefficient"])*n
    return {"lapse_velocity_upper":n,"trace_velocity_upper":z,
        "Hubble_modulus_upper":modulus_upper(d["hat_Hubble"]),
        "conformal_log_velocity_upper":modulus_upper(d["Lt"])+modulus_upper(d["Ln"])*n,
        "logarithmic_partial_upper":{name:logarithmic_partial_upper(d[name]) for name in PARTIAL_CAPS}}


@cache
def data():
    logs={name:sum(v*c for v,c in zip((1,N_DOT,Z_DOT),caps,strict=True))
          for name,caps in PARTIAL_CAPS.items()}
    clock=(logs["D"]+logs["on_constraint_rescaled_lapse_Hessian"]
           +logs["clock_physical_speed_squared"]/2+N_DOT/sp.Rational(99,100)
           +logs["h"]+H_UPPER+2*logs["M"])/2
    matter=E_LOG_DOT+H_UPPER/2
    frequency=4*(N_DOT/sp.Rational(99,100)+E_LOG_DOT+H_UPPER
                 +logs["clock_physical_speed_squared"]/2)
    return {"actual_complex_time_radius":domains.T,
        "actual_box_radii":box.RADII,"actual_lapse_box_center":box.LAPSE_CENTER,
        "declared_lapse_velocity_upper":N_DOT,"declared_trace_velocity_upper":Z_DOT,
        "declared_Hubble_modulus_upper":H_UPPER,"declared_conformal_log_velocity_upper":E_LOG_DOT,
        "declared_partial_logarithmic_moduli":PARTIAL_CAPS,"total_logarithmic_moduli":logs,
        "clock_action_amplitude_log_velocity_triangle_upper":clock,
        "matter_action_amplitude_log_velocity_triangle_upper":matter,
        "matter_density_velocity_triangle_upper":3*H_UPPER/8,
        "signed_mode_frequency_velocity_triangle_upper":frequency,
        "declared_clock_amplitude_log_velocity_upper":sp.Rational(1,400),
        "declared_matter_amplitude_log_velocity_upper":sp.Rational(1,10000),
        "declared_density_velocity_upper":sp.Rational(3,500000),
        "declared_signed_mode_frequency_velocity_upper":sp.Rational(1,100),
        "complex_modulus_bounds_not_complex_ordering":True}


@cache
def gates():
    b,d=native_bounds(),data()
    out={"actual_lapse_velocity_below_declared":b["lapse_velocity_upper"]<N_DOT,
         "actual_trace_velocity_below_declared":b["trace_velocity_upper"]<Z_DOT,
         "actual_Hubble_modulus_below_declared":b["Hubble_modulus_upper"]<H_UPPER,
         "actual_conformal_log_velocity_below_declared":b["conformal_log_velocity_upper"]<E_LOG_DOT,
         "clock_amplitude_log_velocity_below_one_over_four_hundred":
             d["clock_action_amplitude_log_velocity_triangle_upper"]<sp.Rational(1,400),
         "matter_amplitude_log_velocity_below_one_over_ten_thousand":
             d["matter_action_amplitude_log_velocity_triangle_upper"]<sp.Rational(1,10000),
         "density_derivative_keeps_conserved_matter_charge":
             d["matter_density_velocity_triangle_upper"]==sp.Rational(3,500000),
         "signed_frequencies_have_small_actual_derivatives":
             d["signed_mode_frequency_velocity_triangle_upper"]<sp.Rational(1,100)}
    for name,caps in PARTIAL_CAPS.items():
        for j,(value,cap) in enumerate(zip(b["logarithmic_partial_upper"][name],caps,strict=True)):
            out[f"actual_partial_log_modulus_{name}_{j}_below_cap"]=value<=cap
    return {name:bool(value) for name,value in out.items()}


def checks():
    value=(2+parent.uf)/(3+parent.Nf)
    result=logarithmic_partial_upper(value)
    return {"complex_logarithmic_partial_keeps_numerator_modulus_lower":
                result[0]-1/(2-box.TIME_WINDOW),
            "complex_logarithmic_partial_keeps_denominator_modulus_lower":
                result[1]-1/(3+box.LAPSE_CENTER-box.LAPSE_RADIUS),
            "absent_trace_variable_has_exactly_zero_logarithmic_partial":result[2]}


def controls():
    bad=(True,False,sp.true,sp.false,1.0,sp.Float(1),"1",sp.oo,sp.nan,sp.zoo,None,
         sp.Symbol("unbounded"),sp.sqrt(2),[])
    count=0
    for value in bad:
        for call in (modulus_upper,logarithmic_partial_upper):
            try:
                call(value)
            except ValueError:
                count+=1
            else:
                raise ValueError("An invalid exact jet field was accepted")
    for value in (parent.FIELD.zero,parent.uf,parent.Nf-parent.FIELD.from_expr(box.LAPSE_CENTER)):
        try:
            logarithmic_partial_upper(value)
        except ValueError:
            count+=1
        else:
            raise ValueError("A zero-crossing logarithmic jet was accepted")
    try:
        modulus_upper(1/parent.uf)
    except ValueError:
        count+=1
    else:
        raise ValueError("A pole in the actual modulus box was accepted")
    return {"rejected_inputs":count}
