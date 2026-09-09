"""Continuous C4-to-C0 bounds for the new finite local response only."""
from functools import cache

import sympy as sp
from p8_vector_state import wkb

from . import chart, local

MASS=sp.Integer(1000)
SCALE=sp.Integer(10**400)


def exact_scale(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact positive M*tau scale")
    value=sp.Rational(value)
    if value<=0:
        raise ValueError("Require an exact positive M*tau scale")
    return value


def coefficient_bounds(value,rows):
    fields=tuple(field for row in rows.values() for field in row)
    polynomial=sp.Poly(value,*fields)
    if any(sum(powers)!=1 for powers,c in polynomial.terms() if c!=0):
        raise ValueError("Require a source-linear local response")
    if value.free_symbols-set(fields)-{local.u}:
        raise ValueError("An unknown source or unfixed background coefficient remains")
    tables,checks={},{}
    for label,row in rows.items():
        tables[label]={}
        for j,field in enumerate(row[:5]):
            coefficient=sp.factor(sp.diff(value,field))
            box=wkb.box_bound(coefficient)
            tables[label][j]=box["absolute_upper"]
            checks[label+str(j)]=box["reconstruction"]
        if any(sp.diff(value,field)!=0 for field in row[5:]):
            raise ValueError("The C4 source norm does not control an actual higher derivative")
    return {"coefficient_absolute_envelopes":tables,"box_reconstructions":checks,
            "joint_source_norm_coefficient_upper":sum(sum(row.values()) for row in tables.values())}


@cache
def envelopes():
    physical={component:{order:coefficient_bounds(local.physical_operator(component,order),
                        {"N":local.n,"Z":local.v}) for order in range(3)}
              for component in ("energy","pressure")}
    clock={component:{order:coefficient_bounds(chart.operator(component,order),
                    {"N":local.n,"V":chart.vhat}) for order in range(3)}
           for component in ("N","V")}
    return {"physical_finite_local_response":physical,
            "clock_chart_finite_local_plus_its_fixed_profile":clock}


def at_scale(scale):
    return _at_scale(exact_scale(scale))


@cache
def _at_scale(scale):
    result={}
    for chart_name,outputs in envelopes().items():
        result[chart_name]={}
        for output,rows in outputs.items():
            terms={order:MASS**(4-2*order)*d["joint_source_norm_coefficient_upper"]/(576*scale**2)
                   for order,d in rows.items()}
            result[chart_name][output]={"adiabatic_half_order_bounds":terms,"C4_to_C0_upper":sum(terms.values())}
    return {"M_tau":scale,"fixed_mass_time_product":MASS,"response_bounds":result,
            "source_norm":"Maximum of both source functions and their coordinate-time derivatives 0 through 4 on |u|<=1/2, in the declared physical or clock chart.",
            "normalization":"Dimensionless local Euler/stress response divided by M^2/tau^2; pi^2>9 supplies denominator 576.",
            "full_nonlocal_or_mixed_loop_bound":False}


@cache
def gates():
    d=at_scale(SCALE)
    return {"all_new_local_response_component_bounds_below_one_e_minus_790":all(
                row["C4_to_C0_upper"]<sp.Rational(1,10**790)
                for outputs in d["response_bounds"].values() for row in outputs.values()),
            "same_fixed_mass_and_named_scale":d["fixed_mass_time_product"]==1000 and d["M_tau"]==10**400,
            "zero_order_constant_heat_term_cancels_its_fixed_profile":all(
                row["adiabatic_half_order_bounds"][0]==0
                for row in d["response_bounds"]["clock_chart_finite_local_plus_its_fixed_profile"].values()),
            "finite_local_response_not_the_new_nonlocal_gaussian_kernel":True,
            "coordinate_C4_bound_does_not_establish_a_coupled_causal_inverse":True}
