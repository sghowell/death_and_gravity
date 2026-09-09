"""Continuous fixed-mass Proca stress bounds through five time derivatives."""
from functools import cache

import sympy as sp
from p8_constant_proca import quantum
from p8_vector_regularity import estimates as previous
from p8_vector_state import comparison, wkb

from . import readouts, subtraction

MASS=sp.Integer(1000)
SCALE=sp.Integer(10**400)


def exact_scale(value):
    if isinstance(value,bool) or not isinstance(value,(int,sp.Rational)):
        raise TypeError("Require an exact positive M*tau scale")
    value=sp.Rational(value)
    if value<=0:
        raise ValueError("Require an exact positive M*tau scale")
    return value


def local_derivatives(order):
    readouts.validate("transverse",order)
    return _local_derivatives(order)


@cache
def _local_derivatives(order):
    values={n:sp.factor(sp.diff(row["energy"],wkb.u,order))
            for n,row in quantum.local_coefficients().items()}
    boxes={n:wkb.box_bound(value) for n,value in values.items()}
    return {"actual_new_local_energy_derivatives":values,
            "absolute_envelopes":{n:row["absolute_upper"] for n,row in boxes.items()},
            "reconstructions":{n:row["reconstruction"] for n,row in boxes.items()}}


def physical_bounds(scale):
    return _physical_bounds(exact_scale(scale))


@cache
def _physical_bounds(scale):
    old=previous.physical_bounds(scale,MASS)
    A=comparison.AMAX
    radial=previous.radial_envelope(MASS)
    J5=radial["integrated_beta_times_frequency_sixth_upper"]
    energy={}
    for order in range(6):
        E=max(readouts.row_bounds(kind,order)["normalized_reference_product_envelope"]
              for kind in ("transverse","longitudinal"))
        tails=[subtraction.reference_tail(kind,order) for kind in ("transverse","longitudinal")]
        if any(not row["integrable_tail_certified"] or not row["no_lower_Laurent_power"] for row in tails):
            raise ValueError("A nonintegrable new energy subtraction coefficient survived")
        T=max(row["reference_bracket_tail_integer_upper"] for row in tails)
        difference=9*E*A**(order+1)*MASS**(order-5)*J5/scale**2
        reference=T/(72*MASS**2*scale**2)
        local=sum(value*MASS**(4-2*n) for n,value in local_derivatives(order)["absolute_envelopes"].items())/(576*scale**2)
        energy[order]={"exact_mode_minus_projected_reference":difference,
                       "projected_reference_minus_actual_differentiated_subtraction":reference,
                       "new_ordinary_finite_local_matched_term":local,
                       "total":difference+reference+local}
    pressure=old["normalized_derivative_bounds"]["pressure"]
    delta={order:energy[order]["total"]+old["normalized_derivative_bounds"]["energy"][order]["total"]
           for order in range(6)}
    return {"M_tau":scale,"fixed_canonical_mass_time_product":MASS,
            "unchanged_initial_state_radial_envelope":radial,
            "new_energy_derivative_bounds":energy,
            "identical_pressure_derivative_bounds":pressure,
            "new_minus_old_energy_derivative_upper":delta,
            "new_minus_old_pressure_derivative_exactly_zero":True,
            "normalization":"Derivative j is divided by M^2/tau^(2+j), for j=0,...,5. The same actual continuum state is used; no momentum grid replaces the radial integral."}


@cache
def gates():
    d=physical_bounds(SCALE)
    return {"every_new_energy_and_unchanged_pressure_derivative_below_one_e_minus_770":all(
                max(d["new_energy_derivative_bounds"][j]["total"],d["identical_pressure_derivative_bounds"][j]["total"])<sp.Rational(1,10**770)
                for j in range(6)),
            "every_energy_profile_change_derivative_below_one_e_minus_770":all(
                value<sp.Rational(1,10**770) for value in d["new_minus_old_energy_derivative_upper"].values()),
            "fixed_original_mass_and_named_classical_scale_kept":d["fixed_canonical_mass_time_product"]==MASS and d["M_tau"]==SCALE,
            "no_new_functional_metric_or_mixed_loop_estimate_from_C5_profile_bound":True}
