"""Physical energy/pressure readouts retain all rate and evaluation contacts."""
from functools import cache

import sympy as sp

from . import transport


@cache
def algebra():
    f, rate, xi = transport.rate, transport.r, transport.xi
    A, B, dA, dB, U, dU, k = sp.symbols("readout_A readout_B varied_A varied_B U varied_U k", real=True)
    S, P, Sp, Pp, dS, dP = sp.symbols("S P S_prime P_prime varied_S varied_P", real=True)
    D = k**2+U
    frequency2 = f**2*D
    physical = (dA*f*P+A*f*(dP+xi*Pp+rate*P)+dB*frequency2*S/f
                +B*(2*rate*frequency2+f**2*dU)*S/f
                +B*frequency2*(dS+xi*Sp-rate*S)/f)/2
    converted = f*((dA+rate*A)*P+(dB+rate*B)*D*S+B*dU*S
                   +A*dP+B*D*dS+xi*(A*Pp+B*D*Sp))/2
    return {"physical_variation": physical, "acoustic_variation": converted,
            "contact_only": f*((dA+rate*A)*P+(dB+rate*B)*D*S+B*dU*S
                               +xi*(A*Pp+B*D*Sp))/2,
            "checks": {"full_physical_readout_variation": sp.factor(physical-converted)}}


@cache
def checks():
    out = dict(algebra()["checks"])
    alpha, beta, z, chi, rate, n, v = sp.symbols("alpha beta z pump_source rate_source lapse_source scale_source", real=True)
    source_map = {chi: v+(alpha+beta)*n/4, rate: (1+(beta-alpha)/2)*n-v}
    kinetic = 2*chi-(1-z)*(2*v+alpha*n)-rate
    potential = -(2*chi+rate)
    out["longitudinal_acoustic_action_kinetic_vertex"] = sp.factor(
        kinetic.subs(source_map)-((-1+alpha*z)*n+(1+2*z)*v))
    out["longitudinal_acoustic_action_potential_vertex"] = sp.factor(
        potential.subs(source_map)-((-1-beta)*n-v))
    # The transverse acoustic pump is constant.
    rT = n-v
    out["transverse_acoustic_action_kinetic_vertex"] = sp.expand(-rT-(-n+v))
    out["transverse_acoustic_action_potential_vertex"] = sp.expand(
        -(rT+(1-z)*(2*v+beta*n))-((-1-beta*(1-z))*n+(2*z-1)*v))
    return out


@cache
def actual_readouts():
    from p8_vector_metric_response import source as physical

    from . import source
    result = {}
    for sector in ("T", "L"):
        item = source.data(sector)
        result[sector] = {}
        for component in ("energy", "pressure"):
            value = physical.readout(sector, component)
            A, B, dA, dB = [source.project(value[key]) for key in ("A", "B", "delta_A", "delta_B")]
            f, D, r, xi = item["rate"], item["D"], item["log_rate_variation"], item["prepared_shift"]
            factor = f/(2*source.a**3)
            baseline = factor*sp.ImmutableMatrix([[B*D, 0, A]])
            local = factor*sp.ImmutableMatrix([[(dB+r*B)*D+B*item["delta_U_fixed_u"], 0, dA+r*A]])
            history = xi*baseline*item["matrix"]
            result[sector][component] = {"response_row": baseline,
                                         "fixed_output_contact_row": local.applyfunc(source.clean),
                                         "prepared_history_row": history.applyfunc(source.clean),
                                         "scope": "Rows multiply the acoustic covariance (variance,cross,momentum variance). The varied A/B already retain physical output normalization and second mass contacts."}
    return result


@cache
def actual_checks():
    from p8_vector_metric_response import source as physical

    from . import source
    out = {}
    for sector in ("T", "L"):
        item, old = source.data(sector), physical.data(sector)
        f, D, r, xi = item["rate"], item["D"], item["log_rate_variation"], item["prepared_shift"]
        S = sp.diag(f, 1, 1/f)
        w2 = f**2*D
        dw2 = 2*source.project(old["r"])*w2
        for component in ("energy", "pressure"):
            value = physical.readout(sector, component)
            A, B, dA, dB = [source.project(value[key]) for key in ("A", "B", "delta_A", "delta_B")]
            original = sp.ImmutableMatrix([[B*w2, 0, A]])/(2*source.a**3)
            varied_original = sp.ImmutableMatrix([[dB*w2+B*dw2, 0, dA]])/(2*source.a**3)
            target = actual_readouts()[sector][component]
            out[sector+"_"+component+"_exact_physical_readout_row"] = sp.ImmutableMatrix(
                (original*S.inv()-target["response_row"]).applyfunc(sp.factor))
            converted = varied_original*S.inv()+original*(xi*S.inv()*item["matrix"]
                                  -r*sp.diag(1, 0, -1)*S.inv())
            out[sector+"_"+component+"_full_physical_contact_and_history_row"] = sp.ImmutableMatrix(
                (converted-target["fixed_output_contact_row"]-target["prepared_history_row"]).applyfunc(source.clean))
    return out
