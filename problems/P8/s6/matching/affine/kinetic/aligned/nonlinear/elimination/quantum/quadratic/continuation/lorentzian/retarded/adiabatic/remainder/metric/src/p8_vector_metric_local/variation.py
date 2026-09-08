"""Full metric variation of the dimensional adiabatic readouts."""
from functools import cache

import sympy as sp
from p8_vector_dimensional import local

from . import canonical, jets


def data(sector):
    return _data(jets.kind(sector))


@cache
def _data(sector):
    old = local.reference("transverse" if sector == "T" else "longitudinal")
    d, P2, P4, c1 = [old[key] for key in ("rate", "P2", "P4", "c1")]
    lam = -jets.H[0]*jets.z
    physical = canonical.data(sector)
    r, rg = physical["delta_log_frequency"], physical["delta_log_g"]
    dd, dl = jets.time(rg), jets.time(r)
    dU = jets.time(dd)+2*d*dd
    dP2 = sp.expand(-dU/2-jets.time(dl)/4+lam*dl/4)
    dP4 = sp.expand(-P2*dP2-jets.time(jets.time(dP2))/4
                    +sp.Rational(5, 4)*(dl*jets.time(P2)+lam*jets.time(dP2))
                    +(jets.time(dl)/2-3*lam*dl)*P2
                    +(jets.time(lam)/2-sp.Rational(3, 2)*lam**2)*dP2)
    b2 = jets.time(P2)-2*lam*P2
    db2 = jets.time(dP2)-2*dl*P2-2*lam*dP2
    dc1 = dd+dl/2
    coefficients, baseline = {}, {}
    for label in ("N", "Z"):
        A, B = physical["weights"][label]
        dA, dB = physical["varied_weights"][label]
        rows = {0: A+B, 1: (A-B)*P2+A*c1**2,
                2: (A-B)*P4+B*P2**2+A*(c1*b2-c1**2*P2)}
        changes = {
            0: dA+dB,
            1: (dA-dB)*P2+(A-B)*dP2+dA*c1**2+2*A*c1*dc1,
            2: (dA-dB)*P4+(A-B)*dP4+dB*P2**2+2*B*P2*dP2
               +dA*(c1*b2-c1**2*P2)
               +A*(dc1*b2+c1*db2-2*c1*dc1*P2-c1**2*dP2)}
        baseline[label] = rows
        coefficients[label] = {j: sp.expand(changes[j]+(1-2*j)*r*rows[j]) for j in (0, 1, 2)}
    return {"delta_log_frequency": r, "delta_log_g": rg, "delta_U": sp.expand(dU),
            "delta_P2": dP2, "delta_P4": dP4, "baseline": baseline, "coefficients": coefficients}


def combined(output, order):
    output, order = jets.output(output), jets.order(order)
    return sp.expand((jets.D-1)*data("T")["coefficients"][output][order]+data("L")["coefficients"][output][order])


@cache
def checks():
    out = {}
    for sector in ("T", "L"):
        old = local.reference("transverse" if sector == "T" else "longitudinal")
        new = data(sector)
        P2, P4 = old["P2"], old["P4"]
        r, dU, dP2, dP4 = [new[key] for key in ("delta_log_frequency", "delta_U", "delta_P2", "delta_P4")]
        dl, lam = jets.time(r), -jets.H[0]*jets.z
        dS1, dS2 = dP2-2*r*P2, dP4-4*r*P4
        dL1 = jets.time(dS1)-2*lam*dS1
        b2 = jets.time(P2)-2*lam*P2
        out[sector+"_full_metric_Riccati_order_two"] = sp.expand(-2*dS1-4*r*P2-dU-jets.time(dl)/2+lam*dl/2)
        out[sector+"_full_metric_Riccati_order_four"] = sp.expand(
            -2*(dS2+P2*dS1)-2*r*(P2**2+2*P4)
            -(jets.time(dL1)-2*lam*dL1)/2+(lam*dL1+b2*dl)/2)
        for label in ("N", "Z"):
            for j in (0, 1, 2):
                value = new["coefficients"][label][j]
                out[sector+"_"+label+"_linear_"+str(j)] = sp.expand(
                    sum(field*sp.diff(value, field) for field in jets.n+jets.v)-value)
        out[sector+"_energy_baseline"] = sp.expand(new["baseline"]["N"][2]+local.energy_coefficients()[sector_name(sector)][2])
        out[sector+"_pressure_baseline"] = sp.expand(new["baseline"]["Z"][2]-jets.D*local.pressure_coefficients()[sector_name(sector)][2])
    return out


def sector_name(sector):
    return "transverse" if jets.kind(sector) == "T" else "longitudinal"
