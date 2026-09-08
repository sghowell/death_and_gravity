"""Independent physical Hamiltonian vertices and moving canonical map."""
from functools import cache

import sympy as sp

from . import jets


def data(sector):
    return _data(jets.kind(sector))


@cache
def _data(sector):
    D, z, a, b, a2, b2, n, v = (jets.D, jets.z, jets.alpha[0], jets.beta[0],
                               jets.alpha2[0], jets.beta2[0], jets.n[0], jets.v[0])
    if sector == "T":
        r = (1+b*(1-z)/2)*n-z*v
        rg = ((D-2)*v-n)/2
        weights = {"N": (-sp.Integer(1), -1-b*(1-z)),
                   "Z": (D-2, 2*z-D+2)}
        varied = {"N": (n, (1-b2*(1-z)+b**2*(1-z)**2)*n-2*b*z*(1-z)*v),
                  "Z": (sp.Integer(0), -4*z*(1-z)*v-2*b*z*(1-z)*n)}
    else:
        r = (1+(b-a*z)/2)*n-z*v
        rg = ((D-2+2*z)*v+(a*z-1)*n)/2
        weights = {"N": (-1+a*z, -1-b),
                   "Z": (D-2+2*z, -(D-2))}
        varied = {"N": ((1+a2*z-a**2*z*(2-z))*n-2*a*z*(1-z)*v, (1-b2+b**2)*n),
                  "Z": (-4*z*(1-z)*v-2*a*z*(1-z)*n, sp.Integer(0))}
    return {"delta_log_frequency": r, "delta_log_g": rg,
            "weights": weights, "varied_weights": varied}


@cache
def checks():
    N, a, q, m2 = sp.symbols("physical_lapse physical_scale physical_q mass_squared", positive=True)
    e = sp.Symbol("metric_parameter", real=True)
    am = 1+jets.alpha[0]*(N-1)+jets.alpha2[0]*(N-1)**2/2
    bm = 1+jets.beta[0]*(N-1)+jets.beta2[0]*(N-1)**2/2
    out = {}
    def scale(value):
        return a*sp.diff(value, a)-2*q*sp.diff(value, q)
    def change(value):
        return jets.n[0]*sp.diff(value, N)+jets.v[0]*scale(value)
    def project(value):
        return sp.factor(value.subs(N, 1).subs(q, m2*jets.z/(1-jets.z)))
    for sector in ("T", "L"):
        if sector == "T":
            g2, bare = a**(jets.D-2)/N, N**2*(q+m2*bm)
        else:
            g2 = a**jets.D*m2*am*q/(N*(q+m2*am))
            bare = N**2*bm*(q/am+m2)
        target = data(sector)
        out[sector+"_physical_log_g"] = sp.factor(project(change(g2)/(2*g2))-target["delta_log_g"])
        out[sector+"_physical_log_frequency"] = sp.factor(project(change(bare)/(2*bare))-target["delta_log_frequency"])
        for label, derivative in (("N", lambda value: sp.diff(value, N)), ("Z", scale)):
            A = sp.factor(derivative(g2)/g2)
            B = sp.factor(-derivative(g2*bare)/(g2*bare))
            for j, value in enumerate((A, B)):
                out[sector+"_"+label+"_physical_weight_"+str(j)] = sp.factor(project(value)-target["weights"][label][j])
                out[sector+"_"+label+"_physical_contact_"+str(j)] = sp.factor(project(change(value))-target["varied_weights"][label][j])
        perturbed = g2.subs({N: 1+e*jets.n[0], a: a*(1+e*jets.v[0]), q: q/(1+e*jets.v[0])**2}, simultaneous=True)
        rate = (jets.time(perturbed)+jets.H[0]*scale(perturbed))/(2*perturbed)
        out[sector+"_physical_rate_variation"] = sp.factor(project(sp.diff(rate, e).subs(e, 0))-jets.time(target["delta_log_g"]))
    return out
