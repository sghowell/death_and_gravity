"""Independent D-spatial-dimensional lapse and isotropic metric variations."""
from functools import cache

import sympy as sp

from . import local


@cache
def variations():
    D = local.dimension
    N, a, q, mass2 = sp.symbols("N a_scale q m_squared", positive=True)
    w, dw, sigma, ds, temporal = sp.symbols("w w_dot sigma sigma_dot W0", real=True)
    am, bm = 1+local.alpha*(N-1), 1+local.beta*(N-1)
    LT = a**(D-2)*(dw**2/N-N*(q+mass2*bm)*w**2)/2
    LL = a**D*(q*(ds-temporal)**2/N+mass2*am*temporal**2/N-N*mass2*bm*q*sigma**2)/2
    eliminated = q*ds/(q+mass2*am)
    rhoT = sp.factor(-sp.diff(LT, N).subs(N, 1)/a**D)
    rhoL = sp.factor(-sp.diff(LL, N).subs(temporal, eliminated).subs(N, 1)/a**D)
    omega2, z = mass2+q, q/(mass2+q)
    expectedT = (dw**2+(q+mass2*(1+local.beta))*w**2)/(2*a**2)
    expectedL = (mass2*q/omega2*(1-local.alpha*z)*ds**2+mass2*q*(1+local.beta)*sigma**2)/2
    k2 = sp.Symbol("k_com_squared", positive=True)
    pT = sp.factor(sp.diff(LT.subs({N: 1, q: k2/a**2}), a).subs(k2, q*a**2)/(D*a**(D-1)))
    pL = sp.factor(sp.diff(LL.subs({N: 1, q: k2/a**2}), a).subs(k2, q*a**2)
                   .subs(temporal, eliminated.subs(N, 1))/(D*a**(D-1)))
    expectedPT = ((D-2)*dw**2+((4-D)*q-(D-2)*mass2)*w**2)/(2*D*a**2)
    expectedPL = (mass2*q/omega2*(D-2+2*z)*ds**2-(D-2)*mass2*q*sigma**2)/(2*D)
    gT2, gL2 = a**(D-2), a**D*mass2*q/omega2

    def log_rate(square):
        return sp.factor((a*local.H*sp.diff(square, a)-2*local.H*q*sp.diff(square, q))/(2*square))

    return {"D_dimensional_temporal_constraint": sp.factor(sp.diff(LL, temporal).subs(temporal, eliminated)),
            "D_dimensional_transverse_lapse_energy": sp.factor(rhoT-expectedT),
            "D_dimensional_longitudinal_lapse_energy": sp.factor(rhoL-expectedL),
            "D_dimensional_transverse_pressure": sp.factor(pT-expectedPT),
            "D_dimensional_longitudinal_pressure": sp.factor(pL-expectedPL),
            "D_dimensional_transverse_canonical_rate": sp.factor(log_rate(gT2)-local.reference("transverse")["rate"]),
            "D_dimensional_longitudinal_canonical_rate": sp.factor(log_rate(gL2)-local.reference("longitudinal")["rate"].subs(local.z, z))}
