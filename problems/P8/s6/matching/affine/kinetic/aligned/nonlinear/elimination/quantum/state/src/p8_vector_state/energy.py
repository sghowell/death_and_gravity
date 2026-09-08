"""Actual clock-sensitive lapse variation and physical canonical mode energy."""
from functools import cache

import sympy as sp

N, a, q, m2 = sp.symbols("N a_scale physical_momentum_squared m0_squared", positive=True)
alpha, beta = sp.symbols("a_N b_N", real=True)
H = sp.Symbol("H", real=True)
v, dv = sp.symbols("canonical_mode canonical_mode_dot", real=True)


@cache
def variations():
    w, dw, s, ds, t0 = sp.symbols("w w_dot sigma sigma_dot W0", real=True)
    am, bm = 1+alpha*(N-1), 1+beta*(N-1)
    LT = a*dw**2/(2*N)-N*a*(q+m2*bm)*w**2/2
    LL = a**3*(q*(ds-t0)**2/(2*N)+m2*am*t0**2/(2*N)-N*m2*bm*q*s**2/2)
    temporal = q*ds/(q+m2*am)
    rhoT = sp.factor(-sp.diff(LT, N).subs(N, 1)/a**3)
    rhoL = sp.factor((-sp.diff(LL, N)/a**3).subs(t0, temporal).subs(N, 1))
    rhoL_eliminate_first = sp.factor(-sp.diff(LL.subs(t0, temporal), N).subs(N, 1)/a**3)
    omega2 = q+m2
    expectedT = (dw**2+(q+m2*(1+beta))*w**2)/(2*a**2)
    expectedL = (m2*q/omega2*(1-alpha*q/omega2)*ds**2+m2*q*(1+beta)*s**2)/2
    k2 = sp.Symbol("k_com_squared", positive=True)
    # Pressure varies the physical scale factor at fixed comoving momentum.
    pT = sp.factor(sp.diff(LT.subs({N: 1, q: k2/a**2}), a).subs(k2, q*a**2)/(3*a**2))
    pL = sp.factor(sp.diff(LL.subs({N: 1, q: k2/a**2}), a).subs(k2, q*a**2)
                   .subs(t0, temporal.subs(N, 1))/(3*a**2))
    return {"w": w, "dw": dw, "sigma": s, "dsigma": ds, "W0": t0,
            "rhoT": rhoT, "rhoL": rhoL, "pT": pT, "pL": pL,
            "full_temporal_constraint": sp.factor(sp.diff(LL, t0).subs(t0, temporal)),
            "lapse_variation_and_temporal_elimination_commute": sp.factor(rhoL-rhoL_eliminate_first),
            "clock_sensitive_transverse_lapse_energy": sp.factor(rhoT-expectedT),
            "clock_sensitive_longitudinal_lapse_energy": sp.factor(rhoL-expectedL)}


@cache
def canonical():
    d = variations()
    omega2, z = q+m2, q/(q+m2)
    rateT, rateL = H/2, H*(sp.Rational(1, 2)+z)
    gL2 = a**3*m2*q/omega2
    mapT = {d["w"]: v/sp.sqrt(a), d["dw"]: (dv-rateT*v)/sp.sqrt(a)}
    mapL = {d["sigma"]: v/sp.sqrt(gL2), d["dsigma"]: (dv-rateL*v)/sp.sqrt(gL2)}
    rhoT = ((dv-rateT*v)**2+(omega2+beta*m2)*v**2)/(2*a**3)
    rhoL = ((1-alpha*z)*(dv-rateL*v)**2+omega2*(1+beta)*v**2)/(2*a**3)
    pT = ((dv-rateT*v)**2+(q-m2)*v**2)/(6*a**3)
    pL = ((1+2*z)*(dv-rateL*v)**2-omega2*v**2)/(6*a**3)
    return {"z": z, "omega2": omega2, "rateT": rateT, "rateL": rateL,
            "rhoT": rhoT, "rhoL": rhoL, "pT": pT, "pL": pL,
            "physical_transverse_canonical_energy": sp.factor(d["rhoT"].subs(mapT)-rhoT),
            "physical_longitudinal_canonical_energy": sp.factor(d["rhoL"].subs(mapL)-rhoL),
            "physical_transverse_canonical_pressure": sp.factor(d["pT"].subs(mapT)-pT),
            "physical_longitudinal_canonical_pressure": sp.factor(d["pL"].subs(mapL)-pL)}


@cache
def actual_coefficients():
    from p8_vector_bubble import pole
    d = pole.actual_mass_direction()
    return {"h": d["h"], "a_N": d["alpha"], "b_N": d["beta"]}


@cache
def zero_order_pole():
    """Independent zero-point energy pole including the lapse derivative.

    I=1/2 integral d^3p/(2pi)^3 sqrt(p²+m²) has pole
    -m^4/(64pi² epsilon_DR). Differentiating in m² fixes the
    integral of 1/omega; momentum moments then reduce algebraically.
    """
    I, Iprime = sp.symbols("scalar_zero_point scalar_zero_point_mass_derivative", real=True)
    transverse = I+beta*m2*Iprime
    longitudinal = I+(beta-alpha)*I/2+alpha*m2*Iprime
    combined = sp.expand(2*transverse+longitudinal)
    relative_pole = sp.factor(combined.subs({I: 1, Iprime: 2/m2}))
    from p8_aligned_quantum import potential
    d = actual_coefficients()
    jets = potential.clock_jets()
    actual = relative_pole.subs({alpha: d["a_N"], beta: d["b_N"]})
    expected = jets["pole_weight"]["value"]+jets["pole_weight"]["N_first"]
    expected = expected.subs(jets["h"], d["h"])
    epsilon = sp.Symbol("epsilon_DR", positive=True)
    dimensional_scalar = (m2**(2-epsilon)*sp.gamma(epsilon-2)
                          /(2*(4*sp.pi)**(sp.Rational(3, 2)-epsilon)*sp.gamma(-sp.Rational(1, 2))))
    return {"relative_pole": relative_pole,
            "direct_spatial_dimensional_zero_point_residue": sp.simplify(
                sp.limit(epsilon*dimensional_scalar, epsilon, 0)+m2**2/(64*sp.pi**2)),
            "full_lapse_local_pole_value_plus_N_derivative": sp.factor(actual-expected),
            "omitting_clock_variation_pole_difference": sp.factor(actual-3)}


@cache
def checks():
    d, c = variations(), canonical()
    out = {name: d[name] for name in ("full_temporal_constraint", "lapse_variation_and_temporal_elimination_commute",
                                    "clock_sensitive_transverse_lapse_energy", "clock_sensitive_longitudinal_lapse_energy")}
    out.update({name: value for name, value in c.items() if name.startswith("physical_")})
    out["full_lapse_local_pole_value_plus_N_derivative"] = zero_order_pole()["full_lapse_local_pole_value_plus_N_derivative"]
    out["direct_spatial_dimensional_zero_point_residue"] = zero_order_pole()["direct_spatial_dimensional_zero_point_residue"]
    return out
