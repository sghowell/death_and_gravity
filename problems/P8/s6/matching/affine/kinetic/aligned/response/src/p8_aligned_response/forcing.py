"""Derive, rather than prescribe, the leading canonical source and readouts."""
from functools import cache

import sympy as sp
from p8_affine_aligned import modes

u = modes.u


@cache
def temporal():
    t, S, sigma, sigmad = sp.symbols("t S2 sigma sigma_dot", real=True)
    q, zeta = sp.symbols("q zeta", positive=True)
    r = zeta*q
    L = (t-S)**2/2-q*sigma**2/2+r*(sigmad-t)**2/2
    solution = (S+r*sigmad)/(1+r)
    target = r*(sigmad-S)**2/(2*(1+r))-q*sigma**2/2
    return {"t": t, "S": S, "sigma": sigma, "sigmad": sigmad, "q": q, "zeta": zeta,
            "L": L, "solution": solution, "reduced": target,
            "Euler_residual": sp.factor(sp.diff(L, t).subs(t, solution)),
            "action_residual": sp.factor(L.subs(t, solution)-target),
            "temporal_readout_residual": sp.factor(solution-S-r*(sigmad-S)/(1+r))}


@cache
def canonical():
    g = sp.Function("g")(u)
    source = sp.Function("S2")(u)
    v = sp.Function("v")(u)
    omega0 = sp.Symbol("q_plus_mass", positive=True)
    sigma = v/g
    original_Euler = sp.diff(g**2*(sp.diff(sigma, u)-source), u)+g**2*omega0*sigma
    force = sp.diff(g**2*source, u)/g
    canonical_Euler = sp.diff(v, u, 2)+(omega0-sp.diff(g, u, 2)/g)*v-force
    expected = g*sp.diff(source, u)+2*sp.diff(g, u)*source
    return {"g": g, "source": source, "force": force,
            "equation_residual": sp.factor(original_Euler/g-canonical_Euler),
            "force_residual": sp.expand(force-expected),
            "force_prime": g*sp.diff(source, u, 2)+3*sp.diff(g, u)*sp.diff(source, u)+2*sp.diff(g, u, 2)*source,
            "force_second": g*sp.diff(source, u, 3)+4*sp.diff(g, u)*sp.diff(source, u, 2)
                            +5*sp.diff(g, u, 2)*sp.diff(source, u)+2*sp.diff(g, u, 3)*source}


@cache
def checks():
    aux, data = temporal(), canonical()
    source = data["source"]
    return {"full_temporal_Euler": aux["Euler_residual"],
            "full_temporal_action": aux["action_residual"],
            "temporal_physical_readout": aux["temporal_readout_residual"],
            "actual_canonical_equation": data["equation_residual"],
            "actual_canonical_force": data["force_residual"],
            "canonical_force_prime": sp.expand(sp.diff(data["force"], u)-data["force_prime"]),
            "canonical_force_second": sp.expand(sp.diff(data["force"], u, 2)-data["force_second"]),
            "source_preparation_gives_zero_force": data["force"].subs({source: 0, sp.diff(source, u): 0}, simultaneous=True),
            "source_preparation_gives_zero_force_prime": data["force_prime"].subs(
                {source: 0, sp.diff(source, u): 0, sp.diff(source, u, 2): 0}, simultaneous=True),
            "physical_spatial_readout": sp.factor(modes.canonical()["q"]
                /modes.canonical()["longitudinal"]["weight"]
                -(1+modes.canonical()["r"])/(modes.canonical()["a"]**3*modes.ZETA))}
