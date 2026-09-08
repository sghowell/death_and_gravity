"""Mass-only Hamiltonian insertion and temporal-constraint contact."""
from functools import cache

import sympy as sp
from p8_vector_variation import modes as prior

a, m, q = sp.symbols("scale mass momentum_squared", positive=True)
alpha, beta, n = sp.symbols("alpha beta mass_source", real=True)
v, p, sigma, momentum = sp.symbols("v p sigma canonical_momentum", real=True)


def kind(value):
    if type(value) is not str or value not in ("T", "L"):
        raise ValueError("Require the physical T or L polarization sector")
    return value


def clean(value):
    return value.applyfunc(sp.factor) if isinstance(value, sp.MatrixBase) else sp.factor(value)


@cache
def hamiltonians():
    am, bm = 1+alpha*n, 1+beta*n
    return {"T": momentum**2/(2*a)+a*(q+m**2*bm)*sigma**2/2,
            "L": momentum**2/(2*a**3*q)+momentum**2/(2*a**3*m**2*am)
                 +a**3*m**2*bm*q*sigma**2/2}


def insertion(sector):
    sector = kind(sector)
    return _insertion(sector)


@cache
def _insertion(sector):
    omega2 = m**2+q
    g2 = a if sector == "T" else a**3*m**2*q/omega2
    replace = {sigma: v/sp.sqrt(g2), momentum: sp.sqrt(g2)*p}
    H = hamiltonians()[sector]
    J = clean(-sp.diff(H, n).subs(n, 0).subs(replace))
    contact = clean(-sp.diff(H, n, 2).subs(n, 0).subs(replace))
    return {"J": J, "contact": contact, "J_matrix": sp.ImmutableMatrix(sp.hessian(J, (v, p))),
            "g_squared": g2}


@cache
def constraint():
    temporal, divergence = sp.symbols("temporal divergence_of_momentum", real=True)
    am = sp.Symbol("a_mass", positive=True)
    part = -temporal*divergence-a**3*m**2*am*temporal**2/2
    solution = -divergence/(a**3*m**2*am)
    reduced = sp.factor(part.subs(temporal, solution))
    return {"solution": solution, "reduced_constraint_energy": reduced,
            "Euler_identity": sp.factor(sp.diff(part, temporal).subs(temporal, solution)),
            "energy_identity": reduced-divergence**2/(2*a**3*m**2*am),
            "local_contact_density": -alpha**2*divergence**2/(a**3*m**2),
            "source_contact_identity": sp.factor(
                -sp.diff(reduced.subs(am, 1+alpha*n), n, 2).subs(n, 0)
                +alpha**2*divergence**2/(a**3*m**2))}


@cache
def checks():
    out = {name: value for name, value in constraint().items() if name.endswith("identity")}
    targets = {"T": -m**2*beta*v**2/2,
               "L": (alpha*q*p**2/(m**2+q)-beta*(m**2+q)*v**2)/2}
    for sector, long_name in (("T", "transverse"), ("L", "longitudinal")):
        data = insertion(sector)
        out[sector+"_mass_insertion"] = clean(data["J"]-targets[sector])
        target_contact = 0 if sector == "T" else -alpha**2*q*p**2/(m**2+q)
        out[sector+"_contact"] = clean(data["contact"]-target_contact)
        old = prior.canonical()[long_name]
        replace = {prior.N: 1, prior.scale: a, prior.mass2: m**2,
                   prior.q: q, prior.am: 1+alpha*n, prior.bm: 1+beta*n}
        g2 = old["g_squared"].subs(replace)
        bare = old["bare_frequency_squared"].subs(replace)
        H = momentum**2/(2*g2)+g2*bare*sigma**2/2
        out[sector+"_frozen_physical_Hamiltonian"] = clean(H-hamiltonians()[sector])
        g0 = sp.sqrt(data["g_squared"])
        C0 = sp.diag(g0, 1/g0)
        physical_M = sp.Matrix([[0, 1/g2], [-g2*bare, 0]])
        variation = C0*physical_M.diff(n).subs(n, 0)*C0.inv()
        symplectic = sp.Matrix([[0, 1], [-1, 0]])
        out[sector+"_fixed_clock_map_source_generator"] = clean(variation+symplectic*data["J_matrix"])
    return out
