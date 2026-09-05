"""Independent quadratic bridges and explicit failures of omitted terms."""

from functools import cache
from itertools import combinations_with_replacement

import sympy as sp
from p8_physical.momentum import residual as vacuum_constraint

from . import background, quadratic
from .vertices import Leg, construction, hamiltonian_kernel, tensor_basis

TRI = ((2, 3, 0), (0, 2, 4), (-2, -5, -4))
QUAD = ((2, 3, 0), (0, 2, 4), (3, -2, 1), (-5, -3, -5))


def fixture(kinds):
    waves = TRI if len(kinds) == 3 else QUAD
    if len(kinds) not in (3, 4):
        raise ValueError("Use three or four fixture legs")
    return tuple(Leg(k, kind, tensor_basis(k)[0] if kind in ("t", "t_dot", "pi") else None)
                 for k, kind in zip(waves, kinds))


def encode(legs):
    return [{"wave": list(map(str, leg.wave)), "kind": leg.kind,
             "polarization": [list(map(str, row)) for row in leg.polarization]
             if leg.polarization else None} for leg in legs]


@cache
def quadratic_checks(time_point=None):
    wave = (2, 2, 0)
    opposite = tuple(-k for k in wave)
    q = sum(k*k for k in wave)
    kinds = ("s", "m", "p", "P")
    variables = (*quadratic.Q, *quadratic.P)
    out = {}
    for chart in ("unitary", "gamma"):
        target = quadratic.target(time_point, chart, q)
        for i, j in combinations_with_replacement(range(4), 2):
            legs = (Leg(wave, kinds[i]), Leg(opposite, kinds[j]))
            expected = sp.diff(target, variables[i], variables[j])
            out[f"{chart}_{kinds[i]}_{kinds[j]}"] = sp.cancel(
                hamiltonian_kernel(legs, time_point, chart)["kernel"]-expected)
    E = tensor_basis(wave)[0]
    norm = sp.trace(sp.Matrix(E)**2)
    entries = {("t", "t"): q*norm/4, ("pi", "pi"): 4/norm, ("t", "pi"): 0,
               **{(kind, tensor_kind): 0 for kind in kinds for tensor_kind in ("t", "pi")}}
    for pair, target in entries.items():
        legs = tuple(Leg(k, kind, E if kind in ("t", "pi") else None)
                     for k, kind in zip((wave, opposite), pair))
        out[f"tensor_{'_'.join(pair)}"] = sp.cancel(hamiltonian_kernel(legs, time_point)["kernel"]-target)
    for kind in kinds:
        out[f"tadpole_{kind}"] = hamiltonian_kernel((Leg((0, 0, 0), kind),), time_point)["kernel"]
    return out


def boundary_checks():
    a, H, l, adot, ldot = sp.symbols("a H l adot ldot", real=True)
    v, s, dv, ds = sp.symbols("v s dv ds", real=True)
    p, pm = sp.symbols("p pm", real=True)
    F = 3*a**3*l*v*s
    shift = (p+sp.diff(F, v))*dv+(pm+sp.diff(F, s))*ds
    residual = shift-p*dv-pm*ds-sp.diff(F, v)*dv-sp.diff(F, s)*ds
    generator = sp.diff(F, a)*adot+sp.diff(F, l)*ldot
    return {"mixed_boundary_symplectic": sp.expand(residual),
            "mixed_boundary_time_generator": sp.expand(generator.subs({adot: H*a, ldot: -3*H*l})),
            "compact_q_drift": sp.cancel((1-background.rtime**2)/2*sp.diff(
                (1-background.functions()["x"]**2)**3,
                                               background.rtime)
                +6*background.functions()["x"]*(1-background.functions()["x"]**2)**3)}


def negative_controls():
    point = sp.Rational(1, 3)
    bg = background.functions(point)
    tadpole_legs = (Leg((0, 0, 0), "P"),)
    wrong_clock = construction(tadpole_legs, point, subtract_clock_velocity=False)
    clock_error = wrong_clock["hamiltonian"].coefficient(1)
    if sp.cancel(clock_error-bg["l"]) != 0 or clock_error == 0:
        raise ValueError("Missing background matter clock failed its tadpole control")
    pair = (Leg((2, 2, 0), "p"), Leg((-2, -2, 0), "m"))
    correct = construction(pair, point)
    wrong_boundary = construction(pair, point, mixed_boundary=False)
    boundary_error = sp.factor(wrong_boundary["hamiltonian"].coefficient(3)
                               -correct["hamiltonian"].coefficient(3))
    if boundary_error == 0:
        raise ValueError("Omitted mixed matter boundary was not detected")
    # A vacuum York equation is not a constraint of the rolling-matter model.
    vacuum = vacuum_constraint(correct["momentum"]["momentum"], correct["geometry"]["christoffel"])
    vacuum_errors = [sp.factor(value.coefficient(2)) for value in vacuum]
    if not any(vacuum_errors):
        raise ValueError("Matter-free spatial constraint was not detected")
    legs = fixture(("s", "m", "s", "m"))
    full = hamiltonian_kernel(legs, 0)["kernel"]
    truncated = construction(legs, 0, constraint_order=1)
    truncation_error = sp.factor(truncated["hamiltonian"].coefficient(15)-full)
    if truncation_error == 0:
        raise ValueError("Linear-only York solve was not detected at quartic order")
    return {"missing_matter_clock_tadpole": str(clock_error),
            "missing_mixed_boundary_quadratic_error": str(boundary_error),
            "matter_free_constraint_nonzero_residual": list(map(str, vacuum_errors)),
            "linear_only_York_quartic_error": str(truncation_error)}
