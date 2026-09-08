"""Actual per-mode energy/pressure contact matrices before momentum integration."""
from functools import cache

import sympy as sp
from p8_aligned_quantum import kernel

from . import modes

N, a, q, m2 = modes.N, modes.scale, modes.q, modes.mass2
alpha2, beta2 = sp.symbols("alpha_lapse_second beta_lapse_second", real=True)


def Da(value):
    return a*sp.diff(value, a)-2*q*sp.diff(value, q)


def clean(matrix):
    return sp.ImmutableMatrix(matrix.applyfunc(sp.factor))


@cache
def matrices():
    aN = 1+modes.alpha*(N-1)+alpha2*(N-1)**2/2
    bN = 1+modes.beta*(N-1)+beta2*(N-1)**2/2
    out = {}
    z, omega2 = q/(q+m2), q+m2
    targets = {
        "transverse": (sp.diag(q+m2*(1+modes.beta), 1)/a**3,
                       sp.diag((q-m2)/3, sp.Rational(1, 3))/a**3),
        "longitudinal": (sp.diag(omega2*(1+modes.beta), 1-modes.alpha*z)/a**3,
                         sp.diag(-omega2/3, (1+2*z)/3)/a**3)}
    for kind, data in ((name, modes.canonical()[name]) for name in ("transverse", "longitudinal")):
        g2 = data["g_squared"].subs({modes.am: aN, modes.bm: bN}, simultaneous=True)
        bare = data["bare_frequency_squared"].subs({modes.am: aN, modes.bm: bN}, simultaneous=True)
        Hamiltonian = sp.diag(g2*bare, 1/g2)
        rho = clean(Hamiltonian.diff(N)/a**3)
        pressure = clean(-Hamiltonian.applyfunc(Da)/(3*N*a**3))
        g0 = sp.sqrt(g2.subs(N, 1))
        inverse_map = sp.diag(1/g0, g0)
        energy_clock = clean(inverse_map.T*rho.subs(N, 1)*inverse_map)
        pressure_clock = clean(inverse_map.T*pressure.subs(N, 1)*inverse_map)
        contacts = {name: clean((modes.n*matrix.diff(N)+modes.zeta*matrix.applyfunc(Da)).subs(N, 1))
                    for name, matrix in (("energy", rho), ("pressure", pressure))}
        out[kind] = {"Hamiltonian_quadratic_matrix": clean(Hamiltonian),
                     "energy_matrix": rho.subs(N, 1), "pressure_matrix": pressure.subs(N, 1),
                     "energy_contact_matrix": contacts["energy"], "pressure_contact_matrix": contacts["pressure"],
                     "energy_replays_actual_frozen_readout": clean(energy_clock-targets[kind][0]),
                     "pressure_replays_actual_frozen_readout": clean(pressure_clock-targets[kind][1])}
    return out


@cache
def actual_mass_jets():
    p = kernel.geometry.P
    h = sp.Symbol("h", positive=True)
    pN, pNN = -1/(2*h), 3/(2*h)-1/(2*h**2)
    expected = {"a": (4/(9*h), -4/(3*h)-8/(9*h**2)),
                "b": (28/(81*h), -28/(27*h)+sp.Rational(152, 729)/h**2)}
    out = {}
    for kind, mass in kernel.masses().items():
        first = sp.factor(sp.diff(mass, p).subs(p, sp.Rational(1, 2))*pN)
        second = sp.factor(sp.diff(mass, p, 2).subs(p, sp.Rational(1, 2))*pN**2
                           +sp.diff(mass, p).subs(p, sp.Rational(1, 2))*pNN)
        out[kind] = {"N_first": first, "N_second": second,
                     "first_jet_identity": sp.factor(first-expected[kind][0]),
                     "second_jet_identity": sp.factor(second-expected[kind][1])}
    return out


@cache
def covariance_identity():
    # This is the symmetric equal-time covariance of a canonical real mode.
    x, y, z, dx, dy, dz, aa, bb, da, db = sp.symbols(
        "xx xp pp delta_xx delta_xp delta_pp A B delta_A delta_B", real=True)
    covariance = sp.Matrix([[x, y], [y, z]])
    changed = sp.Matrix([[dx, dy], [dy, dz]])
    generator = sp.Matrix([[0, aa], [-bb, 0]])
    delta_generator = sp.Matrix([[0, da], [-db, 0]])
    t = sp.Symbol("variation_parameter")
    literal = ((generator+t*delta_generator)*(covariance+t*changed)
               +(covariance+t*changed)*(generator+t*delta_generator).T).diff(t).subs(t, 0)
    target = generator*changed+changed*generator.T+delta_generator*covariance+covariance*delta_generator.T
    kx, ky, dkx, dky = sp.symbols("readout_xx readout_pp contact_xx contact_pp", real=True)
    readout, contact = sp.diag(kx, ky), sp.diag(dkx, dky)
    observable = sp.trace((readout+t*contact)*(covariance+t*changed))/2
    return {"covariance_retarded_variation_equation": clean(literal-target),
            "observable_contact_and_state_variation": sp.factor(
                sp.diff(observable, t).subs(t, 0)-sp.trace(contact*covariance+readout*changed)/2)}


@cache
def checks():
    out = covariance_identity()
    for kind, data in matrices().items():
        for name in ("energy_replays_actual_frozen_readout", "pressure_replays_actual_frozen_readout"):
            out[kind+"_"+name] = data[name]
    for kind, data in actual_mass_jets().items():
        for name in ("first_jet_identity", "second_jet_identity"):
            out[kind+"_"+name] = data[name]
    return out
