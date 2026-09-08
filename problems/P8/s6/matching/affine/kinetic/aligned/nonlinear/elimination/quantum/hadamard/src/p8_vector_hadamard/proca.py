"""Constrained Minkowski seed and polynomial-growth reconstruction controls."""
from functools import cache

import sympy as sp


@cache
def minkowski_seed():
    k, m = sp.symbols("k m", positive=True)
    omega = sp.sqrt(k**2+m**2)
    metric = sp.diag(-1, 1, 1, 1)
    momentum = sp.Matrix([-omega, 0, 0, k])
    transverse1, transverse2 = sp.Matrix([0, 1, 0, 0]), sp.Matrix([0, 0, 1, 0])
    longitudinal = sp.Matrix([-k/m, 0, 0, omega/m])
    polarizations = (transverse1, transverse2, longitudinal)
    gram = sp.Matrix([[sp.simplify((first.T*metric*second)[0]) for second in polarizations]
                      for first in polarizations])
    projector = sum((polarization*polarization.T for polarization in polarizations), sp.zeros(4))
    return {"three_positive_polarizations_on_the_constraint_surface": gram-sp.eye(3),
            "Minkowski_Proca_polarization_sum": (projector-metric-momentum*momentum.T/m**2).applyfunc(sp.simplify),
            "Lorenz_constraint_for_each_polarization": sp.Matrix([(momentum.T*metric*p)[0] for p in polarizations]).applyfunc(sp.simplify),
            "positive_mass_shell": sp.simplify((momentum.T*metric*momentum)[0]+m**2)}


def reconstruction_checks():
    a, m, k, v, vp, H = sp.symbols("a m k v v_prime H", positive=True)
    q = k**2/a**2
    omega = sp.sqrt(m**2+q)
    z = q/omega**2
    g = sp.sqrt(a**3*m**2*q/omega**2)
    d = H*(sp.Rational(1, 2)+z)
    spatial = k*v/g
    temporal = q*(vp-d*v)/(omega**2*g)
    return {"physical_longitudinal_spatial_multiplier": sp.simplify(spatial-omega*v/(sp.sqrt(a)*m)),
            "physical_longitudinal_temporal_multiplier": sp.simplify(temporal-k*(vp-d*v)/(a**sp.Rational(5, 2)*m*omega))}
