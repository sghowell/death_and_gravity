"""Literal lower-scalar deformation; no frozen coefficient is overwritten."""
from functools import cache

import sympy as sp
from p8_affine_kinetic import scalar as old

u = old.u
epsilon = sp.Symbol("epsilon_margin", positive=True)
x = sp.Symbol("source_x", real=True)
N = sp.Symbol("physical_lapse", positive=True)


@cache
def deformation():
    h = old.background()["h"]
    value = epsilon*(x+1)**2/h**2
    U = ((h-1+N**-2)/h)**(-sp.Rational(3, 4))
    density = N*U*value.subs(x, -N**-2)
    delta_J = 4*epsilon/h**2
    return {"physical_lower_scalar": value, "hat_volume_ratio": U,
            "hat_density": density, "delta_J": delta_J,
            "clock_value": value.subs(x, -1),
            "clock_x_first": sp.diff(value, x).subs(x, -1),
            "clock_phi_first": sp.diff(value, u).subs(x, -1),
            "clock_x_first_total_time_derivative": sp.diff(sp.diff(value, x).subs(x, -1), u),
            "clock_x_second": sp.factor(sp.diff(value, x, 2).subs(x, -1)),
            "density_value": sp.simplify(density.subs(N, 1)),
            "density_lapse_first": sp.simplify(sp.diff(density, N).subs(N, 1)),
            "density_lapse_quadratic": sp.simplify(sp.diff(density, N, 2).subs(N, 1)/2),
            "Hamiltonian_joint_auxiliary_change": sp.diag(-2*delta_J, 0)}


@cache
def checks():
    d = deformation()
    return {name: d[name] for name in ("clock_value", "clock_x_first", "clock_phi_first",
            "clock_x_first_total_time_derivative", "density_value", "density_lapse_first")} | {
            "literal_physical_volume_quadratic_lapse_change": sp.factor(d["density_lapse_quadratic"]-d["delta_J"]),
            "lower_scalar_second_x_derivative": sp.factor(d["clock_x_second"]-2*epsilon/old.background()["h"]**2),
            "joint_auxiliary_change_from_full_Hamiltonian_density":
                sp.diag(sp.simplify(-sp.diff(d["hat_density"], N, 2).subs(N, 1)), 0)-d["Hamiltonian_joint_auxiliary_change"]}
