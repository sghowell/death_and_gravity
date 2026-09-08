"""Dimensional lapse variation and counterterm Ward-identity controls."""
from functools import cache

import sympy as sp

from . import local


@cache
def variations():
    D, H, Hd = local.dimension, local.H, local.Hd
    N, Nd = sp.symbols("N N_dot", real=True)
    acc = (Hd+H**2)/N**2-H*Nd/N**3
    expansion2 = H**2/N**2
    R = 2*D*acc+D*(D-1)*expansion2
    Ricci2 = D**2*acc**2+D*(acc+(D-1)*expansion2)**2
    Riemann2 = 4*D*acc**2+2*D*(D-1)*expansion2**2
    a4_at_four = -R**2/8+sp.Rational(29, 60)*Ricci2-Riemann2/15

    def rho(invariant):
        density = N*invariant
        first = sp.diff(density, N).subs({N: 1, Nd: 0})
        rate = sp.diff(density, Nd).subs({N: 1, Nd: 0})
        return sp.factor(-first+local.derivative(rate)+D*H*rate)

    densities = {0: 3*rho(sp.Integer(1)),
                 1: rho(R), 2: 2*rho(a4_at_four)}
    # Relative to m^4/(64*pi²*epsilon), the Lorentzian pole density
    # is +3 and its lapse energy is -3. The derivative orders have
    # their own m² and m^0 prefactors, kept outside this local table.
    pressure = {order: sp.factor(-value-local.derivative(value)/(D*H)) for order, value in densities.items()}
    # Independent local heat action input: a0=D, a2=(D/6-1)R,
    # and the general one-form-minus-scalar a4 traces. Differentiate
    # the bundle-rank coefficient, not the metric dimension, here.
    a4_rank_derivative = R**2/72-Ricci2/180+Riemann2/180
    finite_heat_densities = {0: -(3*local.ell-sp.Rational(5, 2)),
                            1: -(local.ell-sp.Rational(5, 3))*R,
                            2: -2*local.ell*a4_at_four-4*a4_rank_derivative}
    heat_energy = {order: sp.factor(rho(value).subs(D, 3)) for order, value in finite_heat_densities.items()}
    return {"R": R, "Ricci2": Ricci2, "Riemann2": Riemann2,
            "rho_pole_counterterm_in_D": densities,
            "pressure_pole_counterterm_in_D": pressure,
            "energy_evanescent_counterterm": {order: sp.factor(-2*sp.diff(value, D).subs(D, 3)) for order, value in densities.items()},
            "pressure_evanescent_counterterm": {order: sp.factor(-2*sp.diff(value, D).subs(D, 3)) for order, value in pressure.items()},
            "independent_finite_heat_action_lapse_energy": heat_energy,
            "minimal_scalar_curvature_pole_energy": sp.factor(2*rho(a4_rank_derivative).subs(D, 3))}


@cache
def ordinary_matching():
    metric = variations()
    raw_rho, raw_p = local.integrated(), local.integrated_pressure()
    ordinary = {local.alpha: 0, local.beta: 0}
    out, matched = {}, {}
    for order in range(3):
        rho_pole, pressure_pole = [table[order]["pole"].subs(ordinary) for table in (raw_rho, raw_p)]
        rho_finite, pressure_finite = [table[order]["radial_finite_part"].subs(ordinary) for table in (raw_rho, raw_p)]
        out[str(order)+"_actual_lapse_variation_pole"] = sp.factor(rho_pole-metric["rho_pole_counterterm_in_D"][order].subs(local.dimension, 3))
        out[str(order)+"_covariant_pressure_pole"] = sp.factor(pressure_pole-metric["pressure_pole_counterterm_in_D"][order].subs(local.dimension, 3))
        ward_raw = sp.factor(local.derivative(rho_finite)+3*local.H*(rho_finite+pressure_finite))
        out[str(order)+"_radial_finite_Ward_defect_identity"] = sp.factor(ward_raw-2*local.H*(rho_pole+pressure_pole))
        matched_rho = sp.factor(rho_finite-metric["energy_evanescent_counterterm"][order])
        matched_p = sp.factor(pressure_finite-metric["pressure_evanescent_counterterm"][order])
        out[str(order)+"_counterterm_varied_before_limit_restores_Ward_identity"] = sp.factor(
            local.derivative(matched_rho)+3*local.H*(matched_rho+matched_p))
        out[str(order)+"_independent_finite_heat_action_energy"] = sp.factor(
            matched_rho-metric["independent_finite_heat_action_lapse_energy"][order])
        matched[order] = {"energy": matched_rho, "pressure": matched_p, "radial_Ward_defect": ward_raw}
    return {"identities": out, "ordinary_covariantly_subtracted_local_terms": matched}
