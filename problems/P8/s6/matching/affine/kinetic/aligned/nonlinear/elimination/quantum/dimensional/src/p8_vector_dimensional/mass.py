"""Covariant first-mass-variation pole reconstruction on physical FLRW."""
from functools import cache

import sympy as sp
from p8_aligned_quantum import potential

from . import local, metric


@cache
def reconstruction():
    H, Hd, alpha, beta = local.H, local.Hd, local.alpha, local.beta
    R = 6*(Hd+2*H**2)
    boxR = sp.factor(-local.derivative(local.derivative(R))-3*H*local.derivative(R))
    scalar_energy = metric.variations()["minimal_scalar_curvature_pole_energy"]
    scalar_pressure = sp.factor(-scalar_energy-local.derivative(scalar_energy)/(3*H))
    Znormal = scalar_energy-boxR/6
    Ztrace = -scalar_energy+3*scalar_pressure+sp.Rational(2, 3)*boxR
    mass_trace = alpha+3*beta
    ricci_contraction = 3*alpha*(Hd+H**2)+3*beta*(Hd+3*H**2)
    curvature_mass_pole = sp.factor(2*(R*mass_trace/12-sp.Rational(5, 6)*ricci_contraction))
    fourth_mass_pole = sp.factor((beta*Ztrace+(beta-alpha)*Znormal)/2)
    ordinary = {alpha: 0, beta: 0}
    data = local.integrated()
    return {"R": R, "boxR": boxR,
            "minimal_scalar_pole_energy": scalar_energy, "minimal_scalar_pole_pressure": scalar_pressure,
            "scalar_gradient_normal_pole": Znormal, "scalar_gradient_trace_pole": sp.factor(Ztrace),
            "covariant_second_order_mass_pole": curvature_mass_pole,
            "covariant_fourth_order_mass_pole": fourth_mass_pole,
            "scalar_pole_EOM_trace_identity": sp.factor(Ztrace-boxR/3),
            "actual_second_order_mass_pole_matches_covariant_coincidence": sp.factor(
                data[1]["pole"]-data[1]["pole"].subs(ordinary)-curvature_mass_pole),
            "actual_fourth_order_mass_pole_matches_scalar_gradient_reconstruction": sp.factor(
                data[2]["pole"]-data[2]["pole"].subs(ordinary)-fourth_mass_pole),
            "actual_fourth_order_mass_pole_compact_form": sp.factor(
                fourth_mass_pole-(beta-alpha)*scalar_energy/2-(alpha+beta)*boxR/12)}


@cache
def frozen_flat_finite_match():
    jets = potential.clock_jets()
    h = jets["h"]
    expected = jets["finite_weight"]["value"]+jets["finite_weight"]["N_first"]
    actual = local.integrated()[0]["radial_finite_part"].subs({local.alpha: 4/(9*h), local.beta: 28/(81*h), local.ell: 0})
    return {"actual_finite_flat_lapse_term_preserves_frozen_MSbar_jets": sp.factor(actual-expected)}


@cache
def checks():
    data = reconstruction()
    return {name: data[name] for name in
            ("scalar_pole_EOM_trace_identity", "actual_second_order_mass_pole_matches_covariant_coincidence",
             "actual_fourth_order_mass_pole_matches_scalar_gradient_reconstruction", "actual_fourth_order_mass_pole_compact_form")}
