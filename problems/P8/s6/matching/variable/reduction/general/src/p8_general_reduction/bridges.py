"""Continuous specialization to the frozen source-preserving reduction."""
from functools import cache

import sympy as sp
from p8_variable_reduction import dictionary as old_dictionary
from p8_variable_reduction import stationary as old_stationary

from . import basis, general


@cache
def checks():
    now, old = general.potential(), old_stationary.potential_hessian()
    rule = {general.RATIO: old_stationary.R, general.P: old_stationary.BETA1,
            general.Z: old_stationary.BETA1}
    rule.update(dict(zip(now["entries"], old["entries"], strict=True)))
    out = {"beta1_full_quadratic": sp.expand(now["quadratic"].subs(rule)-old["quadratic"]),
           "beta1_full_cubic": sp.expand(now["cubic"].subs(rule)-old["cubic"]),
           "beta1_kappa": sp.factor(general.stationary()["kappa"].subs({
               general.MF2: old_stationary.M**2, general.RATIO: old_stationary.R,
               general.P: old_stationary.BETA1})-old["kappa"])}
    prior = old_dictionary.coefficients()
    symbols = prior["symbols"]
    phi, r, k, s = (symbols[key] for key in ("phi", "r", "kappa", "q"))
    leading = general.leading_action()
    mapping = {
        general.MG2: symbols["M_squared"], general.MF2: symbols["M_squared"],
        general.RATIO: r, general.S: s, general.SPHI: sp.diff(s, phi),
        general.SPHIPHI: sp.diff(s, phi, 2), general.KP: sp.diff(k, phi),
        general.KPP: sp.diff(k, phi, 2), general.KPPP: sp.diff(k, phi, 3),
        leading["r_phi"]: sp.diff(r, phi), leading["X"]: symbols["X"],
        general.BETA[0]: symbols["beta0"], general.BETA[1]: symbols["beta1"],
        general.BETA[2]: 0, general.BETA[3]: 0,
    }
    for key in ("F2", "F2_X", "A1", "A2", "A3", "A4", "A5", "F", "K"):
        out["full_variable_IBP_"+key] = sp.expand(leading[key].subs(mapping, simultaneous=True)-prior[key])
    jets = old_stationary.center_jets()
    center = basis.actual_center_scaling()
    c, _ = center["symbols"]
    old_c = old_stationary.C
    out["actual_center_a_is_twice_own_f_inverse"] = sp.factor(
        center["relative_correction_a"].subs(c, old_c)-2*jets["own_f_inverse_parameter_bar_u0_u2"][0])
    out["actual_center_a_second_jet"] = sp.factor(
        center["relative_correction_a_theta2"].subs(c, old_c)-2*jets["own_f_inverse_parameter_bar_u0_u2"][1])
    out["actual_center_log_ratio_second_jet"] = sp.factor(
        center["s_theta"].subs(c, old_c)-jets["r_squared_u0_u2_u4"][1]/8)
    return out
