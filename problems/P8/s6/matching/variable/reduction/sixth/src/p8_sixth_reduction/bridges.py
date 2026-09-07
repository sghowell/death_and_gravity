"""Independent literal derivation and continuous frozen-profile bridges."""
from functools import cache

import sympy as sp
from p8_variable_reduction import stationary as prior

from . import literal, tensor


def literal_to_primary():
    return {literal.A: tensor.INV[0], literal.AP: tensor.INV[1],
            literal.A2: tensor.INV[2], literal.A3: tensor.INV[3], literal.A4: tensor.INV[4],
            literal.RV: sp.sqrt(tensor.R2), literal.L: tensor.LOG[0],
            literal.L1: tensor.LOG[1], literal.L2: tensor.LOG[2], literal.L3: tensor.LOG[3],
            literal.Q1: tensor.Q[1], literal.Q2: tensor.Q[2], literal.Q3: tensor.Q[3]}


@cache
def checks():
    primary, raw = tensor.derive(), literal.diagonal_action()
    normal, point = literal.normal_and_euler(), literal.independent_center()
    mapping = literal_to_primary()
    out = {}
    for index, value in enumerate(raw["h"]):
        out[f"coordinate_Schouten_component_{index}"] = sp.factor(
            value.subs(mapping)-primary["h_mixed_diagonal"][index])
    for (i, j), expected in primary["raw_coefficients"].items():
        out[f"literal_ADM_Newton_raw_a{i}{j}"] = sp.factor(raw["raw_coefficients"][f"a{i}{j}"].subs(mapping)-expected)
    for n in (1, 2, 3):
        out[f"independent_normal_B{n}"] = sp.factor(normal[f"B{n}"].subs(mapping)-primary[f"B{n}"])
    for n, value in normal["euler_coefficients"].items():
        out[f"independent_full_Euler_{n}"] = sp.factor(value.subs(mapping)-primary["euler_coefficients"][n])
    center = tensor.center()
    for n, value in point["B"].items():
        out[f"independent_series_center_B{n}"] = sp.factor(value.subs(point["c"], tensor.C)-center["B_dimensionless"][n])
    for n, value in point["E"].items():
        out[f"independent_series_center_E{n}"] = sp.factor(value.subs(point["c"], tensor.C)-center["sixth_E_dimensionless"][n])
    old, current = prior.profile(), tensor.center_jets()
    for key in ("b1", "b4", "r_cubed", "kappa_bar", "log_r_u_over_u"):
        out["continuous_frozen_"+key] = sp.factor(current[key].subs({current["v"]: old["v"], tensor.C: old["c"]})-old[key])
    out["independent_rational_root_representation"] = sp.factor(
        point["r_cubed"].subs({point["v"]: old["v"], point["c"]: old["c"]})-old["r_cubed"])
    jets = prior.center_jets()
    now = current["substitution"]
    out["own_inverse_center_factor_two"] = sp.factor(now[tensor.INV[0]]-2*jets["own_f_inverse_parameter_bar_u0_u2"][0])
    out["own_inverse_second_jet_factor_two"] = sp.factor(now[tensor.INV[2]]-2*jets["own_f_inverse_parameter_bar_u0_u2"][1])
    out["stationary_ratio_log_second_jet"] = sp.factor(now[tensor.LOG[1]]-jets["r_squared_u0_u2_u4"][1]/8)
    constant = tensor.calibration()["sixth_highest_symbol"].subs({
        tensor.R2: prior.R**2, tensor.INV[0]: prior.M**2*prior.R/prior.BETA1})
    out["constant_ratio_source_preserving_high_symbol"] = sp.factor(
        constant-prior.constant_ratio_calibration()["six_derivative_flat_curvature_coefficient"])
    return out
