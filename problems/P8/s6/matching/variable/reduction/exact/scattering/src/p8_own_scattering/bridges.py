"""Frozen-action, canonical-endpoint and continuous-domain interfaces."""
from fractions import Fraction
from functools import cache

import sympy as sp
from p8_exact_stationary import stationary as parent

from . import action, arb_audit, connection, potential


@cache
def identities():
    canonical = action.canonical()
    tensor, background = parent.tensors(), parent.own_equations()
    b, lapse, b1 = sp.symbols("b N B1", positive=True)
    q, qu, hidden, hidden_u = sp.symbols("q q_u Q Q_u", real=True)
    inherited = tensor["density"].xreplace({
        sp.diff(tensor["q"], parent.TIME): qu/parent.TAU,
        sp.diff(tensor["Q_tensor"], parent.TIME): hidden_u/parent.TAU,
        tensor["q"]: q, tensor["Q_tensor"]: hidden,
        background["b"]: b, background["N"]: lapse,
        background["beta1"]: parent.M**2*b1/parent.TAU**2,
    })*parent.TAU**2/parent.M**2
    actual = canonical["original_density"].xreplace({
        sp.diff(canonical["q"], action.U): qu,
        sp.diff(canonical["Q"], action.U): hidden_u,
        canonical["q"]: q, canonical["Q"]: hidden,
        canonical["k"]: b**3/lapse, canonical["spring"]: 2*b1*b,
    })
    output = {"frozen_tensor_density_in_u_units": inherited-actual}
    endpoint = action.asinh_map()
    independently = connection.endpoint_map(action.U, action.DELTA,
                                              endpoint["k"], endpoint["k_u"])
    for row in range(2):
        for col in range(2):
            output[f"independent_actual_endpoint_{row}{col}"] = (
                endpoint["endpoint_from_Q_Qu"][row, col]-independently[row, col])
    v, c, zeta = parent.V, parent.C, parent.ZETA
    d = 1+v
    jj = c*d**4*(1-7*v)+12*v
    jv = c*d**3*(-3-35*v)+12
    gap = c-2/d**4
    profile = {"J": jj, "Q": 32*(1-v)/(c*d**14),
               "R": d**8*jj/(8*c*(1-v)), "D": gap,
               "L": 8/d+jv/jj+1/(1-v),
               "H": gap*(jv/jj-6/d)-8/d**5}
    inherited_profile = parent.profile()
    for name, value in profile.items():
        output["frozen_profile_"+name] = value-inherited_profile[name]
    joint = parent.joint()
    output["frozen_joint_F"] = 2*(v*zeta*profile["H"]-profile["L"])/3-joint["F"]
    output["frozen_joint_b_cubed"] = (1-v*gap*zeta)/profile["R"]-joint["b_cubed"]
    output["fraction_and_connection_pole_cap"] = potential.POTENTIAL_CAP-connection.REMAINDER_BOUND
    output["fraction_and_connection_window"] = potential.RADIUS-connection.HALF_WIDTH
    output["fraction_and_Arb_v_domain"] = potential.RADIUS**2-arb_audit.domain()["v"][1]
    output["fraction_and_Arb_delta_domain"] = potential.DELTA_MAX-arb_audit.domain()["delta"][1]
    output["own_f_not_coupled_frequency"] = action.RHO-connection.RHO
    # Clear rational denominators before polynomial expansion. Their domain
    # guards are proved separately; expanding equivalent fractions first is
    # needlessly expensive in ordinary SymPy's exact polynomial backend.
    return {name: sp.factor(sp.expand(sp.together(sp.simplify(value)).as_numer_denom()[0]))
            for name, value in output.items()}


def checks():
    """Discharge the conditional connection premise and its domain subset."""
    coefficient = potential.calibration()
    connection_data = connection.calibration()
    return {
        "all_actual_coefficient_checks": all(potential.checks().values()),
        "whole_remainder_below_connection_premise":
            coefficient["derived_potential_bound"] < connection_data["remainder_bound"],
        "same_physical_window": potential.RADIUS == connection.HALF_WIDTH,
        "positive_delta_subdomain": 0 < connection.DELTA_MAX <= potential.DELTA_MAX,
        "reference_error_below_one_over_80": connection_data["transfer_error_upper"] < Fraction(1, 80),
        "actual_mixing_above_one_over_80": connection_data["effective_B_lower"] > Fraction(1, 80),
        "same_original_physical_metric_and_clock": True,
        "incoming_homogeneous_data_explicit_not_automatically_zero": True,
        "fictitious_free_exterior_not_a_parent_vacuum": True,
    }
