"""Entire current parent clock germs, primitive and physical ADM coefficients."""

from functools import cache

import sympy as s
from p8_affine_vacuum_domain import family as original
from p8_vacuum_affine_heavy_curved_state import quantum, state
from p8_vacuum_affine_heavy_scalar_parent import family as heavy
from p8_vacuum_affine_quantum_retuning import profile

u, X = heavy.u, heavy.X
N = s.Symbol("positive_background_lapse", positive=True)
h = (1 + u * u) ** 3
H = 4 * u / (1 + u * u)
ell = 1 / (10 * (1 + u * u) ** 6)
FAMILIES = ("M", "U", "C3", "Cchi", "B", "Fhat")
RHO = s.Function("fixed_total_normalized_reference_energy", real=True)(u)
PRESSURE = s.Function("fixed_total_normalized_reference_pressure", real=True)(u)


@cache
def fixed_functions():
    old, full, current = original.data(), heavy.coefficients(), quantum.fixed_profile()
    heavy_map = {state.TIME: u, state.X: X}
    rho = profile.rho + current["rho"].subs(state.TIME, u) / state.KAPPA
    pressure = profile.P + current["pressure"].subs(state.TIME, u) / state.KAPPA
    A, B = -pressure, -(rho + pressure) / 2
    tree = original.original.data()["original_retuned_tree_scalar"]
    return {
        "R_full": full["R"],
        "F_full": full["F"] + current["DeltaF"].subs(heavy_map, simultaneous=True),
        "R_clock": 1 + (X - 1) / h,
        "F_clock": tree + A + B * (X - 1),
        "rho": rho,
        "pressure": pressure,
        "A": A,
        "B": B,
        "vacuum_pressure_total": profile.PV
        + (current["p_v_H"] + current["p_v_Phi"]) / state.KAPPA,
        "original": old,
        "heavy": full,
        "tree": tree,
    }


@cache
def source_germs():
    d = fixed_functions()
    old, full = d["original"], d["heavy"]
    remainder_F = s.exp(-(u**4)) * (1 - old["bump"]) * (
        original.previous.vacuum_lower_function() - d["tree"]
    ) - (d["A"] + d["B"] * (X - 1) + d["vacuum_pressure_total"])
    remainder_R = -(X - 1) * (1 - old["bump"]) / h
    denominator = X**original.N + (1 - X) ** original.N
    marker = s.Symbol("complete_switch_inverse_denominator")
    checks = {}
    for name, local in (("R", remainder_R), ("F", remainder_F)):
        expression = (
            d[name + "_full"]
            - d[name + "_clock"]
            - full["delta_" + name]
            - (1 - old["T"]) * local
        )
        lifted = expression.xreplace({1 / denominator: marker})
        if lifted.xreplace({marker: 1 / denominator}) != expression:
            raise ValueError(
                "The exact common-denominator replacement did not roundtrip"
            )
        checks["entire_" + name + "_clock_difference_factorization"] = s.expand(lifted)
    checks["literal_complete_complement_switch_identity"] = s.cancel(
        (1 - old["T"]) * denominator - (1 - X) ** original.N
    )
    checks["complete_heavy_R_eighth_factor"] = s.cancel(
        full["delta_R"] / heavy.switch() - heavy.GAMMA * heavy.K0 * X**2
    )
    checks["complete_heavy_F_eighth_factor"] = s.cancel(
        full["delta_F"] / heavy.switch()
        + heavy.LAMBDA * heavy.K0 * X**2
        - (heavy.GAMMA / 3 + heavy.CONTACT / 24) * heavy.K0 * u**4
    )
    checks["entire_fixed_profile_clock_name_binding"] = s.expand(
        d["F_clock"]
        - (d["tree"] - PRESSURE - (RHO + PRESSURE) * (X - 1) / 2).subs(
            {RHO: d["rho"], PRESSURE: d["pressure"]}, simultaneous=True
        )
    )
    return {
        "entire_R_difference": (1 - old["T"]) * remainder_R + full["delta_R"],
        "entire_F_difference": (1 - old["T"]) * remainder_F + full["delta_F"],
        "fixed_total_normalized_energy": d["rho"],
        "fixed_total_normalized_pressure": d["pressure"],
        "all_three_vacuum_constants": d["vacuum_pressure_total"],
        "fixed_total_coefficient_function_bindings": {
            RHO: d["rho"],
            PRESSURE: d["pressure"],
        },
        "clock_jet_rule": "Both entire differences have (X-1)^8 factors, with a smooth nonvanishing switch denominator at X1. Fixed smooth time derivatives preserve that factor. This licenses the specified finite jets, not replacement of the full off-tube functions.",
        "checks": checks,
        "gates": {
            "full_original_switch_order_1024": original.N == 1024,
            "full_heavy_clock_order_eight": heavy.ORDER == 8,
            "complete_switch_denominator_one_on_clock": denominator.subs(X, 1) == 1,
            "fixed_profiles_not_varied_as_live_quantum_means": True,
        },
    }


@cache
def lapse_jets():
    # These names bind to the entire two fixed profiles in source_germs().
    # They are never live means varied with the probe or new tunable functions.
    R = 1 + (X - 1) / h
    F = (
        original.original.data()["original_retuned_tree_scalar"]
        - PRESSURE
        - (RHO + PRESSURE) * (X - 1) / 2
    )
    Ru, RX = s.diff(R, u), s.diff(R, X)
    U = R ** -s.Rational(3, 4)
    primitive_s = U * 3 * X * Ru * RX / (2 * R)
    first_lapse = (-primitive_s / N**2).subs(X, N**-2)
    Ijets = (s.Integer(0),) + tuple(
        s.factor(s.diff(first_lapse, N, j - 1).subs(N, 1)) for j in range(1, 5)
    )
    primitive = sum(Ijets[j] * (N - 1) ** j / s.factorial(j) for j in range(5))
    functions = {
        "M": R ** s.Rational(1, 4),
        "U": U,
        "C3": R ** s.Rational(3, 4) / 2,
        "Cchi": R ** -s.Rational(1, 4),
        "B": -U * s.sqrt(X) * Ru / 2,
        "Fhat": U * (F + 9 * X * Ru**2 / (16 * R)),
    }
    functions = {key: value.subs(X, N**-2) for key, value in functions.items()}
    functions["B"] -= primitive
    functions["Fhat"] -= s.diff(primitive, u) / N
    rows = {
        key: tuple(s.factor(s.diff(value, N, j).subs(N, 1)) for j in range(5))
        for key, value in functions.items()
    }
    return {
        "rows": rows,
        "primitive_lapse_jets_zero_through_four": Ijets,
        "R_lapse_jets_zero_through_four": tuple(
            s.factor(s.diff(R.subs(X, N**-2), N, j).subs(N, 1)) for j in range(5)
        ),
    }


@cache
def reference_data():
    packet = lapse_jets()
    rows = packet["rows"]

    def divide_N(values):
        return (
            values[0],
            values[1] - values[0],
            values[2] - 2 * values[1] + 2 * values[0],
        )

    D, DN, DNN = divide_N(rows["M"])
    Z, ZN, ZNN = divide_N(rows["U"])
    b0, b1, b2 = rows["B"][:3]
    f0, f1, f2 = rows["Fhat"][:3]
    L0 = -3 * D * H**2 + 3 * b0 * H + f0 + Z * ell**2 / 2
    theta, w = -H * DN + b1 / 2, ell * ZN
    Cnn = (-3 * H**2 * DNN + 3 * H * b2 + 2 * f1 + f2 + ZNN * ell**2 / 2) / 2
    Qn = -3 * H**2 * DN + 3 * H * b1 + f0 + f1 + ZN * ell**2 / 2
    cv = -18 * D * H + 9 * b0
    Vvv = s.Rational(9, 2) * L0 - (s.diff(cv, u) + 3 * H * cv) / 2
    Vvs = -3 * (s.diff(Z * ell, u) + 3 * H * Z * ell)
    J = Cnn + 3 * theta**2 / D - w**2 / (2 * Z)
    delta = 1 / (2 * h)
    E = 1 - 3 * delta
    theta0 = H - u / (1 + u * u) ** 4
    gradient = (
        theta0 * s.diff(E, u)
        - E * s.diff(theta0, u)
        + H * E * theta0
        - theta0**2
        - ell**2 * E**2 / 2
    )
    A, Bprofile = -PRESSURE, -(RHO + PRESSURE) / 2
    expected_J = (
        gradient
        + 1 / (50 * (1 + u * u) ** 6)
        + (21 * delta**2 - 3 * delta) * A / 2
        + (1 - 6 * delta) * Bprofile
    )
    checks = {
        "complete_reference_D": D - 1,
        "complete_reference_Z": Z - 1,
        "complete_reference_Theta": theta - theta0,
        "complete_reference_mixed_momentum": w + ell * E,
        "complete_reference_current_lapse_pivot": J - expected_J,
        "complete_reference_lapse_volume_current": Qn - RHO + 3 * delta * PRESSURE,
        "complete_reference_metric_potential": Vvv - s.Rational(9, 2) * A,
        "complete_reference_matter_Ward_contact": Vvs,
        "complete_reference_spatial_curvature": 2 * rows["C3"][0] - 1,
        "complete_reference_lapse_spatial_curvature": 4
        * (rows["C3"][0] + rows["C3"][1])
        - 2 * E,
        "complete_reference_matter_gradient": rows["Cchi"][0] - 1,
        "primitive_first_lapse_jet": packet["primitive_lapse_jets_zero_through_four"][
            1
        ],
        "primitive_second_lapse_jet": packet["primitive_lapse_jets_zero_through_four"][
            2
        ]
        + 3 * s.diff(h, u) / h**3,
    }
    return {
        "all_six_coefficient_lapse_jets_zero_through_four": rows,
        "primitive_lapse_jets_zero_through_four": packet[
            "primitive_lapse_jets_zero_through_four"
        ],
        "whole_reference_coefficients": {
            "Jc": expected_J,
            "Theta": theta0,
            "ell": ell,
            "E": E,
            "A": A,
            "T": RHO - 3 * delta * PRESSURE,
        },
        "checks": {name: s.factor(value) for name, value in checks.items()},
        "gates": {
            "matter_gradient_not_mistaken_for_volume_factor": rows["Cchi"][1]
            != rows["U"][1],
            "all_required_higher_primitive_jets_retained": all(
                packet["primitive_lapse_jets_zero_through_four"][j] != 0
                for j in (2, 3, 4)
            ),
            "exactly_six_complete_families_with_five_lapse_jets": len(rows) == 6
            and all(len(value) == 5 for value in rows.values()),
        },
    }
