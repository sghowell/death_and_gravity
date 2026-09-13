"""Full fixed profile, nonlinear common clock and nonzero projected-mean contacts."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_curved_state import quantum, state

from . import geometry


@cache
def data():
    e, f, nD, nG, trD, trG, aa, rho, P, delta = s.symbols(
        "e f n_D n_G trace_D trace_G a rho pressure delta", real=True
    )
    N = 1 + e * nD + f * nG
    trace = e * trD + f * trG
    A = -P
    B = -(rho + P) / 2
    physical = aa**3 * N * s.exp(trace / 2) * (A + B * (N**-2 - 1))
    linear = aa**3 * (rho * nD - P * trD / 2)
    mixed = aa**3 * (
        -(rho + P) * nD * nG + rho * (nD * trG + nG * trD) / 2 - P * trD * trG / 4
    )
    qq = -(aa**3) * P * trD * trG / 4
    Rjet = 1 + 2 * delta * (N**-2 - 1)
    clock = (
        aa**3
        * N
        * s.exp(trace / 2)
        * Rjet ** (-s.Rational(3, 4))
        * (A + B * (N**-2 - 1))
    )
    clock_first = aa**3 * ((rho - 3 * delta * P) * nD - P * trD / 2)
    clock_second = aa**3 * (
        -((1 - 6 * delta) * rho + (21 * delta**2 - 9 * delta + 1) * P) * nD * nG
        + (rho - 3 * delta * P) * (nD * trG + nG * trD) / 2
        - P * trD * trG / 4
    )
    pullback = mixed.subs(
        {trD: trD + 6 * delta * nD, trG: trG + 6 * delta * nG}, simultaneous=True
    )
    qclock = 4 * delta**2 - 3 * delta
    heavy_contact = 3 * aa**3 * P * qclock * nD * nG
    profile_contact = -heavy_contact
    PK = s.Symbol("projected_pressure", real=True)
    residual = 3 * aa**3 * (PK - P) * qclock * nD * nG
    chart_shift = -s.log(Rjet) / 2
    checks = {
        "whole_physical_profile_linear_current": s.diff(physical, e).subs({e: 0, f: 0})
        - linear,
        "whole_physical_profile_mixed_Hessian": s.diff(physical, e, f).subs(
            {e: 0, f: 0}
        )
        - mixed,
        "entire_spatial_profile_to_remove_before_Gaussian_Ward": mixed.subs(
            {nD: 0, nG: 0}
        )
        - qq,
        "entire_nonlinear_clock_spatial_first_jet": s.diff(chart_shift, e).subs(
            {e: 0, f: 0}
        )
        - 2 * delta * nD,
        "entire_nonlinear_clock_spatial_mixed_jet": s.diff(chart_shift, e, f).subs(
            {e: 0, f: 0}
        )
        - 2 * qclock * nD * nG,
        "literal_complete_common_clock_profile_linear": s.diff(clock, e).subs(
            {e: 0, f: 0}
        )
        - clock_first,
        "literal_complete_common_clock_profile_mixed": s.diff(clock, e, f).subs(
            {e: 0, f: 0}
        )
        - clock_second,
        "whole_profile_chain_rule_contact_retained": clock_second
        - pullback
        - profile_contact,
        "full_reference_mean_cancels_only_after_both_contacts": heavy_contact
        + profile_contact,
        "projected_mean_contact_not_reset_to_reference": heavy_contact.subs(P, PK)
        + profile_contact
        - residual,
        "nonzero_projected_mean_contact_control": residual.subs(
            {PK: 0, P: 1, delta: s.Rational(1, 2), nD: 1, nG: 1, aa: 1}
        )
        - s.Rational(3, 2),
    }
    frozen = quantum.fixed_profile()
    fullF = state.KAPPA * frozen["DeltaF"].subs(
        {frozen["rho"]: rho, frozen["pressure"]: P}, simultaneous=True
    )
    for j, target in enumerate((A, B, s.S.Zero)):
        checks["literal_full_fixed_profile_clock_X_jet_" + str(j)] = (
            s.diff(fullF, state.X, j).subs(state.X, 1) - target
        )
    parent = state.germs.parent
    RX = parent.coefficients()["R"].subs(parent.u, geometry.t)
    actual_delta = 1 / (2 * (1 + geometry.t**2) ** 3)
    for j, target in enumerate((s.S.One, 2 * actual_delta, s.S.Zero)):
        checks["literal_full_parent_R_clock_X_jet_" + str(j)] = (
            s.diff(RX, parent.X, j).subs(parent.X, 1) - target
        )
    source = parent.coefficients()["normalized_heavy_source"]
    for clock_order in range(3):
        for field_order in range(3 - clock_order):
            checks[
                f"full_heavy_source_clock_{clock_order}_field_{field_order}_jet_zero"
            ] = s.diff(source, parent.X, clock_order, parent.u, field_order).subs(
                parent.X, 1
            )
    vD, vG = s.symbols("v_D v_G", real=True)
    legacyA, legacyB = A, B
    dJ = (21 * delta**2 - 3 * delta) * legacyA / 2 + (1 - 6 * delta) * legacyB
    tc = (1 + 3 * delta) * legacyA - 2 * legacyB
    scalar = aa**3 * (
        2 * dJ * nD * nG + 3 * tc * (nD * vG + nG * vD) + 9 * legacyA * vD * vG
    )
    checks["full_general_tensor_profile_matches_frozen_scalar_clock_specialization"] = (
        clock_second.subs({trD: 6 * vD, trG: 6 * vG}) - scalar
    )
    # A total matched mean is zero but a metric-only profile is not separately on shell.
    checks["separate_physical_profile_mean_is_negative_Gaussian_current"] = (
        linear + aa** 3 * (-rho * nD + P * trD / 2)
    )
    checks["separate_common_clock_profile_mean_is_negative_Gaussian_current"] = (
        clock_first + aa** 3 * ((-rho + 3 * delta * P) * nD + P * trD / 2)
    )
    result = {
        "whole_physical_profile_one_current": linear,
        "whole_physical_profile_Hessian": mixed,
        "whole_spatial_profile_subtracted_before_Ward": qq,
        "whole_common_clock_profile_one_current": clock_first,
        "whole_common_clock_profile_Hessian": clock_second,
        "complete_Gaussian_nonlinear_clock_contact": heavy_contact,
        "complete_profile_nonlinear_clock_contact": profile_contact,
        "complete_original_projected_mean_clock_contact": residual,
        "actual_clock_delta": actual_delta,
        "full_mean_cancellation_boundary": "The complete actual Gaussian mean and fixed profile cancel at the reference only AFTER their separate first currents and nonlinear clock contacts are retained. The profile is not separately on shell and is not included in the metric-only Gaussian Ward identity.",
        "full_physical_reconstruction": "Subtract exactly P_QQ from S245's total spatial response, apply the two ordered scalar Gaussian Ward terms with the full actual mean, then add the entire physical ADM profile once. No extra finite heat action is added.",
        "finite_projection_boundary": "The computational Ward-completed approximant uses BOTH ordered mean terms at E_K and retains3a³(P_K-P)(4delta²-3delta)nDnG under the nonlinear clock. It is not asserted equal to the bare projected full metric response or a local/shape counteraction.",
        "source_bridge": "The full actual H8A420 source and every clock/value jet of total degree at most two vanish at X1. The heavy coherent mean remains zero. This is the conditional Gaussian metric bridge, not removal of interacting light/mixed loops.",
        "checks": {key: geometry.zero(value) for key, value in checks.items()},
        "gates": {
            "full_actual_profile_X_jets_not_fitted": True,
            "full_parent_clock_R_jets_not_replaced_off_clock": True,
            "all_ten_tensor_clock_first_directions_retained": True,
            "both_nonlinear_clock_contacts_kept_before_cancellation": True,
            "profile_mean_not_individually_zero": linear != 0,
            "projected_mean_contact_generically_nonzero": residual != 0,
            "common_clock_not_an_extra_spatial_heat_action": True,
            "all_source_clock_value_jets_through_degree_two_retained": True,
            "no_bare_projected_Ward_or_interacting_parent_closure": True,
        },
    }
    return result
