"""Original CD/M1 comparison, with the free chi action and physical metric fixed.

Operator numbers use the immutable witness's M=tau=1 clock normalization.
The Hubble/null-equation budgets below restore arbitrary physical tau and
the matched effective Planck mass. No field/potential content is silently
identified with the pinned one-scalar composite family.
"""

from fractions import Fraction as Q
from functools import cache

import sympy as sp
from p8 import jets, rational_candidates

from .branch import rational


def null_residual_budget(epsilon_hdot, epsilon_chidot):
    """Exact positive actual-null-residual budget in the stated solution norm.

    The inputs bound |tau² Hdot-4| and |tau chidot/M_target-1/10|
    at the centre. They are not per-operator coefficient estimates.
    Require 0<=epsilon_hdot<4 and 0<=epsilon_chidot<1/10: the latter
    condition is needed before squaring the lower bound on |chidot|.
    """
    epsilon_hdot, epsilon_chidot = map(rational, (epsilon_hdot, epsilon_chidot))
    if not 0 <= epsilon_hdot < 4 or not 0 <= epsilon_chidot < Q(1, 10):
        raise ValueError("Require 0<=epsilon_hdot<4 and 0<=epsilon_chidot<1/10")
    return 2*(4-epsilon_hdot)+(Q(1, 10)-epsilon_chidot)**2


@cache
def derive():
    old = rational_candidates.specification("CD_matter")
    x, phi = rational_candidates.X, jets.t
    tau, planck = sp.symbols("tau M_target", positive=True)
    t = sp.Symbol("physical_T", real=True)
    h = 4*t/(tau**2+t**2)
    chidot = planck/(10*tau*(1+(t/tau)**2)**6)
    eps_h, eps_chi = sp.symbols("epsilon_Hdot epsilon_chidot", nonnegative=True)
    # -2M_target² Hdot = n_clock + chidot² + residual_null.
    # For the exact target, n_clock>=0 implies a negative order-one residual.
    return {"old": old, "X": x, "phi": phi, "tau": tau, "M_target": planck, "T": t,
            "target_F2X": sp.diff(old["F2"], x), "target_A3": old["A3"],
            "zero_relative_F2X": sp.S.Zero, "zero_relative_A3": sp.S.Zero,
            "target_F2X_at_centre": sp.diff(old["F2"], x).subs({phi: 0, x: 1}),
            "target_A3_at_centre": old["A3"].subs({phi: 0, x: 1}),
            "target_H": h, "target_chidot": chidot,
            "target_H_left": sp.simplify(h.subs(t, -tau/2)),
            "target_H_right": sp.simplify(h.subs(t, tau/2)),
            "endpoint_H_error_threshold": sp.Rational(8, 5),
            "exact_null_residual_lower_magnitude": sp.Rational(801, 100),
            "epsilon_Hdot": eps_h, "epsilon_chidot": eps_chi,
            "approximate_null_residual_lower_magnitude": 2*(4-eps_h)+(sp.Rational(1, 10)-eps_chi)**2}


@cache
def checks():
    d = derive()
    old, phi, tau, t, planck = (d[key] for key in ("old", "phi", "tau", "T", "M_target"))
    nclock = sp.Symbol("n_clock", nonnegative=True)
    residual = -2*planck**2*sp.diff(d["target_H"], t)-d["target_chidot"]**2-nclock
    return {
        "old_CD_X_curvature_jet": d["target_F2X_at_centre"]+sp.Rational(1, 2),
        "old_CD_DHOST_jet": d["target_A3_at_centre"]-1,
        "old_free_M1_charge_is_nonzero_constant": sp.cancel(old["a"]**3*old["chi_dot"]-sp.Rational(1, 10)),
        "physical_time_H_left": d["target_H_left"]+8/(5*tau),
        "physical_time_H_right": d["target_H_right"]-8/(5*tau),
        "physical_time_Hdot_at_bounce": sp.diff(d["target_H"], t).subs(t, 0)-4/tau**2,
        "physical_time_free_chi_charge_restoration": sp.cancel(
            (1+(t/tau)**2)**6*d["target_chidot"]-planck/(10*tau)),
        "actual_null_residual_with_exact_free_M1": sp.simplify(
            residual.subs(t, 0)+nclock+d["exact_null_residual_lower_magnitude"]*planck**2/tau**2),
        "approximate_null_budget_reduces_to_exact": d["approximate_null_residual_lower_magnitude"].subs(
            {d["epsilon_Hdot"]: 0, d["epsilon_chidot"]: 0})-d["exact_null_residual_lower_magnitude"],
        "target_background_is_pinned_CD_after_time_restore": sp.cancel(
            d["target_H"]-old["H"].subs(phi, t/tau)/tau),
    }


def controls():
    return {
        "required_F2X_remainder_centre": Q(1, 2),
        "required_A3_remainder_centre": Q(1),
        "static_GR_endpoint_error_threshold": Q(8, 5),
        "constant_vacuum_vs_pinned_bounce_potential_minimum_gap": Q(3, 8)-Q(1, 10000),
        "dropping_free_chi_understates_null_residual_budget": Q(1, 100),
        "NEC_violating_source_can_reverse_GR_Hdot": Q(-1),
        "wrong_physical_metric_or_clock_invalidates_coefficient_budget": "Contract condition, not an invariant off-frame coefficient theorem",
    }
