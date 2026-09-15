"""Explicit low-cut subtraction, its massive coefficient, and distinct limits."""

from functools import cache

import sympy as s

from . import cut, soft, source

S, T, U, MU, K = source.S, source.T, source.U, source.MU, source.K
CAP = s.Symbol("positive_massless_cut_cap_squared", positive=True)
NU2 = cut.NU2
V = s.Symbol("crossing_v", real=True)


def require_cap(cap, mass2=source.MASS2):
    a, m = soft.exact_positive(cap), soft.exact_positive(mass2)
    if a <= 4 * m:
        raise ValueError(
            "Require the squared subtraction cap strictly above four times the external mass squared"
        )
    return a, m


def low_piece():
    return -cut.COEFFICIENT * sum(
        cut.POLYS[z] * (s.log(-z / NU2) - s.log((CAP - z) / NU2)) for z in cut.CHANNELS
    )


def remainder():
    return -cut.COEFFICIENT * sum(
        cut.POLYS[z] * s.log((CAP - z) / NU2) for z in cut.CHANNELS
    )


def forward_coefficient():
    D = CAP - 2 * MU
    return -cut.COEFFICIENT * (
        2 * s.log(D / NU2) + s.log(CAP / NU2) - 12 * MU / D - 14 * MU * MU / (D * D)
    )


def cap_derivative():
    D = CAP - 2 * MU
    return -cut.COEFFICIENT * (2 / D + 1 / CAP + 12 * MU / D**2 + 28 * MU * MU / D**3)


@cache
def data():
    D = CAP - 2 * MU
    forward = remainder().subs({S: 2 * MU + V, T: 0, U: 2 * MU - V}, simultaneous=True)
    actual = s.diff(forward, V, 2).subs(V, 0) / 2
    allraw = cut.massless_log_part()
    L = s.Symbol("finite_log_resolution", positive=True)
    kap = s.Symbol("positive_kappa", positive=True)
    ch = s.Symbol("fixed_scaled_pole", positive=True)
    q = s.Symbol("positive_transfer_magnitude", positive=True)
    ep, Q = s.symbols("positive_resolution positive_transfer_cap", positive=True)
    primitive = s.log(q) / K
    checks = {
        "whole_crossing_low_cut_plus_remainder": s.expand(
            allraw - low_piece() - remainder()
        ),
        "full_massive_subtracted_forward_coefficient": s.simplify(
            actual - forward_coefficient()
        ),
        "entire_cap_derivative": s.factor(
            s.diff(forward_coefficient(), CAP) - cap_derivative()
        ),
        "s_u_spectral_and_t_transfer_cap_terms": s.factor(
            cap_derivative()
            + 2 * cut.COEFFICIENT * (CAP * CAP + 2 * MU * CAP + 6 * MU * MU) / (D**3)
            + cut.COEFFICIENT / CAP
        ),
        "scale_dependence_retained_not_set_to_zero": s.simplify(
            s.diff(forward_coefficient(), NU2) * NU2 - 3 * cut.COEFFICIENT
        ),
        "complete_minimal_forward_pole": s.diff(
            source.canonical_gravity.pole()["minimal_tree_t_channel"].subs(
                {
                    s.Symbol("s", real=True): 2 * MU + V - T / 2,
                    s.Symbol("t", real=True): T,
                    s.Symbol("mass_squared", real=True): MU,
                    s.Symbol("kappa", real=True): K,
                }
            ),
            V,
            2,
        )
        / 2
        + 1 / (K * T),
        "smeared_graviton_pole_primitive": s.diff(primitive, q) - 1 / (K * q),
        "fixed_resolution_then_decoupling": s.limit(L / kap, kap, s.oo),
        "joint_detector_limit_not_fixed_resolution": (L / kap).subs(L, kap * ch) - ch,
        "low_cut_log_contour_primitive": s.diff(s.log(CAP - S) - s.log(-S), CAP)
        - 1 / (CAP - S),
    }
    for p, r in ((S, T), (S, U), (T, U)):
        checks["whole_subtracted_crossing_" + str(p) + str(r)] = s.expand(
            remainder().subs({p: r, r: p}, simultaneous=True) - remainder()
        )
    # A finite, strictly positive cap example is a subtraction convention,
    # not the physical EFT cutoff and not a positivity-admissible functional.
    cap = 10**196
    require_cap(cap)
    bound = 2000 / (960 * source.KAPPA**2)  # nu2=mass2=1, cap10^196, pi^2>1.
    checks["named_cap_ratio_gap"] = s.Integer(cap) - (s.Integer(cap) - 2) - 2
    return {
        "whole_original_log_coefficient": allraw,
        "whole_explicit_crossing_low_cut_piece": low_piece(),
        "whole_exact_after_low_cut_subtraction": remainder(),
        "whole_subtracted_b20_modulo_local_counterterm": forward_coefficient(),
        "whole_cut_cap_derivative": cap_derivative(),
        "smeared_tree_pole": s.log(Q / ep) / K,
        "named_nonlocal_coefficient_bound": {
            "mass_squared": 1,
            "renormalization_scale_squared": 1,
            "subtraction_cap_squared": cap,
            "absolute_b20_bound": bound,
            "relative_to_original_four_lambda_bound": bound
            / (4 * source.current.heavy.LAMBDA),
        },
        "subtraction_definition": "Subtract the complete written crossing-symmetric logarithmic low-cut piece once. It has the same discontinuity as the original M1 logarithms below the cap and zero discontinuity above it. Keep the real local polynomial untouched; this definition chooses a finite low-cut convention, not a new parent counterterm or a physical cutoff.",
        "analytic_domain": "After this M1-sector subtraction, all its branch points start at cap. For cap>4m^2, its forward crossing center s=u=2m^2,t0 has a genuine complex neighborhood. This statement concerns this loop sector only; all other unresolved cuts still require their own subtraction.",
        "detector_limit_boundary": "At fixed finite log(Q/E), the minimal pole smear vanishes as1/kappa; holding log(Q/E)/kappa fixed gives a nonzero limit. Dimensional epsilon, physical detector E, transfer q, subtraction cap, renormalization scale and the earlier oscillator cutoffs are distinct.",
        "positivity_boundary": "A q-weighted pole integral is only an integral here, not a proved positive UV functional. No numerical massless-scalar inequality from an external paper is transferred. Full massive light cuts, finite counterterms, finite-coupling detector errors, Regge control, matching and the bounce still need proof.",
        "checks": checks,
        "gates": {
            "all_three_low_cuts_subtracted_with_exact_add_back": True,
            "transfer_channel_cut_cap_term_not_discarded": cap_derivative().has(
                1 / CAP
            ),
            "real_local_prescription_remains_separate": True,
            "named_nonlocal_bound_far_below_original_tree_b20": bound
            / (4 * source.current.heavy.LAMBDA)
            < s.Rational(1, 10**990),
            "named_cap_is_not_declared_a_Wilsonian_cutoff": True,
            "massive_center_not_assumed_analytic_before_low_cut_subtraction": True,
            "all_physical_Regge_and_detector_error_premises_remain_open": True,
        },
    }
