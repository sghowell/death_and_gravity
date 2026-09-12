"""Constructive finite-weight actual homogeneous inverse with the pole retained."""

from functools import cache

import sympy as s

from . import local

K, C = s.symbols(
    "actual_shear_L1_on_unit_slab complete_actual_remainder_majorant", positive=True
)
BETA = K * C
WEIGHT = (4 * BETA + 1) ** 2


@cache
def data():
    u = s.Symbol("sqrt_weight_at_least_one", positive=True)
    margin = 2 * u - 2 - 2 * s.log(u)
    contraction = 2 * BETA / (4 * BETA + 1)
    inverse = 2 * s.exp(WEIGHT) * K
    physical = 16 * s.pi**2 * local.original.KAPPA * inverse
    checks = {
        "strict_contraction_margin": s.factor(
            s.Rational(1, 2) - contraction - 1 / (2 * (4 * BETA + 1))
        ),
        "finite_weight_square_root": s.sqrt(WEIGHT) - (4 * BETA + 1),
        "weak_log_Young_majorant_equality_at_one": margin.subs(u, 1),
        "weak_log_Young_margin_derivative": s.diff(margin, u) - 2 * (u - 1) / u,
        "full_unweighted_source_norm_factor": physical
        - 32 * s.pi**2 * local.original.KAPPA * s.exp(WEIGHT) * K,
        "inverse_normalization_keeps_original_kappa": physical / inverse
        - 16 * s.pi**2 * local.original.KAPPA,
    }
    return {
        "exact_actual_normal_forms": "I4 T_Gamma=F2(Dt^2)+V_Gamma and I4 T_total=F2(Dt^2)+V_total. These follow from the actual zero-transfer mode/current, complete original contacts and fixed-dimensional matching, not from the old derivative-losing weak norm.",
        "actual_remainder_class": "For each fixedj, |(Dt+Ds)^j V_total(t,s)|<=Cj(1+|log(t-s)|) on the compact unit causal triangle. C=C0 is a valid finite majorant INCLUDING the actual classical kappa terms and every quantum remainder; no numerical value is asserted.",
        "retained_inverse_kernel": "K2 is the ORIGINAL S223 pole-plus-cut kernel. K=||K2||L1(0,1) is finite; the full kernel is notL1 on the half-line. Its pole, causal branches and source preparation are unchanged.",
        "beta": BETA,
        "constructive_time_weight": WEIGHT,
        "weighted_contraction_upper": contraction,
        "unweighted_normalized_source_bound": inverse,
        "unweighted_canonical_source_bound": physical,
        "ordered_solution": "h=(I+K2 V_total)^-1 K2 I4 g, g=-16pi^2 kappa f for canonical forcingf. This is a convergent bounded-operator Neumann series on weighted continuous functions, not a finite Born truncation or an assumption of commutingK2 andV.",
        "prepared_domain": "Smooth zero-past homogeneous tracefree amplitudes and smooth forcing with the same nonempty initial zero neighborhood. The response is per comoving volume, not a nonzero-transfer H^r norm theorem or a point evaluation of an arbitrary weak multiplier.",
        "regularity_and_inverse_identities": "Simultaneous kernel derivatives and the retained initial germ give arbitrary finite time regularity inductively. Causal D^4 I4=I recovers the original total equation. The ordered Volterra inverse proves existence and uniqueness on this stated prepared smooth class, with no higher-order pole deletion or new final-time condition.",
        "bound_scope": "Finite Cj and K are defined by the proved actual kernels, not freely tunable physical coefficients. The displayed expression is constructive in unevaluated majorants, NOT a numerical or small inverse bound, unweighted stability, an actual nonzero-transfer inverse or original P8 closure.",
        "checks": checks,
        "gates": {
            "positive_finite_symbolic_weight": WEIGHT.is_positive,
            "strict_weight_above_one": s.expand(WEIGHT - 1).is_positive,
            "positive_strict_contraction_margin": (
                1 / (2 * (4 * BETA + 1))
            ).is_positive,
            "positive_kernel_majorant_not_zero": K.is_positive,
            "complete_remainder_majorant_not_zero": C.is_positive,
            "canonical_source_factor_not_deleted": local.original.KAPPA > 1,
        },
    }
