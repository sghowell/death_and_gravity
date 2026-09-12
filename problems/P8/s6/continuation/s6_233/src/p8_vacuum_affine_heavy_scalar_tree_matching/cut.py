"""Full rational angular integrals and first elastic comparison, not a full-loop error."""

from functools import cache

import sympy as s

from . import matching, model

a, k, g2, d, x = s.symbols(
    "angular_a angular_k coupling_squared angular_d angle_cosine", real=True
)
A = d + 2 * g2 * a / (a * a - k * k * x * x)
LOG = s.log((a + k) / (a - k))
I1 = 2 * d + 2 * g2 * LOG / k
I2 = (
    2 * d * d
    + (4 * d * g2 / k + 2 * g2 * g2 / (a * k)) * LOG
    + 4 * g2 * g2 / (a * a - k * k)
)


@cache
def data():
    logarithm = s.log((a + k * x) / (a - k * x))
    P1 = d * x + g2 * logarithm / k
    P2 = (
        d * d * x
        + (2 * d * g2 / k + g2 * g2 / (a * k)) * logarithm
        + 2 * g2 * g2 * x / (a * a - k * k * x * x)
    )
    ratio = k * x / a
    series = d + 2 * g2 / a * sum(ratio ** (2 * n) for n in range(4))
    geometric_rest = 2 * g2 / a * ratio**8 / (1 - ratio * ratio)
    S = s.Symbol("physical_s", positive=True)
    kval = (S - 4) / 2
    aval = model.MASS2 + kval
    dval = model.CONTACT + model.G2 / (model.MASS2 - S)
    return {
        "full_rational_angular_form": A,
        "physical_substitution": {"k": kval, "a": aval, "d": dval, "g2": model.G2},
        "complete_angular_integral": I1,
        "complete_squared_angular_integral": I2,
        "threshold_integrals": (2 * d + 4 * g2 / a, 2 * (d + 2 * g2 / a) ** 2),
        "normalization": "t0_tree=beta*I1/(64pi) and rho_first=beta*I2/(64pi), beta=sqrt(1-4/s). The original identical state and phase-space factors remain. At threshold beta0, the displayed integral limits are finite.",
        "all_partial_waves_scope": "For s>4, k>0 and a-k=M_H^2>0. The two exchange denominators have a convergent even angular power series with infinitely many nonzero even Legendre coefficients; the contact and s-channel term affect only l0. No l0/l2 truncation is used for the new model.",
        "first_absorptive_scope": "In the named window s<=10^196<M_H^2, the heavy single pole and the HH pair cut are inaccessible, and phi-parity forbids the mixed phi-H two-body channel. The first elastic absorptive coefficient uses the complete new tree in the light two-body cut. This is not the full real loop amplitude or a bounded all-loop error.",
        "first_elastic_comparison": "For4<s<=10^196, rho_first_original<=rho_first_V2S_T1<(61/60)^2 rho_first_original; both vanish at the massive threshold. Its relative difference is less than121/3600. This follows from the all-angle positive tree matching, not from replacing the new rational amplitude by its original quartic truncation.",
        "relative_first_elastic_error_upper": s.Rational(121, 3600),
        "checks": {
            "first_angular_primitive": s.cancel(s.diff(P1, x) - A),
            "full_squared_angular_primitive": s.cancel(s.diff(P2, x) - A * A),
            "first_integral_threshold_limit": s.cancel(
                s.limit(I1, k, 0) - 2 * d - 4 * g2 / a
            ),
            "squared_integral_threshold_limit": s.cancel(
                s.limit(I2, k, 0) - 2 * (d + 2 * g2 / a) ** 2
            ),
            "angular_exchange_evenness": s.cancel(A - A.subs(x, -x)),
            "complete_angular_geometric_remainder": s.cancel(
                A - series - geometric_rest
            ),
            "physical_positive_denominator_gap": aval - kval - model.MASS2,
            "first_elastic_relative_error_constant": (1 + matching.RELATIVE_BOUND) ** 2
            - 1
            - s.Rational(121, 3600),
        },
        "gates": {
            "new_heavy_pair_cut_above_named_window": 4 * model.MASS2
            > matching.ENERGY**2,
            "new_heavy_single_pole_above_named_window": model.MASS2
            > matching.ENERGY**2,
            "positive_physical_exchange_denominator_gap": model.MASS2 > 0,
            "no_l0_l2_only_cut_for_rational_model": True,
            "first_elastic_comparison_not_full_quantum_error": True,
        },
    }
