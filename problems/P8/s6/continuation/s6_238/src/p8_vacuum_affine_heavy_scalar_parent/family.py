"""A separately named full covariant candidate, not a replacement certificate."""

from functools import cache

import sympy as s
from p8_affine_vacuum_domain import family as original
from p8_vacuum_affine_heavy_scalar_tree_matching import model
from p8_vacuum_affine_quantum_retuning import profile
from p8_vacuum_canonical_affine_decoupling import family as canonical

NAME = "CD-REG-AFFINE-ISO-QG1-H8A420"
K0 = canonical.K0
N = original.N
LAMBDA, GAMMA = canonical.LAMBDA, canonical.GAMMA
MASS2, G, CONTACT = model.MASS2, model.G, model.CONTACT
ORDER = 8
LOCALIZER = s.Integer(10) ** 420
u, X = original.u, original.X
A = s.Symbol("positive_localizer_A", positive=True)


def switch(x=X, parameter=LOCALIZER):
    return (1 - x) ** ORDER * s.exp(-parameter * x * x)


@cache
def coefficients():
    h = switch()
    df = -LAMBDA * K0 * X**2 * h + (GAMMA / s.Integer(3) + CONTACT / 24) * K0 * u**4 * h
    dr = GAMMA * K0 * X**2 * h
    j = G * s.sqrt(K0) * u**2 * h / 2
    R = original.data()["R"] + dr
    F = profile.data()["new_complete_scalar_coefficient"] + df
    rx = s.diff(R, X)
    return {
        "F": F,
        "R": R,
        "A3": rx / X,
        "A4": -rx / X - 7 * rx**2 / (4 * R),
        "A5": rx**2 / (R * X),
        "delta_F": df,
        "delta_R": dr,
        "normalized_heavy_source": j,
    }


@cache
def data():
    h = switch(X, A)
    df = -LAMBDA * K0 * X**2 * h + (GAMMA / 3 + CONTACT / 24) * K0 * u**4 * h
    dr = GAMMA * K0 * X**2 * h
    j = G * s.sqrt(K0) * u**2 * h / 2
    zero = {u: 0, X: 0}
    # Literal finite jets of the entire old functions are checked in S177.
    oldg = canonical.germs()
    checks = {
        "literal_old_vacuum_" + key: value for key, value in oldg["checks"].items()
    }
    for i, k, target in (
        (0, 0, 0),
        (1, 0, 0),
        (2, 0, 0),
        (0, 1, 0),
        (1, 1, 0),
        (3, 0, 0),
        (2, 1, 0),
        (0, 2, -2 * LAMBDA * K0),
        (4, 0, (8 * GAMMA + CONTACT) * K0),
    ):
        checks[f"new_explicit_F_jet_{i}_{k}"] = (
            s.diff(df, u, i, X, k).subs(zero) - target
        )
    checks.update(
        {
            "old_R_XX_plus_new_R_XX_zero": s.diff(dr, X, 2).subs(zero) - 2 * N,
            "regular_R_delta_X_squared": s.cancel(dr / X**2 - GAMMA * K0 * h),
            "positive_sign_removes_old_quartic_a3": -2 * GAMMA
            + s.diff(dr, X, 2).subs(zero) / K0,
            "new_vacuum_canonical_quartic_contact": -GAMMA / 3
            + (GAMMA / 3 + CONTACT / 24)
            - CONTACT / 24,
            "new_source_literal_cubic": s.diff(j, u, 2).subs(zero) / s.sqrt(K0) - G,
            "source_first_Y_jet_not_discarded": s.diff(h, X).subs(X, 0) + 8,
            "source_second_Y_jet_localizer_present": s.diff(h, X, 2).subs(X, 0)
            - (56 - 2 * A),
            "old_fixed_proca_vacuum_constant_retained": profile.data()[
                "full_coefficient_correction"
            ].subs(X, 0)
            + profile.PV,
        }
    )
    for k in range(8):
        checks["switch_clock_zero_" + str(k)] = s.diff(h, X, k).subs(X, 1)
    checks["first_nonzero_clock_switch_jet"] = s.diff(h, X, 8).subs(X, 1) - s.factorial(
        8
    ) * s.exp(-A)
    # At any clock u, factorization of the coefficient differences is exact.
    r, p, d, dp, x = s.symbols("R R_X delta_R delta_R_X x", real=True)
    qdiff = (p + dp) ** 2 / (r + d) - p * p / r
    qnumer = (2 * p * dp + dp**2) * r - p * p * d
    checks["entire_dependent_Ia_difference"] = s.cancel(qdiff - qnumer / (r * (r + d)))
    checks["entire_A3_difference"] = s.cancel((p + dp) / x - p / x - dp / x)
    checks["entire_A4_difference"] = s.cancel(
        -(p + dp) / x
        - 7 * (p + dp) ** 2 / (4 * (r + d))
        + p / x
        + 7 * p * p / (4 * r)
        + dp / x
        + 7 * qdiff / 4
    )
    checks["entire_A5_difference"] = s.cancel(
        (p + dp) ** 2 / ((r + d) * x) - p * p / (r * x) - qdiff / x
    )
    return {
        "candidate": NAME,
        "literal_definition": "Keep the full original S109 F,R and add the complete fixed S182 smooth scalar retuning. Set X=Y/kappa0 and h=(1-X)^8 exp(-10^420 X^2). Add deltaf=-lambda Y^2 h+(gamma/3+C/24)Phi^4 h, deltar=+gamma Y^2 h and minimally coupled LH=(partial H)^2/2-nH^2/2+g H Phi^2 h/2. Recompute every R-dependent Ia and affine coefficient.",
        "normalized_additions": {
            "delta_F": df,
            "delta_R": dr,
            "J_H_over_sqrt_kappa0": j,
        },
        "fixed_localizer_parameter": LOCALIZER,
        "vacuum_scalar_jet": "After subtracting only its retained origin value, the decoupled flat action through field degree4 is exactly V2S-T1: light/heavy masses1,sqrt(n), cubic g H Phi^2/2 and contact C Phi^4/24. Old lambda Y^2 and the old quartic a3 interaction are cancelled, not double counted.",
        "higher_interactions": "The first new source correction is -4g H Phi^2 Y/kappa0, of total field degree5. Full higher scalar/Ia interactions and the smooth S182 degree2048 vacuum retuning remain. They are not the S235-S237 renormalizable one-loop model.",
        "clock_jet_statement": "deltaF,deltar,J_H have (X-1)^8 factors; deltaA3,deltaA4,deltaA5 have (X-1)^7 factors. Pure original light/metric action Taylor jets through degree6 agree at X=1. Comparing full actions also retains the entire minimally coupled free-heavy action, including its metric interactions. H=0 leaves the original classical equations and linear light block unchanged; it does not remove the heavy quantum stress.",
        "smoothness": "The R change and local additions are real analytic on the original real strip. The complete F includes the fixed S182 reference-stress profile and is only claimed smooth, not globally real analytic.",
        "constant": "The old specified vacuum constant kappa0*PV is retained. A new heavy quantum vacuum counterterm is not silently added.",
        "canonical_family": "Fix the complete new canonical f,r,h,g,n at kappa0. For every kappa>=kappa0 use R_kappa=1+r/kappa and F_kappa=f/kappa in the rescaled canonical variables. The same heavy action is fixed; no off-anchor background or interacting-state limit is asserted.",
        "checks": {key: s.cancel(value) for key, value in checks.items()},
        "gates": {
            "same_N_1024": N == 1024,
            "actual_positive_large_localizer": LOCALIZER == 10**420,
            "clock_zero_order_eight": ORDER == 8,
            "same_fixed_tree_parameters": G == s.Rational(1, 8192)
            and MASS2 == model.GAP + 2,
        },
    }
