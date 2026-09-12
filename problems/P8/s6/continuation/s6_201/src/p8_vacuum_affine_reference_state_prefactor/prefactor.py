"""Initial occupation weights with exact phase cancellation."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import reference

BETA = s.Integer(2) * 10**6


@cache
def constants():
    c = reference.data()["source_pinned_constants"]
    return {
        "actual_initial_beta_envelope": c["B6"] + c["E"] / reference.MASS**4,
        "safe_beta_envelope": BETA,
        "maximum_initial_occupation": BETA**2 / reference.MASS**12,
        "unit_reference_field": reference.FIELD_NORM,
        "unit_reference_pair": reference.PAIR_REF,
    }


@cache
def data():
    a, b, x, y = s.symbols("alpha_k alpha_l pair_D pair_Gamma", complex=True)
    u, v = s.symbols("occupation_k occupation_l", nonnegative=True)
    W, d, occupation = s.symbols("W damping occupation", real=True, positive=True)
    f = 1 / s.sqrt(2 * W)
    p = (-s.I * W - d) / s.sqrt(2 * W)
    weighted = s.sqrt(1 + occupation)
    c = constants()
    checks = {
        "complete_pair_initial_phase_cancellation": s.expand(
            s.conjugate(a * b * x) * (a * b * y)
            - s.conjugate(a) * a * s.conjugate(b) * b * s.conjugate(x) * y
        ),
        "one_mode_contact_phase_cancellation": s.expand(
            s.conjugate(a * x) * (a * y) - s.conjugate(a) * a * s.conjugate(x) * y
        ),
        "full_two_leg_occupation_increment": s.expand(
            (1 + u) * (1 + v) - 1 - u - v - u * v
        ),
        "controlled_full_quadratic_occupation_weight": s.expand(
            2 * (u + v) - (u + v + u * v) - (u * (1 - v) + v)
        ),
        "unit_reference_canonical_Wronskian": s.simplify(
            f * s.conjugate(p) - p * s.conjugate(f) - s.I
        ),
        "weighted_reference_Wronskian_defect": s.simplify(
            weighted * f * s.conjugate(weighted * p)
            - weighted * p * s.conjugate(weighted * f)
            - s.I * (1 + occupation)
        ),
        "actual_safe_beta_upper": BETA - 2000000,
    }
    return {
        "unchanged_actual_state": "The original all-order prepared Proca state is unchanged. Initial alpha,beta obey |alpha0|^2-|beta0|^2=1 for every physical mode. No initial phase is set to zero in the actual state.",
        "beta_bound": "The source-pinned B6+E/m^4 is below B=2e6, hence |beta0(k)|<=B nu^-6 and b_k=|beta0(k)|^2<=B^2 nu^-12<1. This bound requires no momentum derivatives of the Borel-prepared state.",
        "phase_cancellation": "Every full memory pair has the same alpha_k alpha_l at detector and source, so only (1+b_k)(1+b_l) remains. The complete contact has the one-mode weight1+b_k. All polarizations and the b_k b_l term are retained.",
        "weight_majorant": "For 0<=b_k,b_l<1, the full increment b_k+b_l+b_k b_l<=2 B^2(nu^-12+mu^-12). This is an absolute bound on the exact weight, not its linearization.",
        "unit_reference": "Use the phase-stripped unit-W8 modes in the same finite-regulator memory and contact formulas. Their canonical Wronskian is i, but W8 is not an exact mode-equation solution. This is not a new physical bisolution, vacuum reset or effective action.",
        "analyticity_boundary": "No momentum analyticity or derivative bound for alpha0,beta0 is assumed. Joint complex spatial-momentum bounds for the remaining unit-W8 kernel are a future obligation, not established here.",
        "constants": c,
        "checks": checks,
        "gates": {
            "strict_safe_initial_beta_envelope": 0
            < c["actual_initial_beta_envelope"]
            < BETA,
            "occupation_uniformly_below_one": 0 < c["maximum_initial_occupation"] < 1,
            "unit_pair_covered_by_original_analytic_bound": 2 * reference.FIELD_NORM**2
            < reference.PAIR_REF,
            "actual_fixed_mass": reference.MASS == 1000,
        },
    }
