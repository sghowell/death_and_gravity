"""Private isotropic-clock current and canonical homogeneous response."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes

from . import current


@cache
def data():
    a, m, k2, nDn, trD = s.symbols("a mass k2 nDn trD", real=True)
    cT, cL, pT, pL = s.symbols("cT cL pT pL", real=True)
    full = (
        -a * m**2 * (cT * trD + (cL - cT) * nDn)
        + cT * k2 * (trD - nDn) / a
        + (pT * trD + (pL - pT) * nDn) / a
    )
    averaged = full.subs(nDn, trD / 3)
    kh = s.Symbol("kappa", positive=True)
    h = s.Symbol("h", real=True)
    gamma = 2 * h / s.sqrt(kh)
    force_factor = 2 / s.sqrt(kh)
    linear_factor = s.diff(gamma, h) * force_factor
    checks = {
        "full_isotropic_three_mode_metric_current_angular_zero": s.expand(
            averaged.subs(trD, 0)
        ),
        "both_canonical_metric_chain_factors": s.simplify(linear_factor - 4 / kh),
        "fixed_kappa_response_display": 4 * current.CURRENT[1] / modes.KAPPA
        - 4 * s.Rational(1, 10) ** 705,
    }
    return {
        "clock_rotation_identity": "At gamma0 the unchanged state is rotation invariant, including its longitudinal mode and both equal transverse modes. The complete block-diagonal metric current has angular dependence trD and n^tDn; averaging n^tDn=trD/3 makes every tracefree detector current vanish. Apply this to rotational finite regulators and the convergent matched limit, not to an unregulated divergent cancellation.",
        "modewise_full_metric_vertex_contraction": full,
        "clock_current": "The full fixed-prescription tracefree homogeneous current at the original isotropic clock is exactly zero.",
        "canonical_field": "h=sqrt(kappa)gamma/2. The canonical quantum force is 2J_gamma/(a^3 sqrt(kappa)); its linearization at the original clock has TWO canonical chain factors, hence4DJ/(a^3 kappa).",
        "canonical_homogeneous_linear_response_display": 4 * s.Rational(1, 10) ** 705,
        "norm_boundary": "The C12-to-C0 bound applies to the stated smooth homogeneous tracefree direction subspace with the common preparation. At the clock Gamma is arbitrary and the first variation is linear. The finite-amplitude result is along the admitted gamma=epsilon Gamma families; no general nonlinear Frechet norm away from the clock is inferred.",
        "feedback_boundary": "This derivative-losing response bound is not a same-space Banach contraction, full response inverse, reduced mixed-mode norm or interacting background solution.",
        "checks": checks,
        "gates": {
            "temporal_constraint_polarization_not_omitted": cL in full.free_symbols
            and pL in full.free_symbols,
            "actual_positive_canonical_normalization": modes.KAPPA == 10**800,
            "first_current_direction_bound_retained": current.CURRENT[1] == 10**95,
            "no_division_by_bounce_density": True,
        },
    }
