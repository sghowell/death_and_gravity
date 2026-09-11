"""Coherent sourced Hadamard family with an unchanged explicit Cauchy covariance."""

from functools import cache

import sympy as s

from .bridge import KAPPA, MASS, ZETA


@cache
def data():
    omega = s.symbols("positive_frequency", positive=True)
    d = s.symbols("rate", real=True)
    v = 1 / s.sqrt(2 * omega)
    p = (-d - s.I * omega) * v
    phase = s.Matrix([v, p])
    C = phase * phase.conjugate().T
    mean = s.Matrix(s.symbols("mean_v mean_p", real=True))
    full = C + mean * mean.T
    Omega = s.Matrix([[0, 1], [-1, 0]])
    checks = {
        "positive_Cauchy_Gram_Hermiticity": C - C.conjugate().T,
        "positive_Cauchy_Gram_CCR": s.simplify(C - C.T - s.I * Omega),
        "coherent_connected_covariance_unchanged": full - mean * mean.T - C,
        "coherent_commutator_unchanged": full - full.T - (C - C.T),
        "W_covariance_normalization": KAPPA * ZETA / (KAPPA * ZETA) - 1,
    }
    return {
        "reference_state": "Exactly the S6.55 cutoff-summed all-order Proca Cauchy state at t0=-1/2, mass=1000, not its older finite-order comparison. Extend with the actual global CD Proca equation.",
        "admissible_histories": "Smooth real compact spacetime perturbations of the CD metric and u=t, supported strictly to the future of a common Cauchy neighborhood of t0; physical g remains globally hyperbolic and the full analytic coefficient domain is respected. S is then compact and smooth.",
        "family": "Keep the reference Cauchy covariance fixed, evolve the source-free Proca fluctuations on each metric, and add the unique retarded classical mean Abar=G_P,ret J with J=sqrt(kappa/zeta)S. This is a state on the sourced affine Proca algebra, identified with the homogeneous algebra by the classical shift.",
        "positivity_and_CCR": "The c-number shift preserves positivity and the field commutator; the connected covariance is exactly the source-free covariance. There is no assertion of a global preferred Fock implementation.",
        "Hadamard": "The source-free covariance is Hadamard by Proca time-slice propagation; the mean is smooth, so its tensor square adds no wavefront singularities.",
        "W_two_point_normalization": "Connected C_W=C_A/(kappa*zeta); the unscaled covariance cannot be assigned directly to W.",
        "primary_support": "Moretti-Murro-Volpe arXiv:2210.09278v3 Proposition4.7, with the Proca wavefront-set gap repaired by Fewster arXiv:2503.12544v2 Theorem5.1. The actual globally hyperbolic metrics, shared Cauchy neighborhoods and fixed positive mass are checked.",
        "mass": MASS,
        "checks": checks,
    }
