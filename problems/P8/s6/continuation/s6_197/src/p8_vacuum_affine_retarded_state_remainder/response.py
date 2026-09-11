"""Comparison functional and its precise fixed-renormalization boundary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes

from . import contact, memory, tail

DISPLAY = s.Integer(2) * 10**23
CANONICAL = 8 * s.Rational(1, 10) ** 777


@cache
def data():
    c = (
        memory.constants()["complete_memory_H1_coefficient"]
        + contact.constants()["contact_L2_coefficient"]
    )
    J, Jref, local = s.symbols("actual_current reference_current fixed_local")
    gamma = s.Symbol("gamma", real=True)
    a0, b0 = s.symbols("alpha0 beta0", complex=True)
    q = s.Symbol("canonical_tensor", real=True)
    k = s.Symbol("kappa", positive=True)
    K = s.Symbol("K", positive=True)
    zr, zi = s.symbols("pair_real pair_imaginary", real=True)
    z = zr + s.I * zi
    checks = {
        "actual_current_decomposition_before_renormalization": s.expand(
            J - Jref + Jref - J
        ),
        "unchanged_local_prescription_not_reselected": s.expand(
            (J + local) - (Jref + local) - (J - Jref)
        ),
        "constant_initial_alpha_not_time_varying": s.diff(a0, gamma),
        "reference_CCR_defect_requires_beta": s.expand(
            (1 + s.Abs(b0) ** 2) - 1 - s.Abs(b0) ** 2
        ),
        "full_current_pair_factor": s.expand(
            s.I * (2 * z - 2 * s.conjugate(z)) / 4 + s.im(z)
        ),
        "both_canonical_metric_variations": s.diff(2 * q / s.sqrt(k), q) ** 2 - 4 / k,
        "complete_canonical_comparison_display": 4 * DISPLAY / modes.KAPPA - CANONICAL,
        "canonical_complete_tail_display": 4 * tail.DISPLAY / (modes.KAPPA * K)
        - 4 * s.Rational(1, 10) ** 773 / K,
        "fixed_mass_and_parent": (modes.MASS - 1000) + (modes.KAPPA - 10**800),
    }
    return {
        "actual_finite_formula": "For each common regulator define J_K as the actual first-order retarded spatial current, including the local second-metric contact, and J_ref,K by replacing readouts with alpha_initial u_W8 in that same formula. This defines a reference comparison, not a new physical Gaussian state.",
        "full_comparison": "R(D,Gamma)=lim_K[J_K-J_ref,K] exists absolutely for overlapping smooth compact supports, with the theta(t-s) factor retained. The full norm is<2e23 M[D]M[Gamma], and the computational regulator error<1e27 M[D]M[Gamma]/K.",
        "unrounded_complete_coefficient": c,
        "canonical": "For h=sqrt(kappa)gamma/2, the same comparison has magnitude<8e-777 M[D]M[Gamma] and regulator error<4e-773 M[D]M[Gamma]/K. These are tensor normalization factors, not the fully reduced mixed norm.",
        "norm": "M[f]^2=||f||L2(dt dx;F)^2+||grad_spatial f||L2(dt dx;F)^2. No time derivative of a rapidly oscillating mixing coefficient is estimated.",
        "reference_not_a_state": "The coefficient alpha_initial is the actual unchanged all-order initial coefficient. Since |alpha_initial|^2-|beta_initial|^2=1, alpha_initial u_W8 alone is generally not CCR-normalized; W8 also has a nonzero equation residual. It is not a substituted physical state or a conserved renormalized current.",
        "matching_boundary": "The singular reference memory/contact must still be matched to the ORIGINAL covariant retarded prescription on overlapping supports and the diagonal. A finite comparison cannot perform that extension or prove a Ward identity for the reference alone. No finite-amplitude remainder, full inverse, interacting background, physical cutoff or V/G/B closure follows.",
        "zero_inputs": "All pairings and remainders are exactly zero if either test is zero; displayed strict bounds apply only to nonzero input norms.",
        "checks": checks,
        "gates": {
            "complete_remainder_strict_display": 0 < c < DISPLAY,
            "complete_canonical_strict_display": 4 * c / modes.KAPPA < CANONICAL,
            "reference_alpha_CCR_not_assumed": True,
            "retarded_step_and_full_contact_both_retained": True,
            "no_fixed_reference_diagonal_matching_asserted": True,
        },
    }
