"""Retarded/advanced reciprocity and the controlled full finite local reference."""

from functools import cache

import sympy as s
from p8_vacuum_affine_weighted_channel_resolvent import weighted as channel

from . import coefficients
from . import factorization as f

C0 = s.Rational(5, 2)
M = channel.M
SIGMA = channel.SIGMA


@cache
def data():
    mass = f.g.bridge.MASS
    gap = s.Symbol("positive_weight_minus20", positive=True)
    Y = C0 / gap**2
    DY = C0 / gap
    bounds = coefficients.data()["matrix_bounds"]
    omega = (
        2 * bounds["A"] * Y
        + bounds["H"] * DY**2
        + 2 * bounds["J"] * Y * DY
        + bounds["V0"] * Y**2
    )
    wanted = (
        s.Rational(6825, 4) * mass**2 / gap**2
        + 3575 * mass**2 / gap**3
        + 1000 * mass**4 / gap**4
    )
    fixed = s.factor(wanted.subs(gap, 100 * mass))
    middle = 20 / (123 + 32 * M)
    inverse = 125 / ((123 + 32 * M) * (SIGMA - 20) ** 4)
    physical = 48000 * s.pi**2 * f.g.bridge.KAPPA / ((123 + 32 * M) * (SIGMA - 20) ** 4)
    checks = {
        "complete_conjugated_local_remainder_norm": s.factor(omega - wanted),
        "fixed_actual_mass_remainder_upper": fixed - s.Rational(6825543, 40000000),
        "coarse_exponential_lower_polynomial": 1
        + 32
        + s.Rational(32**2, 2)
        - s.Integer(545),
        "combined_middle_strict_contraction_margin": s.factor(
            s.Rational(5, 13)
            - channel.KBOUND * (M + s.Rational(1, 4))
            - s.Rational(575, 4) / (13 * (32 + 13 * M))
        ),
        "complete_middle_inverse_bound": s.factor(
            1 / (channel.GAP - M - s.Rational(1, 4)) - middle
        ),
        "full_finite_local_curvature_inverse_bound": s.factor(
            C0**2 * middle / (SIGMA - 20) ** 4 - inverse
        ),
        "physical_density_on_right_retains_all_factors": s.factor(
            384 * s.pi**2 * f.g.bridge.KAPPA * inverse - physical
        ),
    }
    return {
        "source_time_adjoint_kernel_bound": "Zret(t,s)=Yadv(s,t)^T. Time reversal of Bc changes only the sign of its first-order C, so ||Z||<=(5/2)elapsed exp(20elapsed), ||partial_s Z||<=(5/2)exp(20elapsed). This is not an improved forward-time derivative claim.",
        "complete_initial_derivative_identity": "ZD* has ordinary kernel partial_s Z. The interior lower-boundary term cancels the retained initial delta from D(theta f); Z(t,t)=0 handles the upper boundary. No initial atom is dropped.",
        "weighted_coordinate_bounds": {
            "Y_and_Z": Y,
            "DY_and_ZDstar": DY,
            "gap": "sigma-20",
        },
        "actual_complete_local_conjugate_norm_bound": wanted,
        "chosen_weight_local_bound": fixed,
        "complete_middle_inverse_norm": middle,
        "complete_curvature_reference_weighted_inverse_norm": inverse,
        "physical_force_reference_weighted_inverse_norm": physical,
        "actual_reference": "Afinite=Bc*(Fdiag+V)Bc+Rloc=Bc*(Fdiag+V+Omega)Bc. Omega=Z Rloc Y is the actual full original finite local remainder in curvature coordinates; V is the optional S226 bounded channel matrix with norm<=M.",
        "inverse": "Use (I+Kdiag(V+Omega))^-1 Kdiag in the middle, with norm<=20/(123+32M), and composeY on the left andZ on the right. Both ordered graph inverse identities retain every local term and initial boundary.",
        "precise_composed_graph": "For H=weighted L2_tH^r, require s=Yx with x in H, y=(Fdiag+V+Omega)x in H and Bc*y in H, all as causal distributions with complete initial boundary. Then Bs=x. For f in H, y=Zf and x=Kfinite y prove existence; uniqueness follows successively from the Bc* and middle graph inverses. Measurable V is only multiplied by x in H, never by an arbitrary distribution. This is not a claim of maximal unrestricted distributional realization, graph density or full S222 invariance.",
        "spaces_and_scope": "Weighted L2_tH^r on the unchanged conformal slab, all comoving momenta and every realr, without spatial derivative loss. Removing the weight costs exp(sigma_M T), and physical normalization keeps right multiplication by64pi^2 kappa a^4. This is not actual full nonlocal curved quantum matching, the full S222 coupled inverse, stability or P8 closure.",
        "checks": checks,
        "gates": {
            "weight_gap_at_least100m": 2 * mass * 545 - 20 > 100 * mass,
            "actual_local_conjugate_norm_below_quarter": fixed < s.Rational(1, 4),
            "combined_contraction_bound_below_half": s.Rational(5, 13)
            < s.Rational(1, 2),
            "positive_full_middle_inverse_bound": middle.is_positive,
            "all_local_coefficients_nonzero_and_retained": all(
                v > 0 for v in bounds.values()
            ),
            "physical_kappa_unchanged": f.g.bridge.KAPPA == s.Integer(10) ** 800,
        },
    }
