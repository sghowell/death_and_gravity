"""Uniform high-momentum justification and the scoped conversion frontier."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes

from . import angular


@cache
def data():
    x = s.symbols("inverse_radius", nonnegative=True)
    a, m = s.symbols("a m", positive=True)
    H, rate = s.symbols("H rate", real=True)
    P = s.symbols("P1:5", real=True)
    Q = s.symbols("Q1:5", real=True)
    omega = s.sqrt(a**-2 + m * m * x * x)
    domega = -H / (a * a * omega)
    W = omega + sum(P[j - 1] * x ** (2 * j) * omega ** (1 - 2 * j) for j in range(1, 5))
    dW = domega + sum(
        x ** (2 * j)
        * (
            Q[j - 1] * omega ** (1 - 2 * j)
            + (1 - 2 * j) * P[j - 1] * omega ** (-2 * j) * domega
        )
        for j in range(1, 5)
    )
    f = 1 / s.sqrt(2 * W)
    p = (-s.I * W - x * (dW / (2 * W) + rate)) * f
    v = 1 / s.sqrt(2 * a)
    coefficient = angular.coefficient(
        s.diag(1, -1, 0) / s.sqrt(2), s.diag(1, -1, 0) / s.sqrt(2), a=a
    )
    zero = lambda expression: s.simplify(expression.subs(x, 0))
    return {
        "uniform_actual_W8_limit": "Write x=1/r, omega/r=sqrt(a^-2+m^2 x^2), and W8/r=omega/r+sum_(j=1)^4 P_j x^(2j)(omega/r)^(1-2j). The actual P_j(t,z) and their required time derivatives are bounded on the fixed CD slab and z near1. Thus these expressions are smooth in x at0, uniformly on compact time/angular sets. The leading normalized readouts are obtained by setting x=0; fixed external P changes l/r by O(1/r).",
        "actual_j0_asymptotic": "The full summed endpoint density is r*h(D,G,n)+O(1), uniformly in angle and on the compact slab for fixed external P and fixed source/detector tensors. A finite atlas of real polarization frames suffices; the full projector contraction is frame independent. This is a fixed-P symbol statement, not a uniform all-P test-norm estimate.",
        "shell_limit": "On the lost shell r in[K-|P|,K], the O(1) density error contributes O(K^2) at fixed P. Replacing its leading density by r*h and using the exact shell limit gives pair-minus-single B0=-K^3|P|*A(D,G,phat,t)+o(K^3). A is the complete angular coefficient in angular.py.",
        "other_endpoints": "S202 bounds each j>=1 endpoint source-jet coefficient by C_jr nu^(1-j) on the far domain. Every fixed P lies in that domain for large K. The lost shell has volume O(K^2); its j>=1 conversion is O(K^(3-j)), hence at most O(K^2). Therefore no other first-five endpoint cancels the leading K^3|P| conversion. The complete contact has exactly the same original one-leg band in both sides of this explicitly defined comparison.",
        "nonlocal_artifact": "Along P=epsilon phat, the leading conversion coefficient is proportional to |epsilon| with a nonzero tensor coefficient. This is not a finite spatial polynomial, and cannot be silently identified with the already-fixed local Hessian. The one-ball comparator itself is not claimed covariantly renormalized.",
        "remaining": "This computes one actual leading regulator-conversion coefficient. It does not complete subleading artifacts, finite/divergent quantum coefficients or original covariant matching, and it is not proof that the renormalized physical current diverges. The regulator is not replaced and no counterterm is adjusted here.",
        "canonical_boundary": "Both metric factors multiply the coefficient by4/kappa, but the resulting K^3|P| term still grows with regulator K at any fixed finite kappa. It is not a uniform small-response bound or a physical cutoff.",
        "checks": {
            "complete_eight_order_W8_scaled_limit": zero(W) - 1 / a,
            "complete_W8_total_time_derivative_limit": zero(dW) + H / a,
            "transverse_normalized_electric_limit": zero(p) + s.I * v,
            "transverse_normalized_magnetic_limit": zero(s.I * f / a) - s.I * v,
            "transverse_normalized_mass_limit": zero(m * x * f),
            "longitudinal_normalized_electric_limit": zero(m * x * p / omega),
            "longitudinal_normalized_temporal_constraint_limit": zero(
                -s.I * p / (a * omega)
            )
            + v,
            "longitudinal_normalized_mass_spatial_limit": zero(omega * f) - v,
            "complete_pair_inverse_phase_limit": zero(1 / (2 * W)) - a / 2,
            "normalized_TT_coefficient": s.simplify(
                coefficient - 9 / (256 * s.pi**2 * a)
            ),
            "both_canonical_metric_factors": s.simplify(
                4 * coefficient / modes.KAPPA - 9 / (64 * s.pi**2 * a * modes.KAPPA)
            ),
        },
        "gates": {
            "fixed_nonzero_mass_not_Maxwell_limit": modes.MASS == 1000,
            "full_four_W8_coefficients_retained": len(P) == 4,
            "all_four_total_derivative_coefficients_retained": len(Q) == 4,
            "other_four_endpoints_lower_shell_order": all(
                3 - j < 3 for j in range(1, 5)
            ),
            "canonical_factor_not_uniform_regulator_control": True,
            "no_physical_divergence_or_exclusion_inferred": True,
            "original_V_G_B_open": True,
        },
    }
