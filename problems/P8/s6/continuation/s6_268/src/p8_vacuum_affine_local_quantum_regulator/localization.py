"""Strict coherent effects, exact dimension-dependent tails and state estimates."""

from functools import cache

import sympy as s


def positive_integer(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer)) or value < 1:
        raise ValueError("Require a finite positive integer configuration dimension")
    return int(value)


def tail_polynomial(dimension, variable):
    dimension = positive_integer(dimension)
    return s.Poly(sum(variable**k / s.factorial(k) for k in range(dimension)), variable)


def ball_eigenvalue(dimension, occupation, variable):
    dimension = positive_integer(dimension)
    if (
        isinstance(occupation, bool)
        or not isinstance(occupation, (int, s.Integer))
        or occupation < 0
    ):
        raise ValueError("Require a nonnegative integer total occupation")
    return (
        1
        - s.exp(-variable)
        * tail_polynomial(dimension + int(occupation), variable).as_expr()
    )


@cache
def data():
    t = s.Symbol("half_squared_whitened_phase_radius", real=True)
    checks = {}
    tails = {}
    for d in (1, 2, 8, 16, 256):
        P = tail_polynomial(d, t)
        checks["whole_exact_gamma_tail_derivative_dimension_" + str(d)] = (
            P.diff() - P + s.Poly(t ** (d - 1) / s.factorial(d - 1), t)
        ).as_expr()
        tails[str(d)] = s.exp(-t) * P.as_expr()
    for d in (1, 2, 8):
        for m in range(4):
            P, next_P = tail_polynomial(d + m, t), tail_polynomial(d + m + 1, t)
            checks[
                "strict_total_occupation_eigenvalue_gap_" + str(d) + "_" + str(m)
            ] = (next_P - P - s.Poly(t ** (d + m) / s.factorial(d + m), t)).as_expr()
    # Exact Husimi width: vacuum Wigner covariance I/2, Husimi covariance I.
    z = s.Symbol("one_whitened_coordinate", real=True)
    alpha = s.Symbol("positive_symbol_decay", positive=True)
    integral = s.integrate(
        s.exp(-z * z / 2 - alpha * z * z) / s.sqrt(2 * s.pi), (z, -s.oo, s.oo)
    )
    checks["whole_Husimi_width_not_Wigner_width"] = s.simplify(
        integral - (1 + 2 * alpha) ** -s.Rational(1, 2)
    )
    for d in (1, 2, 8):
        # Radial phase volume in normalized measure, x=R^2/2.
        checks["whole_normalized_phase_ball_trace_" + str(d)] = s.expand(
            s.pi**d * (2 * t) ** d / (s.factorial(d) * (2 * s.pi) ** d)
            - t**d / s.factorial(d)
        )
    M, B, width, rho = s.symbols(
        "difference_bound interaction_bound positive_time initial_sqrt_tail",
        nonnegative=True,
    )
    duhamel = M * rho * width + M * B * width**2 / 2
    return {
        "whole_phase_effect": "E_chi=Q_V(chi) for smooth compact 0<=chi<=1 equal1 on a nonempty core. E_chi is a positive trace-class operator, tr(E_chi)=(2pi)^(-d) integral chi. For a Gaussian window, every nonzero vector has a Bargmann transform that cannot vanish on an open phase set. Both E_chi and I-E_chi have strictly positive quadratic forms. Compactness gives 0<||E_chi||<1 and ||I-E_chi||=1. E_chi is not a projection and no nonzero exact chart-supported state exists in this coherent-effect sense.",
        "whole_whitening_convention": "For the full pure unit-CCR covariance V, choose symplectic S with V=SS^T/2. With z=S w the same seed Husimi density has covariance I in 2d real phase dimensions. The original compact translation group preserves this norm; an invariant symplectic square root can be chosen. R below is the WHITENED phase radius, not the raw canonical or S266 invariant radius.",
        "whole_sharp_ball_spectrum": "For chi=1_{||w||<=R}, t=R^2/2, E_chi is diagonal in the squeezed total-number decomposition. At total occupation m, lambda_m=1-exp(-t) sum_{k=0}^{d+m-1} t^k/k!. Its multiplicity is binomial(d+m-1,m), lambda_0 is the largest eigenvalue, and every eigenvalue lies strictly between0 and1 for finite R>0. The sharp ball is a spectral diagnostic, not the smooth Hamiltonian cutoff.",
        "whole_exact_seed_leakage_tails": tails,
        "whole_tail_variable": t,
        "whole_smooth_cutoff_bracket": "If chi=1 on radius R1 and chi=0 outside radius R2, then lambda_0(R1)<=<psi0,E_chi psi0><=lambda_0(R2). For radial chi, psi0 remains an eigenvector. Outside the core, the initial Husimi mass is the displayed t_d(R1), not the narrower Wigner tail and never exactly zero.",
        "whole_dimension_sensitive_bounds": "For t>0, 1-t_d(R)<=t^d/d!, hence t_d(R)>=1-t^d/d!. For t>d, t_d(R)<=exp(-t+d+d*log(t/d)) by the gamma Chernoff bound. No fixed small domain is automatically probable as the mode dimension increases, and no uniform IR/UV regulator removal follows.",
        "whole_state_specific_symbol_difference": "If bounded delta_a vanishes on the whitened core radius R and |delta_a|<=M, the isometric compression inequality gives ||Q_V(delta_a)psi0||^2<=<psi0,Q_V(|delta_a|^2)psi0><=M^2 t_d(R). This is an operator-on-state norm bound, not an assertion of classical joint q,p support. For first-Weyl-calibrated cutoff choices, delta_a must include ALL cutoff derivative contacts; it vanishes only where both cutoffs identically equal1.",
        "whole_interacting_Duhamel_comparison": "For two bounded interaction generators A,B with delta=A-B=Q_V(delta_a), ||delta||<=M, ||B(u)||<=B0, and ||delta(u)psi0||<=M sqrt(tail) uniformly, their propagators from0 satisfy ||(U_A(t)-U_B(t))psi0||<=M sqrt(tail)|t|+M B0 |t|^2/2. This follows from exact Duhamel and ||U_B(s)psi0-psi0||<=B0|s|. It requires both bounded extensions on the same Hilbert space and does not compare with an undefined unregularized singular Hamiltonian. Initial leakage alone is not a later-time support theorem.",
        "whole_Duhamel_bound_polynomial": duhamel,
        "checks": checks,
        "gates": {
            "tail_dimension_includes_all_eight_physical_channels": True,
            "exact_gamma_tails_keep_every_term_through_dimension256": tail_polynomial(
                256, t
            ).degree()
            == 255,
            "full_Husimi_width_is_twice_Wigner_covariance": True,
            "compact_phase_effect_is_not_exact_state_projection": True,
            "initial_tail_not_assumed_to_control_full_interacting_time": True,
            "no_unevaluated_radius_promoted_to_small_probability": True,
        },
    }
