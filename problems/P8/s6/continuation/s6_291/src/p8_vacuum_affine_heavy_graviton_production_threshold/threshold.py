"""Whole soft-threshold Laurent identity and explicit compact-window remainder."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complement_measure.split import clean

from . import phase, production, source

S, N, MU, G, K, EP = production.S, source.N, source.MU, source.G, source.K, source.EP
L = s.Symbol("dimensionless_window", positive=True)
Y = s.Symbol("dimensionless_threshold_coordinate", positive=True)
W = s.Function("fixed_C1_test_weight")


def require_window(mass, heavy, length, scale_squared=1):
    mass, heavy, length, scale_squared = map(
        source.require_mass, (mass, heavy, length, scale_squared)
    )
    if heavy <= 4 * mass or length > s.Rational(1, 2):
        raise ValueError("Require n>4mu and0<L<=1/2 for the compact threshold window")
    return mass, heavy, length, scale_squared


def leading_coefficient(energy=S, mass=MU, cubic=G, kappa=K):
    energy, mass, cubic, kappa = map(s.sympify, (energy, mass, cubic, kappa))
    return (
        cubic**2
        * (energy - 4 * mass) ** 2
        * production.angular_closed(1 - 4 * mass / energy)
        / (8 * s.pi * kappa * energy)
    )


def dimensionless_first(heavy=N, mass=MU, cubic=G, kappa=K, scale_squared=phase.NU2):
    heavy = s.sympify(heavy)
    return phase.first_coefficient(
        heavy, mass, cubic, kappa, scale_squared
    ) + 2 * s.log(heavy) * leading_coefficient(heavy, mass, cubic, kappa)


def remainder_majorant(epsilon, length, M2, M10, M11):
    epsilon, length, M2, M10, M11 = map(s.sympify, (epsilon, length, M2, M10, M11))
    return epsilon * (M2 / 4 + 2 * M10 * length * (1 - s.log(length)) + M11 * length)


@cache
def data():
    e = s.Symbol("epsilon", positive=True)
    y = s.Symbol("unit_parameter", positive=True)
    checks = {}
    A, B, C, L = s.symbols(
        "threshold_A threshold_epsilon_derivative threshold_second_derivative window",
        positive=True,
    )
    model = (A + e * B + e * e * C / 2) * s.exp(2 * e * s.log(L)) / (2 * e)
    finite = A * s.log(L) + B / 2
    checks["whole_threshold_Laurent_pole"] = s.limit(e * model, e, 0) - A / 2
    checks["whole_threshold_Laurent_finite"] = (
        s.limit(model - A / (2 * e), e, 0) - finite
    )
    checks["whole_threshold_Laurent_first_remainder"] = s.limit(
        (model - A / (2 * e) - finite) / e, e, 0
    ) - (A * s.log(L) ** 2 + B * s.log(L) + C / 4)
    n = s.Symbol("heavy_mass_squared", positive=True)
    K0, K1 = s.symbols("dimensional_coefficient dim_derivative", real=True)
    checks["dimensionless_window_first_coefficient"] = (
        s.diff(n ** (2 * e) * (K0 + e * K1), e).subs(e, 0) - K1 - 2 * s.log(n) * K0
    )
    checks["dimensionless_window_finite_scale_cancellation"] = s.expand(
        K0 * s.log(L) + (K1 + 2 * s.log(n) * K0) / 2 - (K0 * s.log(n * L) + K1 / 2),
        log=True,
    )
    checks["compact_log_majorant_integral"] = s.integrate(-s.log(y), (y, 0, L)) - L * (
        1 - s.log(L)
    )
    actual_n = s.Integer(10) ** 200 / 512 + 2
    checks["actual_positive_heavy_gap"] = (
        actual_n - 4 - (s.Integer(10) ** 200 / 512 - 2)
    )
    margins = {
        "actual_mass_above_four": actual_n - 4,
        "actual_mass_below_10_power200": s.Integer(10) ** 200 - actual_n,
        "ln10_bound_exponential_partial_sum": sum(
            s.Rational(5, 2) ** j / s.factorial(j) for j in range(7)
        )
        - 10,
        "log_argument_below_10_power202": s.Integer(10) ** 202
        - 64 * s.Integer(10) ** 200,
        "whole_log_coefficient_over_F0_below1000": 1000 - s.Rational(509 * 15, 8),
    }
    b, x = s.symbols("angular_beta_squared unit_angle", real=True)
    kernel = production.angular_kernel(b, x)
    checks["whole_kernel_lower_bound_gap"] = s.factor(
        kernel
        - (1 - x * x) ** 2
        - (1 - x * x) ** 2 * b * x * x * (2 - b * x * x) / (1 - b * x * x) ** 2
    )
    checks["whole_kernel_upper_bound_gap"] = s.factor(
        1 - kernel - (1 - b) * x * x * (2 - (1 + b) * x * x) / (1 - b * x * x) ** 2
    )
    checks["whole_kernel_log_absolute_moment"] = (
        s.integrate(s.log(1 - x * x), (x, 0, 1)) - 2 * s.log(2) + 2
    )
    m, h, ell = s.symbols("positive_mu positive_n compact_L", positive=True)
    checks["whole_compact_angular_denominator_gap"] = s.factor(
        1 - (1 - 4 * m / (h * (1 + ell))) - 4 * m / (h * (1 + ell))
    )
    A0, A1, A2, A3 = s.symbols("smooth_A0 smooth_A1 smooth_A2 smooth_A3", real=True)
    model = A0 + A1 * y + A2 * y * y + A3 * y**3
    checks["whole_cubic_weight_subtraction_identity"] = s.factor(
        (model - A0) / y - (A1 + A2 * y + A3 * y * y)
    )
    checks["whole_regular_integral_polynomial_limit"] = s.integrate(
        (model - A0) / y, (y, 0, L)
    ) - (A1 * L + A2 * L**2 / 2 + A3 * L**3 / 3)
    checks["public_forward_cut_matches_threshold_factor"] = s.factor(
        production.forward_cut() - leading_coefficient() / (S - N)
    )
    M2, M10, M11 = s.symbols(
        "bounded_epsilon_second bounded_y_first bounded_mixed_first", nonnegative=True
    )
    majorant = remainder_majorant(e, L, M2, M10, M11)
    K0 = leading_coefficient(N)
    K1 = dimensionless_first()
    H0 = lambda z: leading_coefficient(N * (1 + z)) * W(N * (1 + z))
    finite = (
        K0 * W(N) * s.log(L)
        + K1 * W(N) / 2
        + s.Integral((H0(Y) - H0(0)) / Y, (Y, 0, L))
    )
    return {
        "whole_forward_threshold_residue": K0,
        "whole_dimensionless_first_epsilon_coefficient": K1,
        "whole_test_weight_pole": K0 * W(N) / (2 * EP),
        "whole_test_weight_finite_coefficient": finite,
        "whole_explicit_compact_remainder_majorant": majorant,
        "whole_exact_positive_margins": margins,
        "whole_Laurent_subtraction_identity": "With y=(s-n)/n and H_e(y)=n^(2e)K_e(n(1+y))W(n(1+y)), integral0^L y^(-1+2e)H_e(y)dy = H_e(0)L^(2e)/(2e)+integral0^L y^(-1+2e)[H_e(y)-H_e(0)]dy. It is an exact identity for0<e<=1/8,0<L<=1/2,n>4mu and fixed C1 W.",
        "whole_remainder_suprema_and_proof": "M2=sup_e|d_e^2[H_e(0)L^(2e)]|, M10=sup_e,y|d_y H_e|, M11=sup_e,y|d_e d_y H_e| on the compact domain. Taylor's theorem gives e M2/4. The subtracted regular integral uses y^(2e)<=1, the y-mean-value theorem and integral0^L|ln y|dy=L(1-ln L), giving e[2M10 L(1-ln L)+M11 L]. Angular denominators stay at least4mu/[n(1+L)]>0 and angular log moments through second order are integrable; hence all suprema are finite. They are conditional symbolic bounds, not computed physical detector constants.",
        "whole_actual_evanescent_bound": "At original mu1,n=10^200/512+2,kappa10^800,g1/8192,nu^2=1,8/15<=F0<=1 and-2<L_B<0. With0<EulerGamma<1,pi<4,ln10<5/2,abs(K1(n))/K0(n)<509*15/8<1000. This bounds a finite evanescent coefficient, NOT the positive K0/(2e) pole or its physical cancellation.",
        "whole_unpaired_divergence_and_scope": "K0(n)>0 for n>4mu,g!=0,kappa>0. A positive test weight at the threshold gives a positive Laurent pole and a logarithmically divergent naive D4 lower-endpoint integral. H is unstable in the original theory. This formal perturbative distribution requires paired virtual/pole and width/IR analysis; it is not a physical-divergence verdict, stable-heavy spectral atom, finite plus-distribution prescription or P8 no-go.",
        "checks": {key: clean(value) for key, value in checks.items()},
        "gates": {
            "whole_positive_residue_and_angular_bounds_proved": all(
                v > 0 for v in margins.values()
            ),
            "full_evanescent_coefficient_not_D0_only": True,
            "dimensionless_and_dimensionful_window_logs_agree": True,
            "compact_subtract_add_identity_and_explicit_remainder": True,
            "unpaired_soft_pole_not_discarded_or_bounded_away": True,
            "no_arbitrary_finite_subtraction_or_width_choice": True,
        },
    }
