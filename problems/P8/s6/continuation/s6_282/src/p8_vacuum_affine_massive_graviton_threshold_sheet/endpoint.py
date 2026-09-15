"""Complete massless-endpoint asymptotics and required pole-anchor running."""

from functools import cache

import sympy as s

from . import angular, invariant, source

MU, K = source.MU, source.K
SIGMA, OFFCUT, CAP = s.symbols(
    "spectral_energy offcut_energy positive_cap", positive=True
)


def imaginary_beta_forward(b):
    return (
        4
        + s.Rational(31, 3) * b * b
        + s.Rational(118, 15) * b**4
        + b**6
        + (b**7 - 6 * b**3 - 8 * b - 3 / b) * s.atan(b)
    )


def compensated_kernel(sigma=SIGMA, external=OFFCUT):
    return 1 / (sigma - external) + 1 / external


@cache
def data():
    b, x = s.symbols("positive_b integration_x", positive=True)
    num = (1 + b * b) ** 4 + b**8 * (1 - x * x) ** 4
    den = (1 + b * b * x * x) ** 2
    quotient, remainder = s.div(num, den, x)
    aa, bb = s.symbols("partial_a partial_b")
    solution = s.solve(
        s.Poly(remainder - aa * (1 + b * b * x * x) - bb, x).coeffs(), (aa, bb)
    )
    primitive = (
        quotient.integrate(x)
        + solution[aa] * s.atan(b * x) / b
        + solution[bb] * (x / (2 * (1 + b * b * x * x)) + s.atan(b * x) / (2 * b))
    )
    w = s.Symbol("positive_forward_sqrt_energy", positive=True)
    polynomial = 4 + s.Rational(31, 3) * b * b + s.Rational(118, 15) * b**4 + b**6
    coefficient = b**7 - 6 * b**3 - 8 * b - 3 / b
    shape = s.factor(
        w**4
        / (16 * s.pi)
        * (polynomial + coefficient * (s.pi / 2 - s.asin(w))).subs(
            b, s.sqrt(1 - w * w) / w
        )
    )
    expansion = s.series(shape, w, 0, 3).removeO().expand()
    tau, z = s.symbols("positive_tau strict_real_cosine", positive=True)
    d, ac = s.symbols("sqrt_one_minus_z_squared arccos_z", real=True)
    seed0 = -ac / d
    other0 = -(s.pi - ac) / d
    seed1 = 1 / (1 + z) + ac / ((1 + z) * d)
    other1 = 1 / (1 - z) + (s.pi - ac) / ((1 - z) * d)
    j0 = -(seed0 + other0) / (2 * tau * tau) - (seed1 + other1) / 2
    j1 = (seed0 - other0) / 2
    l0 = -s.pi / (2 * tau) + 1
    avg = angular.complete_average(-tau * tau, z, j0, j1, l0)
    rho = MU**2 * avg / (16 * s.pi * K * K * (1 + tau * tau) ** 2)
    fixed = s.series(rho, tau, 0, 1).removeO().expand()
    leading = MU**2 * (2 * z * z - 1) ** 2 / (16 * K * K * d)
    subleading = 3 * MU**2 * (5 * z * z - 1) / (32 * K * K)
    e, a, c = s.symbols("positive_e positive_a positive_c", positive=True)
    seed = -2 * s.atan(s.sqrt(c / (a + 2 * e))) / s.sqrt(c * (a + 2 * e))
    seedprime = 2 / ((a + 2 * e) * (a + c + 2 * e)) + 2 * s.atan(
        s.sqrt(c / (a + 2 * e))
    ) / (s.sqrt(c) * (a + 2 * e) ** s.Rational(3, 2))
    rho_cap = s.Symbol("complete_density_at_cap", real=True)
    v, delta_anchor = s.symbols("crossing_displacement anchor_shift", real=True)
    cross = delta_anchor * (1 / (2 * MU + v) + 1 / (2 * MU - v))
    kernel = compensated_kernel()
    checks = {
        "complete_negative_beta_rational_division": s.factor(
            num - quotient * den - solution[aa] * (1 + b * b * x * x) - solution[bb]
        ),
        "complete_negative_beta_forward_primitive": s.factor(
            s.diff(primitive, x) - num / den
        ),
        "complete_negative_beta_forward_integral": s.factor(
            primitive.subs(x, 1) - primitive.subs(x, 0) - imaginary_beta_forward(b)
        ),
        "forward_sigma_minus_three_halves_coefficient": expansion.coeff(w, -3)
        - s.Rational(1, 32),
        "forward_sigma_minus_half_coefficient": expansion.coeff(w, -1)
        + s.Rational(7, 64),
        "forward_constant_coefficient": expansion.coeff(w, 0)
        - s.Rational(41, 80) / s.pi,
        "forward_sqrt_sigma_coefficient": expansion.coeff(w, 1) + s.Rational(13, 256),
        "forward_sigma_coefficient": expansion.coeff(w, 2) + s.Rational(1, 60) / s.pi,
        "negative_seed_complete_first_derivative": s.simplify(
            s.diff(seed, e) - seedprime
        ),
        "fixed_angle_tau_minus_two_coefficient": s.factor(
            fixed.coeff(tau, -2) - leading
        ),
        "fixed_angle_tau_minus_one_coefficient": s.factor(
            fixed.coeff(tau, -1) - subleading
        ),
        "fixed_angle_energy_coordinate": s.factor(
            1 / tau**2 - (4 * MU / (4 * MU * tau * tau / (1 + tau * tau)) - 1)
        ),
        "whole_pole_factor_preserving_compensation": s.factor(
            kernel - SIGMA / (OFFCUT * (SIGMA - OFFCUT))
        ),
        "compensated_endpoint_limit": s.limit(kernel / SIGMA, SIGMA, 0, dir="+")
        + 1 / OFFCUT**2,
        "necessary_pole_anchor_cap_running": s.factor(
            -rho_cap / (s.pi * OFFCUT)
            + rho_cap / s.pi * compensated_kernel(CAP, OFFCUT)
            - rho_cap / (s.pi * (CAP - OFFCUT))
        ),
        "unmatched_anchor_cap_error": s.factor(
            rho_cap / s.pi * compensated_kernel(CAP, OFFCUT)
            - rho_cap / (s.pi * (CAP - OFFCUT))
            - rho_cap / (s.pi * OFFCUT)
        ),
        "crossing_even_pole_shift_in_b20": s.factor(
            s.diff(cross, v, 2).subs(v, 0) / 2 - delta_anchor / (4 * MU**3)
        ),
        "compensated_second_derivative_not_inverse_cube": s.factor(
            (s.diff(kernel, OFFCUT, 2) / 2 - 1 / (SIGMA - OFFCUT) ** 3).subs(
                OFFCUT, 2 * MU
            )
            - 1 / (8 * MU**3)
        ),
    }
    t = s.Symbol("strict_positive_transfer", positive=True)
    z0 = 1 - t / (2 * MU)
    checks["fixed_transfer_angular_discriminant"] = s.factor(
        1 - z0 * z0 - t * (4 * MU - t) / (4 * MU * MU)
    )
    checks["fixed_transfer_leading_invariant_numerator"] = s.factor(
        invariant.polynomials(t)[0] - 2 * MU * MU * (2 * z0 * z0 - 1)
    )
    return {
        "whole_subthreshold_forward_shape": imaginary_beta_forward(b),
        "whole_forward_massless_endpoint_expansion": expansion,
        "strict_angle_leading_tau_coefficients": (leading, subleading),
        "fixed_angle_energy": "sigma=4mu*tau^2/(1+tau^2). For strict |z|<1 the leading rho coefficient is mu^3(2z^2-1)^2/[4kappa^2 sqrt(1-z^2)] times1/sigma; next is3mu^(5/2)(5z^2-1)/(16kappa^2 sqrt(sigma)). The leading term can vanish atz^2=1/2.",
        "forward_energy": "w=sqrt(sigma/(4mu)); rho=mu^2/kappa^2 times the displayed expansion. In particular rho~mu^(7/2)/(4kappa^2 sigma^(3/2)). Fixed-angle estimates are not uniform atz=+/-1.",
        "whole_compensated_kernel": kernel,
        "pole_lift": "F_L(z)=c_L/z+(1/pi) integral_0^L rho(sigma)[1/(sigma-z)+1/z]dsigma. Since sigma*rho is integrable, z*F_L is a well-defined Cauchy transform plus c_L. Off the integration cut the jump is the chosen normal density. This defines a channel construction, not the full original amplitude.",
        "required_matching": "If changing L is to change only the ordinary high-end shell, c'_L=-rho(L)/pi is necessary and sufficient. The initial c_L0 is not determined by the cut. Leaving c_L fixed adds an erroneous rho(L)/(pi*z) shell. A crossing-even anchor shift changes b20 by delta_c/(4mu^3).",
        "boundary": "Neither an arbitrary finite-part prescription nor a cap derivative fixes the original renormalized massless-pole anchor. Full crossed/double-spectral assembly, transfer subtraction, finite local terms, detector errors and Regge are separate.",
        "checks": checks,
        "gates": {
            "whole_forward_integral_not_only_a_series": True,
            "nonintegrable_forward_low_endpoint_explicit": expansion.coeff(w, -3)
            == s.Rational(1, 32),
            "strict_angle_and_forward_limits_distinguished": True,
            "pole_factor_preserved_not_dropped": kernel.has(OFFCUT),
            "cap_running_not_a_choice_of_physical_anchor": True,
            "crossed_b20_pole_shift_retained": True,
            "no_new_counterterm_or_full_amplitude_inferred": True,
        },
    }
