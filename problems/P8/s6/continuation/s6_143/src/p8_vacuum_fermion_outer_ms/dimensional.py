"""Dimensional on-shell insertion and spectral normalization before MS subtraction."""

from functools import cache

import sympy as sp


@cache
def data():
    e = sp.Symbol("epsilon", real=True)
    v, m, beta, u, s, Tad, I = sp.symbols("v m beta u s tadpole pole")
    x = (1 + beta * u) / 2
    cut_delta = sp.expand(v * x * (1 - x) - m * m)
    B = sp.Function("B")
    f = (4 * m * m - s) * B(s) - 2 * Tad
    exactR = f - f.subs(s, 1) - (s - 1) * sp.diff(f, s).subs(s, 1)
    proposed = (4 * m * m - s) * (B(s) - B(1)) - (s - 1) * (4 * m * m - 1) * sp.diff(
        B(s), s
    ).subs(s, 1)
    scalar_bubble = sp.exp(sp.EulerGamma * e) * sp.gamma(e)
    zero_mass_inverse = scalar_bubble * (4 - 2 / (e - 1))
    normal = sp.sqrt(sp.pi) / (2 * sp.gamma(sp.Rational(3, 2) - e))
    return {
        "dimensional_scalar_bubble": "e^(gamma epsilon) Gamma(epsilon) mu^(2epsilon)/Q integral_0^1 [m^2-x(1-x)s]^(-epsilon) dx",
        "fermion_inverse_over_2NY": "(4m^2-s) B_D(m,m;s)-2 Tad_D(m)",
        "exact_dimensional_OS_remainder": proposed,
        "spectral_density_over_C": "e^(gamma epsilon) mu^(2epsilon) 4^epsilon sqrt(pi)/(2Gamma(3/2-epsilon)) v^(1-epsilon) (1-4m^2/v)^(3/2-epsilon)",
        "regulated_F": "integral_(4m^2)^infinity rho_D(v)/(v-1)^2 B_D(v,1;0) dv",
        "common_dimensional_convention": "MS-bar loop normalization e^(gamma epsilon); mu=mF. Keep the exact D-dimensional physical mass/residue subtraction through the outer finite part.",
        "scope": "A dimensional-regulator representation, not a claim of a new physical noninteger-dimensional spectrum or a normalized exact Kallen-Lehmann measure.",
        "checks": {
            "cut_parameter_geometry": sp.expand(
                cut_delta - v * beta**2 * (1 - u * u) / 4
            ).subs(m * m, v * (1 - beta**2) / 4),
            "gamma_reflection_spectral_cancellation": sp.gammasimp(
                sp.gamma(e) * sp.gamma(1 - e) * sp.sin(sp.pi * e) / sp.pi
            )
            - 1,
            "four_dimensional_density_normalization": normal.subs(e, 0) - 1,
            "exact_OS_remainder_dictionary": sp.expand(exactR - proposed),
            "affine_UV_pole_removed_by_OS_subtraction": sp.expand(
                ((6 * m * m - s) * I) - ((6 * m * m - 1) * I) - (s - 1) * (-I)
            ),
            "one_loop_MS_bubble_finite_zero_at_mu_m": sp.limit(
                scalar_bubble - 1 / e, e, 0
            ),
            "same_first_mass_inverse_pole": sp.limit(e * zero_mass_inverse, e, 0) - 6,
            "same_first_mass_inverse_finite_value": sp.limit(
                zero_mass_inverse - 6 / e, e, 0
            )
            - 2,
            "threshold_beta_power_from_scalar_trace": (1 - 2 * e) + 2 - (3 - 2 * e),
            "two_loop_scale_factor_at_mu_m": sp.powsimp(4**e * 4 ** (-2 * e))
            - 4 ** (-e),
        },
    }
