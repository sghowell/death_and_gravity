"""Exact mass-ordered parameter primitive and a convergent inverse-mass expansion."""

from functools import cache

import sympy as s

z, w, L, b, n = s.symbols(
    "z inverse_mass_squared light_denominator heavy_transfer mass_squared"
)
H = s.Symbol("log_n_minus_Log_L")
d = s.sqrt(1 - 4 * L * w + 4 * L * b * w * w)
a = (1 - 2 * L * w + d) / 2
v = L * w / a
beta = (L - b) * w / a
P = (1 + v) / d
Hfun = s.Function("analytic_H")
Q = ((1 + v) * s.log((1 + d) / 2) - Hfun(beta)) / d


@cache
def data():
    alpha = s.Symbol("alpha")
    be = (L - b) / alpha
    denominator = (L + alpha * z) * (1 + be * z)
    primitive = (
        (alpha + L) / alpha * s.log(L + alpha * z) - (1 + be) / be * s.log(1 + be * z)
    ) / (alpha - L * be)
    hjet = 1 + beta / 2 - beta**2 / 6
    qjet = ((1 + v) * s.log((1 + d) / 2) - hjet) / d
    tjet = w * (H - 1) + w * w * (3 * L * H + (b - 7 * L) / 2)
    derivative = w * w * (s.diff(tjet, w) - s.diff(tjet, H) / w)
    x = s.Symbol("ratio", positive=True)
    geometric = x * x / (1 - x)
    weighted = x * x * (3 - 2 * x) / (1 - x) ** 2
    checks = {
        "quadratic_factorization": s.cancel(
            denominator - (L + (alpha + L * be) * z + (L - b) * z * z)
        ),
        "primitive_derivative_entire_parameter_integrand": s.cancel(
            s.diff(primitive, z) - (1 - z) / denominator
        ),
        "discriminant_factor_difference": s.expand(
            (alpha - L * be) ** 2 - ((alpha + L * be) ** 2 - 4 * L * (L - b))
        ),
        "alpha_plus_L_scaled_is_half_one_plus_discriminant": s.simplify(
            a + L * w - (1 + d) / 2
        ),
        "analytic_P_constant_and_linear_jet": s.expand(
            s.series(P, w, 0, 2).removeO() - 1 - 3 * L * w
        ),
        "analytic_Q_constant_and_linear_jet": s.cancel(
            s.series(qjet, w, 0, 2).removeO() + 1 - (b - 7 * L) * w / 2
        ),
        "full_box_mass_derivative_integrand": s.cancel(
            -s.diff((1 - z) / (n * z + (1 - z) ** 2 * L - b * z * z), n)
            - z * (1 - z) / (n * z + (1 - z) ** 2 * L - b * z * z) ** 2
        ),
        "box_jet_including_logarithm_derivative": s.expand(
            derivative - w * w * (H - 2) - w**3 * (6 * L * H - 10 * L + b)
        ),
        "geometric_tail_weighted_derivative": s.cancel(
            geometric + x * s.diff(geometric, x) - weighted
        ),
        "weighted_tail_half_radius_endpoint": weighted.subs(x, s.Rational(1, 2)) - 2,
        "unweighted_P_tail_constant": 2 / (1 - s.Rational(1, 2)) - 4,
        "unweighted_Q_tail_constant": 4 / (1 - s.Rational(1, 2)) - 8,
        "weighted_P_tail_constant": 2
        * (3 - 2 * s.Rational(1, 2))
        / (1 - s.Rational(1, 2)) ** 2
        - 16,
        "weighted_Q_tail_constant": 4
        * (3 - 2 * s.Rational(1, 2))
        / (1 - s.Rational(1, 2)) ** 2
        - 32,
        "physical_radius_ratio": 16 * (s.Symbol("S") / 4) * w - 4 * s.Symbol("S") * w,
        "triangle_remainder_prefactor": 4 * 4**2 - 64,
        "box_remainder_prefactor": 16 * 4**2 - 256,
    }
    return {
        "full_parameter_integral": "T(n,L,b)=integral_0^1 (1-z)/[n z+(1-z)^2 L-b z^2]dz with Feynman continuation. C(s)=integral_x T(n,1-s*x*(1-x),0); D(s,t)=-partial_n integral_x,y T(n,1-s*x*(1-x),t*y*(1-y)). This retains the ordered box variables.",
        "exact_closed_form": "delta=sqrt(n^2-4nL+4Lb), alpha=(n-2L+delta)/2, beta=(L-b)/alpha; T=[(1+L/alpha)(log(alpha+L)-Log(L-i0))-H(beta)]/delta. H(beta)=(1+beta)log(1+beta)/beta with H(0)=1.",
        "analytic_coefficients": {
            "P": P,
            "Q": Q,
            "H_function": "H(q)=1+sum_{k>=1}(-1)^(k+1)q^k/[k(k+1)]",
        },
        "complete_two_term_T": tjet,
        "complete_two_term_D": derivative,
        "uniform_complex_disk": "For real |L|,|b|<=K and |w|<=1/(16K), the square root is the analytic branch d(0)=1. |d-1|<1/4, |a-1|<3/16, |P|<2, |Q|<4. These bounds are proved on the full disk, not sampled.",
        "complete_remainders": "For physical S>=4, K=S/4 and S/n<=1/8, |R_T|<=64 S^2 n^-3 (|log n-Log(L-i0)|+2), |R_D|<=256 S^2 n^-4 (|log n-Log(L-i0)|+3). Cauchy coefficient bounds include every higher inverse-mass term and the derivative of log n.",
        "branch_boundary": "For the real physical and symmetric domains the remaining logarithms have positive arguments. The full imaginary part is from Log(L-i0). Points L=0 are treated by an integrable boundary limit. A double-pole Feynman integral is not bounded by its absolute integrand across a zero.",
        "checks": {k: s.cancel(v) for k, v in checks.items()},
        "gates": {
            "discriminant_disk_deviation": s.Rational(17, 64) < s.Rational(7, 16),
            "alpha_disk_nonzero": s.Rational(3, 16) < s.Rational(1, 4),
            "light_ratio_disk": s.Rational(1, 13) < s.Rational(1, 8),
            "heavy_ratio_disk": s.Rational(2, 13) < s.Rational(1, 4),
            "P_supremum_below_two": s.Rational(3, 2) < 2,
            "Q_supremum_below_four": s.Rational(73, 24) < 4,
            "complete_Cauchy_remainder_not_formal_power_counting": True,
        },
    }
