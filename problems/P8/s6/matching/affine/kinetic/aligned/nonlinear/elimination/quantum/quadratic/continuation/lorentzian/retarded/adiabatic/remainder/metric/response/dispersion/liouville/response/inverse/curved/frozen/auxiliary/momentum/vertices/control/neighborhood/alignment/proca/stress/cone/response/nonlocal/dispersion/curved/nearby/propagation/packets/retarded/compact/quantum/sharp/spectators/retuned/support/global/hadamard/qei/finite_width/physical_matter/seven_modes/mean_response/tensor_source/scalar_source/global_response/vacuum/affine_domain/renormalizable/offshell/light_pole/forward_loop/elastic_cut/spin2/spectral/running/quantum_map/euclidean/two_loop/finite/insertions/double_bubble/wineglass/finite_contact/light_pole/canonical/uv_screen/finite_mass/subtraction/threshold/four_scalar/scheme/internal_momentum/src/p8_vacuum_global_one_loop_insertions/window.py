"""Finite reference-momentum inverse control for both one-loop sectors."""

from functools import cache

import sympy as sp

from . import halfplane


def defect_bound(scalar_slope_upper, fermion_factor_upper, log_upper, kinetic_lower):
    r, C, ell, k = map(
        halfplane.exact,
        (scalar_slope_upper, fermion_factor_upper, log_upper, kinetic_lower),
    )
    if r < 0 or C < 0 or ell < 0 or k <= 0:
        raise ValueError(
            "Require nonnegative r,C,log bound and positive kinetic lower bound"
        )
    epsilon = (r + C * ell) / k
    return {
        "fractional_inverse_defect_upper": epsilon,
        "normalized_inverse_lower": 1 - epsilon,
        "normalized_inverse_upper": sp.Integer(1),
        "strictly_positive_on_declared_log_window": bool(epsilon < 1),
        "scope": "Selected complete one-loop Euclidean inverse only, assuming the declared finite log window; not a derived Wilsonian cutoff or all-order bound.",
    }


@cache
def data():
    M, g, Q, d, x = sp.symbols("M g Q d x", positive=True)
    m = sp.symbols("m", positive=True)
    A = x * (1 - x)
    Delta1 = x * M + (1 - x) ** 2
    b = A / Delta1
    h = d * b
    Rscalar = g / Q * (h - sp.log(1 + h))
    r, fp, fR, PiR, t = sp.symbols("r fp fR PiR t")
    kappa = 1 + r - fp
    inverse = 1 + t + (fR - PiR) / kappa
    z = sp.symbols("z", positive=True)
    K = 4 * m * m - 1
    reference_gap = m * m - (1 + (m**4 + 1) / K)
    checks = {
        "scalar_OS_parameter_gap": sp.expand(x * M + 1 - x - A - Delta1),
        "scalar_spacelike_gap_factorization": sp.factor(
            Delta1 + A * d - Delta1 * (1 + h)
        ),
        "scalar_remainder_zero_value": Rscalar.subs(d, 0),
        "scalar_remainder_zero_slope": sp.diff(Rscalar, d).subs(d, 0),
        "scalar_remainder_positive_second_derivative": sp.factor(
            sp.diff(Rscalar, d, 2) - g * b * b / (Q * (1 + h) ** 2)
        ),
        "scalar_remainder_below_slope_times_d": sp.factor(
            g * d * b / Q - Rscalar - g * sp.log(1 + h) / Q
        ),
        "scalar_OS_slope_dictionary": sp.factor(
            sp.diff(-g * sp.log(Delta1 - A * z) / Q, z).subs(z, 0) - g * b / Q
        ),
        "both_sector_defects_have_same_inverse_sign": sp.factor(
            (1 + t - inverse) * kappa - PiR + fR
        ),
        "kinetic_normalization_both_slopes": kappa - 1 - r + fp,
        "reference_window_log_argument_gap": sp.factor(
            reference_gap - (3 * m**4 - 5 * m * m) / K
        ),
        "actual_reference_energy_squared_dictionary": (10**400) ** 2 - (10**200) ** 4,
    }
    return {
        "scalar_spacelike_OS_parameter_kernel": Rscalar,
        "scalar_slope_integrand": g * b / Q,
        "scalar_nonnegative_remainder_envelope": "0<=Pi_scalar,R(-t)<=(t+1) r, r=Pi_scalar'(1).",
        "canonical_kinetic_normalization": kappa,
        "complete_one_loop_Euclidean_inverse": inverse,
        "relative_defect_envelope": "0<=(Pi_scalar,R(-t)-f_fermion,R(-t))/(kappa*(t+1)) <= [r+C log(1+(t+1)/K)]/kappa.",
        "reference_window": "0<=t<=Lambda_ref^2, Lambda_ref=10^400, mF=10^200. This is a specified test range, not an established Wilsonian cutoff.",
        "log_window_proof": "K=4mF^2-1 and 1+(Lambda_ref^2+1)/K<mF^2, since 3mF^4>5mF^2. Thus log(1+(t+1)/K)<log(mF^2)=400 log(10).",
        "scope": "The selected combined one-loop inverse has no extra Euclidean zero on this finite reference window when the displayed defect is below one. This establishes neither a full Lorentzian spectrum nor reflection positivity.",
        "checks": checks,
    }
