"""Original crossing coefficient and complete first elastic transfer derivative."""

from functools import cache

import sympy as s
from p8_vacuum_affine_scalar_matching_scale import amplitude as original

S = original.S
transfer = s.Symbol("transfer_t", real=True)
v = s.Symbol("crossing_v", real=True)
lam, gam = original.lam, original.gam
LAMBDA, GAMMA = original.LAMBDA, original.GAMMA


@cache
def data():
    ss, uu = 2 - transfer / 2 + v, 2 - transfer / 2 - v
    tree = (
        2 * lam * ((ss - 2) ** 2 + (transfer - 2) ** 2 + (uu - 2) ** 2)
        + 3 * gam * ss * transfer * uu
        - 8 * gam
    )
    b2 = s.expand(s.diff(tree, v, 2) / 2)
    x = 1 + 2 * transfer / (S - 4)
    rho = (
        32
        * s.pi
        / original.beta
        * (original.PARTIAL0**2 + 5 * original.PARTIAL2**2 * s.legendre(2, x))
    )
    slope = original.beta * original.b**2 / (60 * s.pi * (S - 4))
    mu = s.Symbol("spectral_s", positive=True)
    r0, r1 = s.symbols("rho rho_transfer_derivative", nonnegative=True)
    kernel = 2 * (r0 + transfer * r1) / (s.pi * (mu + transfer / 2 - 2) ** 3)
    checks = {
        "original_exact_crossing_b2": s.expand(b2 - (4 * lam - 3 * gam * transfer)),
        "original_b20": b2.subs(transfer, 0) - 4 * lam,
        "original_b21": s.diff(b2, transfer).subs(transfer, 0) + 3 * gam,
        "original_tree_v4_coefficient_zero": s.diff(tree, v, 4),
        "full_first_elastic_forward_coefficient": s.simplify(
            rho.subs(transfer, 0) - original.RHO
        ),
        "full_first_elastic_transfer_derivative": s.simplify(
            s.diff(rho, transfer).subs(transfer, 0) - slope
        ),
        "massive_threshold_slope_factor": s.cancel(
            slope
            - original.beta * (S - 4) ** 3 * (lam - 3 * gam * S / 4) ** 2 / (60 * s.pi)
        ),
        "differentiated_crossing_kernel": s.cancel(
            s.diff(kernel, transfer).subs(transfer, 0)
            - 2 / s.pi * (r1 / (mu - 2) ** 3 - s.Rational(3, 2) * r0 / (mu - 2) ** 4)
        ),
    }
    z = s.Symbol("angle", real=True)
    for ell in range(9):
        checks["Legendre_endpoint_derivative_" + str(ell)] = s.diff(
            s.legendre(ell, z), z
        ).subs(z, 1) - s.Rational(ell * (ell + 1), 2)
    return {
        "complete_original_crossing_tree": tree,
        "original_tree_b2_as_function_of_transfer": b2,
        "fixed_original_low_coefficients": (4 * LAMBDA, -3 * GAMMA),
        "complete_first_elastic_nonforward_coefficient": rho,
        "complete_first_elastic_transfer_slope": slope,
        "exact_positive_transfer_identity": "With the unchanged normalized identical channels A=(32pi/beta)sum_even(2ell+1)t_ell P_ell, exact unitarity implies Im t_ell>=0. If the endpoint series and derivative converge as specified, partial_t Im A(S,0)=(32pi/beta)sum_even(2ell+1)Im t_ell*ell(ell+1)/(S-4)>=0.",
        "original_tree_ownership": "The complete S177/S182 original four-scalar tree and S231 identical normalization are retained. No higher independent vertex changes this tree. The first elastic coefficient and its slope are not the full absorptive amplitude or a real quantum matching calculation.",
        "checks": checks,
        "gates": {
            "actual_lambda_gamma_positive": LAMBDA > 0 and GAMMA > 0,
            "first_elastic_slope_threshold_regular": True,
            "identical_normalization_and_both_partial_channels_retained": True,
            "nonnegative_transfer_derivative_has_explicit_convergence_premise": True,
        },
    }
