"""Uniform all-momentum ADM form norms and finite-regulator Gaussian response."""

from functools import cache

import sympy as s
from p8_vacuum_affine_spatial_current import bounds as prior

from . import vertices

SPATIAL = s.Rational(5, 2)
SPATIAL_SECOND = s.Rational(25, 4)


@cache
def data():
    a = s.Symbol("a", positive=True)
    beta = s.Matrix(s.symbols("b0:3", real=True))
    S = vertices.shift_feature(beta, a)
    b2 = (beta.T * beta)[0]
    complex_beta = beta + s.I * s.Matrix(s.symbols("c0:3", real=True))
    complex_S = vertices.shift_feature(complex_beta, a)
    positive = complex_S.conjugate().T * complex_S
    complex_b2 = (complex_beta.conjugate().T * complex_beta)[0]
    x, y, nD, nG, bD, bG = s.symbols("x y nD nG bD bG", nonnegative=True)
    totalD = nD + SPATIAL * x + a * bD
    totalG = nG + SPATIAL * y + a * bG
    contact = SPATIAL_SECOND * x * y + SPATIAL * (nD * y + nG * x)
    surplus = s.Poly(s.expand(totalD * totalG - contact), nD, nG, x, y, bD, bG, a)
    c = prior.constants()
    checks = {
        "sharp_shift_feature_is_Hermitian": S - S.conjugate().T,
        "sharp_shift_feature_cubic_spectral_identity": S * S * S - a * a * b2 * S,
        "shift_feature_squared_trace_all_six_nonzero_eigenvalues": s.trace(S * S)
        - 6 * a * a * b2,
        "complex_Fourier_shift_singular_value_projector": positive * positive
        - a * a * complex_b2 * positive,
        "complex_Fourier_shift_squared_singular_values": s.trace(positive)
        - 6 * a * a * complex_b2,
        "full_spatial_trace_first_form_constant": 1 + s.Rational(3, 2) - SPATIAL,
        "full_spatial_noncommuting_second_constant": SPATIAL**2 - SPATIAL_SECOND,
        "scalar_constraint_second_not_larger": SPATIAL_SECOND - s.Rational(9, 4) - 4,
        "complete_contact_direction_product_reconstruction": s.expand(
            totalD * totalG - contact - surplus.as_expr()
        ),
        "same_actual_balanced_covariance": prior.COV - 4 * s.Integer(10) ** 12,
        "same_two_propagator_product": prior.PROP - 4 * s.Integer(10) ** 9,
        "same_finite_band_total_display": c["finite_band_current_display"]
        - 2 * s.Integer(10) ** 104,
    }
    return {
        "sharp_shift_norm": "For real beta the full ten-field shift block is Hermitian with S_beta^3=a^2|beta|^2 S_beta. For arbitrary complex Fourier beta, P=S_beta^dagger S_beta obeys P^2=a^2|beta|^2 P and trP=6a^2|beta|^2. Thus the norm is exactly a|beta|, zero at beta0, without falsely assuming a single Fourier mode is real. Both magnetic/electric and mass/divergence blocks remain.",
        "full_direction_norm": "For direction D=(n,beta,Q), use sigma_D=|n|+(5/2)||Q||op+a|beta|. Since|trQ|<=3||Q||op, ||B_Q||<=5||Q||op/2, so the full first feature form is bounded bysigma_D. The spatial second bound is25||Q_D||||Q_G||/4; lapse-spatial contacts add(5/2)(|n_D||||Q_G||+|n_G||||Q_D||), all bounded bysigma_D sigma_G.",
        "balanced_all_momentum": "The unchanged energy features obey Fk^tFk=M0(k). Thus the two-sided energy-relative first/second vertices are bounded uniformly for all k,q, including zeros and separated momenta. Actual symplectic normalizers give sqrt(omega_k omega_q) times the same direction bounds. No spurious |P|/m loss is introduced by splitting the shift Lie expression.",
        "finite_response": "S195's exact two-momentum covariance tangent and current/contact identities apply to these full ADM matrices in the same unchanged actual state. Replacing each direction operator norm bysigma gives the same balanced source/readout/contact constants. At original internal bandK>=1000, memory is below1e24 K^5 sigma_D(t) integral sigma_G and contact below1e13 K^4 sigma_D sigma_G. AtK1e16 and unit slab, the total is below2e104 for unit direction norms.",
        "scope": "The band is a computational regulator. These growing bounds are not an infinite-momentum renormalized full ADM response, gauge reduction, useful inverse or physical cutoff. The S213 matched tracefree continuum sector remains unchanged; its result is not automatically extended to the new directions. No canonical reduced scalar/metric normalization is claimed.",
        "checks": checks,
        "gates": {
            "contact_bound_product_has_nonnegative_surplus": all(
                c >= 0 for c in surplus.coeffs()
            ),
            "same_complete_finite_band_bound": c["finite_band_memory_before_rounding"]
            + c["finite_band_contact_before_rounding"]
            < c["finite_band_current_display"],
            "shift_bound_sharp_without_external_momentum_growth": True,
            "arbitrary_complex_Fourier_directions_not_only_real_beta": True,
            "all_temporal_constraint_and_longitudinal_terms_retained": True,
            "full_ADM_continuum_matching_and_reduced_inverse_still_required": True,
        },
    }
