"""Finite time-boundary extraction is not yet covariant spatial subtraction."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes

from . import bounds, ibp


@cache
def data():
    p, m = s.symbols("p mass", real=True, positive=True)
    kernel = 1 / (m + s.sqrt(m * m + p * p))
    z = s.Symbol("z", real=True)
    checks = {
        "exact_finite_and_unresolved_boundary_split": s.expand(
            ibp.boundary(6) - ibp.boundary(5) - ibp.G[0] * ibp.iterates()[5]
        ),
        "canonical_finite_reference_bound": 4 * bounds.DISPLAY / modes.KAPPA
        - 4 * s.Rational(1, 10) ** 752,
        "canonical_finite_reference_tail": 4 * bounds.TAIL / modes.KAPPA
        - 4 * s.Rational(1, 10) ** 748,
        "canonical_complete_known_finite_piece": 4 * (2 * 10**48) / modes.KAPPA
        - 8 * s.Rational(1, 10) ** 752,
        "canonical_complete_known_finite_tail": 4 * (2 * 10**52) / modes.KAPPA
        - 8 * s.Rational(1, 10) ** 748,
        "inverse_sum_phase_has_nonpolynomial_spatial_series": s.diff(
            1 / (m + s.sqrt(m * m + z * z)), z, 6
        ).subs(z, 0)
        + s.Rational(225, 8) / m**7,
        "unresolved_first_five_time_endpoints": len(ibp.iterates()) - 2 - 5,
    }
    return {
        "exact_common_regulator": "J_actual,K=R_state,K+C_ref,K+sum_(j=0)^4 B_j,K+F_K. Here R_state is the complete S197 memory/contact difference, C_ref is the FULL reference metric contact, B_j are the exact equal-time terms, and F=B5+bulk6.",
        "separate_projectors": "Every B_j and the bulk retains the common two-created-mode overlap. C_ref retains its original one-momentum band and full spatial convolution. These regions are not identified before cancellations or matching.",
        "finite_progress": "F has a complete all-internal/all-external momentum limit with source time derivatives only through6 and an explicit1/K regulator error. The actual-state remainder is already controlled separately by S197.",
        "complete_known_piece": "With M[D]^2=||D||L2^2+||grad D||L2^2 and N61[Gamma]^2=sum_(j=0)^6||partial_t^j Gamma||L2^2+||grad Gamma||L2^2, the complete known finite piece R_state+F is below2e48 M[D]N61[Gamma], with regulator error below2e52 M[D]N61[Gamma]/K. Both canonical displays are8e-752 and8e-748/K. This does not bound the remaining contact-plus-five-endpoint sector.",
        "canonical": "The canonical finite-reference display is4e-752 ||D||L2 S6[Gamma], with tail4e-748 ||D||L2 S6[Gamma]/K. This is not a bound on the unresolved contact-plus-five-endpoint sector.",
        "spatial_warning": "Equal-time source jets do not imply spatial locality: even the elementary frozen inverse-sum frequency1/(m+sqrt(m^2+p^2)) is not a finite polynomial in p. One must derive the complete high-internal-momentum spatial expansion and match the full tensor/contact combination to the original covariant prescription.",
        "reference_warning": "The constant-alpha W8 readout remains a comparison device, not a new exact physical state. No finite local counterterm is adjusted, no boundary coefficient is dropped, and no independent reference Ward identity is asserted.",
        "remaining": "Original fixed spatial diagonal matching, full mixed response and inverse, finite-amplitude interacting background, stability, physical cutoff and V/G/B remain OPEN.",
        "checks": checks,
        "gates": {
            "unresolved_sector_contains_actual_full_contact": True,
            "time_local_not_automatically_space_local": kernel != 0,
            "same_actual_mass_and_kappa": modes.MASS == 1000 and modes.KAPPA == 10**800,
            "no_odd_boundary_discard_without_full_tensor_contraction": True,
            "combined_known_finite_piece_display": bounds.DISPLAY + 2 * 10**23
            < 2 * 10**48,
            "combined_known_finite_tail_display": bounds.TAIL + 10**27 < 2 * 10**52,
        },
    }
