"""All-external-momentum finite-band response bounds and their boundary."""

from functools import cache
from math import factorial

import sympy as s
from p8_vacuum_affine_cd_source_noise import modes
from p8_vacuum_flat_dirac_hadamard.symbols import rational

COV = s.Integer(4) * 10**12
PROP = s.Integer(4) * 10**9
PAIR = s.Integer(10) ** 23
READOUT = s.Integer(10) ** 24
CONTACT = s.Integer(3) * 10**13
BAND = s.Integer(10) ** 16


def require_band(K):
    K = rational(K)
    if K < 1000:
        raise ValueError("Require a finite rational computational band K>=1000")
    return K


@cache
def constants():
    e_upper = sum(s.Rational(1, factorial(n)) for n in range(9)) + s.Rational(
        1, factorial(9)
    ) / (1 - s.Rational(1, 10))
    pair = 2 * COV * PROP * 2
    readout = 3 * 2 * PAIR
    contact = 3 * 2 * COV
    low_memory = s.Rational(4, 54) * READOUT * BAND**5
    low_contact = s.Rational(2, 54) * CONTACT * BAND**4
    return {
        "e_series_upper": e_upper,
        "full_covariance_pair_source_and_transport": pair,
        "full_six_dimensional_pair_readout": readout,
        "full_same_time_metric_contact": contact,
        "finite_band_memory_before_rounding": low_memory,
        "finite_band_contact_before_rounding": low_contact,
        "finite_band_current_display": 2 * s.Integer(10) ** 104,
        "canonical_projected_current_display": 8 * s.Rational(1, 10) ** 696,
    }


@cache
def data():
    c = constants()
    K, x = s.symbols("K x", positive=True)
    D = s.diag(1, -1, 0)
    mass = s.Integer(1000)
    zero_mode_cov = s.diag(s.eye(3) / (2 * mass), mass * s.eye(3) / 2)
    # Projected cos(p x) has zero first vertex but cos^2 has constant1/2.
    zero_contact = s.diag(mass**2 * D**2 / 2, D**2 / 2)
    actual = -s.trace(zero_contact * zero_mode_cov) / 2
    checks = {
        "all_band_radial_measure": s.integrate(x * x, (x, 0, K)) - K**3 / 3,
        "zero_mode_high_transfer_contact": actual + mass * s.trace(D * D) / 4,
        "proper_canonical_two_chain_factors": 4
        * c["finite_band_current_display"]
        / modes.KAPPA
        - c["canonical_projected_current_display"],
        "selected_band_is_computational_only": require_band(BAND) - 10**16,
    }
    return {
        "preparation_and_norm": "Use the actual unchanged isotropic CD state, uniform balanced covariance bound4e12, two-propagator product norm<=exp(22|t-s|)<4e9, and omega_k(t)<2nu_k on the unit slab, nu_k=sqrt(m^2+|k|^2/a(t0)^2). No growth estimate of exp(|k||t-s|) is used.",
        "pair_covariance_bound": "||delta Sigma_kq(t)||op <1e23 sqrt(nu_k nu_q) integral_(t0)^t ||Gamma_(k-q)(s)||op ds for the balanced off-diagonal covariance tangent. This is an operator norm, not a momentum-integrable ultraviolet majorant.",
        "ordered_pair_current_bound": "The propagated integrand has absolute value <1e24 nu_k nu_q ||D_(q-k)(t)||op integral ||Gamma_(k-q)(s)||op ds. The local contact integrand is <3e13 nu_k ||D||op||Gamma||op for a specified Fourier pair.",
        "finite_band": "For the common sharp Galerkin band |k|,|q|<=K with K>=1000 and q=k-p, the propagation overlap measure is at most K^3/(6pi^2)<K^3/54 and nu_k,nu_q<2K. The current's memory term is <1e24 K^5 integral||Gamma_p|| and its contact is <1e13 K^4||Gamma_p(t)||, uniformly in external p.",
        "selected_display": "At K=1e16, the full projected current is bounded by2e104 sup_I||Gamma_p|| for unit detector polarization. Multiplication by4/(a^3 kappa) yields the projected canonical linear-response display8e-696.",
        "large_transfer_contact": "If |p|>2K the propagation overlap is empty, but the second metric contact need not vanish: the pointwise product of opposite external Fourier modes has zero total momentum. The explicit zero-mode cos(p x) fixture has current derivative -m tr(D^2)/4, not zero.",
        "UV_boundary": "This is an exact finite-regulator response and an all-momentum vertex estimate, not a covariantly renormalized full spatial response. The displayed internal-momentum majorants grow; neither they nor S194's homogeneous subtraction are an infinite spatial-tail proof. The band is not a physical cutoff or a stability/inverse theorem.",
        "constants": c,
        "checks": checks,
        "gates": {
            "elementary_e_bound": c["e_series_upper"] < s.Rational(68, 25),
            "common_exponential_transport_bound": s.Rational(68, 25) ** 22 < PROP,
            "actual_initial_covariance_to_slab_bound": 108 * PROP < COV,
            "two_source_terms_and_both_propagators": c[
                "full_covariance_pair_source_and_transport"
            ]
            < PAIR,
            "trace_dimension_and_frequency_factors": c[
                "full_six_dimensional_pair_readout"
            ]
            < READOUT,
            "complete_contact_display": c["full_same_time_metric_contact"] < CONTACT,
            "finite_memory_display": c["finite_band_memory_before_rounding"] < 10**104,
            "finite_contact_display": c["finite_band_contact_before_rounding"] < 10**77,
            "full_finite_sum_display": c["finite_band_memory_before_rounding"]
            + c["finite_band_contact_before_rounding"]
            < c["finite_band_current_display"],
            "high_transfer_contact_not_zero": actual != 0,
        },
    }
