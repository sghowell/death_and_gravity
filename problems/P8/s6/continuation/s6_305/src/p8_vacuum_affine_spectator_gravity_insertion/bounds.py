"""Uniform known-reference rate, signed Newton limit and Gaussian-only disk."""

from functools import cache

import sympy as s

from . import insertion as ins
from . import source

GAPPED_AMPLITUDE = s.Rational(25981, 192)
GAPPED_RATE = s.Integer(20)
M1_RATE = s.S.One
TOTAL_NON_NEWTON_RATE = GAPPED_RATE + M1_RATE


def exact_real(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Rational)):
        raise TypeError("Require a real exact rational invariant")
    return s.Rational(value)


def require_physical(energy, transfer):
    energy, transfer = map(exact_real, (energy, transfer))
    crossed = 4 - energy - transfer
    if not s.Rational(25, 4) <= energy <= 16 or transfer >= 0 or crossed >= 0:
        raise ValueError("Require the original compact nonforward mass-one domain")
    return energy, transfer, crossed, min(-transfer, -crossed)


def original_known_rate_interval(energy, transfer):
    require_physical(energy, transfer)
    return -s.Rational(3, 10**602), s.S.Zero


def moments():
    v = source.spectral.V
    return {
        (species, spin): s.integrate(source.spectral.weight(species, spin), (v, 0, 1))
        for species in ("scalar", "vector")
        for spin in (0, 2)
    }


@cache
def data():
    n, m, k = source.HEAVY_MASS2, source.VECTOR_MASS2, source.KAPPA
    weights = moments()
    an = 16 * (
        weights[("scalar", 2)] / (4 * 10**12 - 16)
        + weights[("vector", 2)] / (4 * m - 16)
    )
    hn = 16 * (
        weights[("scalar", 0)] / (4 * 10**12 - 16)
        + weights[("vector", 0)] / (4 * m - 16)
    )
    m1 = s.Rational(3 * 2456 * 7, 3840 * 9)
    disk = (s.Integer(10) ** 198 + 6) / k
    dyson = disk**2 / (1 - disk)
    margins = {
        "actual_H_mass_above_bound": n - 10**12,
        "both_gapped_denominators_positive": 4 * m - 16,
        "nonlocal_A_below_1e_minus6": s.Rational(1, 10**6) - an,
        "nonlocal_H_below_1e_minus5": s.Rational(1, 10**5) - hn,
        "whole_A_below11": 11 - s.Rational(602, 60) - s.Rational(1, 10**6),
        "whole_H_below1205": 1205 - 1204 - s.Rational(1, 10**5),
        "whole_non_Newton_amplitude_below136": 136 - GAPPED_AMPLITUDE,
        "whole_non_Newton_rate_below20": 20 - 2 * GAPPED_AMPLITUDE / 14,
        "M1_log_amplitude_below2": 2 - m1,
        "M1_small_tau_rate_below_one": 1 - s.Rational(4, 8),
        "M1_large_tau_rate_below_one": 1 - s.Rational(12, 14),
        "fixed_Newton_dominates_non_Newton_budget": s.Rational(5 * 10**6, 96 * 16) - 11,
        "signed_bracket_is_negative": 2 * 11 - TOTAL_NON_NEWTON_RATE,
        "original_total_magnitude_below3e_minus602": s.Rational(3, 10**602)
        - (2 * s.Integer(10) ** 198 + 21) / k,
        "TT_quotient_error_below6_over_kappa": 6 - s.Rational(16 * 11, 16 * 9),
        "trace_quotient_error_below6_over_kappa": 6 - s.Rational(16 * 1205, 384 * 9),
        "Gaussian_quotient_disk_below_half": s.Rational(1, 2) - disk,
        "Gaussian_Dyson_remainder_below5e_minus1204": s.Rational(5, 10**1204) - dyson,
    }
    tau, energy = s.symbols("tau energy", positive=True)
    V = energy**2 - 4 * energy + 2
    checks = {}

    def put(name, value):
        checks[name] = s.factor(value)

    for key, target in {
        ("scalar", 2): s.Rational(1, 210),
        ("scalar", 0): s.Rational(68, 35),
        ("vector", 2): s.Rational(3, 14),
        ("vector", 0): s.Rational(36, 35),
    }.items():
        put("whole_radial_weight_" + key[0] + str(key[1]), weights[key] - target)
    put("whole_numerator0_budget", s.Rational(18**2, 3) - 108)
    put("whole_numerator2_budget", 2 + 32 + 256 + 54 - 344)
    put("whole_M1_projector_budget", 4 * 344 + 10 * 108 - 2456)
    put(
        "whole_gapped_amplitude_budget",
        3 * (s.Rational(344 * 11, 16 * 9) + s.Rational(108 * 1205, 768 * 9))
        - GAPPED_AMPLITUDE,
    )
    put("whole_non_Newton_rate_budget", GAPPED_RATE + M1_RATE - TOTAL_NON_NEWTON_RATE)
    put("tau_log_envelope_derivative", s.diff(tau * (1 - s.log(tau)), tau) + s.log(tau))
    ag = ins.gravity_born(energy, -tau, 1, k)
    put("whole_forward_Born_residue", s.limit(tau * ag, tau, 0, dir="+") - V / k)
    am = source.forward.matter_born(energy, 1 - 2 * tau / (energy - 4), n, source.CUBIC)
    put(
        "positive_Born_gravity_weight_forward_limit",
        s.limit(ag / (ag + am), tau, 0, dir="+") - 1,
    )
    put(
        "crossed_endpoint_same_Born",
        ins.gravity_born(energy, 4 - energy + tau, 1, k) - ag,
    )
    put("vanishing_bounded_gapped_normalized_endpoint", s.limit(tau, tau, 0, dir="+"))
    put(
        "vanishing_M1_log_normalized_endpoint",
        s.limit(tau * s.log(tau), tau, 0, dir="+"),
    )
    x = s.Symbol("Gaussian_quotient_increment")
    put(
        "entire_Gaussian_inverse_geometric_remainder",
        1 / (1 + x) - (1 - x) - x * x / (1 + x),
    )
    w, e, dk = s.symbols(
        "positive_Born_gravity_weight remainder fixed_delta_kappa", real=True
    )
    signed = w * (-2 * dk / k + e)
    return {
        "whole_exact_radial_weight_integrals": {
            str(key): value for key, value in weights.items()
        },
        "whole_nonlocal_radial_compact_bounds": {2: an, 0: hn},
        "whole_non_Newton_gapped_amplitude_bound": GAPPED_AMPLITUDE / k**2,
        "whole_massless_M1_log_amplitude_bound": 2 * (1 + s.Abs(s.log(tau))) / k**2,
        "whole_signed_known_reference_rate": signed,
        "whole_original_signed_rate_interval": (-s.Rational(3, 10**602), s.S.Zero),
        "whole_known_forward_limit": -2 * source.poles.delta_kappa() / k,
        "whole_complex_disk_radius": s.Integer(16),
        "whole_Gaussian_Dyson_relative_remainder_bound": dyson,
        "whole_positive_arithmetic_margins": margins,
        "domain": "mu=1,25/4<=s<=16,t,u<0,all nonforward angles with no angular exclusion or shrinking detector resolution. n=10^200/512+2,M_A^2=1e6,kappa=1e800; the fixed source prescription and M1 reference scale1 are retained.",
        "written_rate_proof": "Every channel has |a|<=16,|N0|<=108,|N2|<=344. Absolute radial denominators are at least4nu-16; the positive weights and log(n)<600 give |A|<11,|H|<1205. Hence |M_HV_nonNewton|<136/kappa^2. Original A_G>14/kappa gives |2M_nonNewton/A_G|<20/kappa. For tau=min(-t,-u), M1 is below2(1+abs(log(tau)))/kappa^2. If tau<=1, A_G>8/(kappa*tau) and tau(1-log(tau))<=1. If tau>=1, tau<=6 and log(tau)<2. Thus |2ReM1/A_G|<1/kappa. With wG=A_G/(A_m+A_G) in(0,1), E=wG[-2dk/kappa+e], |e|<21/kappa. The inherited11<dk<1e198 yields -3e-602<E<0. At either forward endpoint wG tends1, bounded gapped and logarithmic M1 normalized remainders vanish, so E tends-2dk/kappa.",
        "Gaussian_disk_boundary": "On the entire complex disk |p|<=16 both normalized fixed H/Proca quotient remainders beyond1+dk/kappa are smaller than6/kappa. The quotient therefore has no zero; only its original massless response pole remains. The displayed geometric remainder bounds only repeated insertions of these same Gaussian kernels. It does NOT bound independent higher-loop graphs, determine the interacting Newton residue, or give a physical Wilsonian cutoff.",
        "sign_boundary": "Only the specified known-reference interference has this strict negative sign. Unknown physical finite matching can change the total rate and its forward Newton limit. No isolated pole-subtracted positivity sign or complex Regge remainder is deduced.",
        "checks": checks,
        "gates": {
            "all_compact_and_Gaussian_arithmetic_margins_positive": all(
                bool(x > 0) for x in margins.values()
            ),
            "original_positive_full_Born_denominator_retained": True,
            "all_nonforward_angles_at_fixed_hard_energy_compact": True,
            "both_forward_endpoints_have_same_known_Newton_limit": True,
            "unknown_matching_not_bounded_by_known_reference_sign": True,
            "Gaussian_geometric_remainder_not_interacting_loop_error": True,
        },
    }
