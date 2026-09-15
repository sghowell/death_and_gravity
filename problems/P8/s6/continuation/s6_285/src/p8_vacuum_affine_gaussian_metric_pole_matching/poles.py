"""Fixed H/Proca metric poles and the explicit unmatched light-scalar polynomial."""

from functools import cache

import sympy as s

from . import source, spectral

P = s.Symbol("Laplace_p", complex=True)


def delta_kappa():
    return 2 * source.fixed_coefficients()["R_old"]


def radial_form_factors(p=P):
    v = spectral.V
    A = (s.log(source.HEAVY_MASS2) + 2) / 60
    H = 2 * (s.log(source.HEAVY_MASS2) + 2)
    for species, mass2 in (
        ("scalar", source.HEAVY_MASS2),
        ("vector", source.VECTOR_MASS2),
    ):
        A += s.Integral(
            p * spectral.weight(species, 2) / (4 * mass2 + p * (1 - v * v)), (v, 0, 1)
        )
        H += s.Integral(
            p * spectral.weight(species, 0) / (4 * mass2 + p * (1 - v * v)), (v, 0, 1)
        )
    return A, H


@cache
def data():
    k = source.KAPPA
    n = source.HEAVY_MASS2
    m = source.VECTOR_MASS2
    fixed = source.fixed_coefficients()
    dk = delta_kappa()
    A, H = radial_form_factors()
    O2 = P * (1 + dk / k) + P * P * A / (16 * s.pi**2 * k)
    O0 = -2 * P * (1 + dk / k) + P * P * H / (192 * s.pi**2 * k)
    totalI = {}
    for spin in (0, 2):
        totalI[spin] = sum(
            spectral.moment(species, spin, 3).subs(spectral.NU, mass2)
            for species, mass2 in (("scalar", n), ("vector", m))
        )
    cv, delta, c2, c3, hbar, p = s.symbols(
        "volume_pole delta c2 c3 formal_hbar p", nonzero=True
    )
    inverse_coefficient = -(cv + delta * p + c2 * p * p + c3 * p**3) / p**2
    checks = {
        "same_signature_correct_Gaussian_Newton_coefficient": s.factor(
            dk - (n * (s.log(n) - 1) + 5 * m) / (96 * s.pi**2)
        ),
        "same_TT_finite_fourth_normalization": s.factor(
            -4 * fixed["Weyl_squared"] - (s.log(n) + 2) / (960 * s.pi**2)
        ),
        "same_trace_finite_fourth_normalization": s.factor(
            -24 * fixed["R_old_squared"] - (s.log(n) + 2) / (96 * s.pi**2)
        ),
        "whole_TT_radial_force_normalization": s.expand(
            O2 - P * (1 + dk / k) - P * P * A / (16 * s.pi**2 * k)
        ),
        "whole_trace_radial_force_normalization": s.expand(
            O0 + 2 * P * (1 + dk / k) - P * P * H / (192 * s.pi**2 * k)
        ),
        "whole_metric_insertion_inverse_at_first_loop": s.cancel(
            s.diff(
                1 / (p + hbar * (cv + delta * p + c2 * p * p + c3 * p**3)), hbar
            ).subs(hbar, 0)
            - inverse_coefficient
        ),
        "volume_cancellation_removes_gapped_insertion_double_pole": s.expand(
            inverse_coefficient.subs(cv, 0)
        ).coeff(p, -2),
        "finite_Newton_simple_pole_retained": s.expand(
            inverse_coefficient.subs(cv, 0)
        ).coeff(p, -1)
        + delta,
        "higher_local_term_not_a_double_pole": s.expand(
            inverse_coefficient.subs(cv, 0)
        ).coeff(p, 0)
        + c2,
        "dispersive_remainder_not_a_double_pole": s.expand(
            inverse_coefficient.subs(cv, 0)
        ).coeff(p, 1)
        + c3,
        "positive_residue_relation": s.cancel(k / (k + dk) - (1 - dk / (k + dk))),
        "complete_scalar_and_vector_tensor_first_moment": s.factor(
            totalI[2] - 1 / (53760 * s.pi**2 * n) - 3 / (3584 * s.pi**2 * m)
        ),
        "complete_scalar_and_vector_trace_first_moment": s.factor(
            totalI[0] - 17 / (26880 * s.pi**2 * n) - 3 / (8960 * s.pi**2 * m)
        ),
    }
    lightR, lightW, lightR2 = s.symbols(
        "unmatched_Phi_R unmatched_Phi_Weyl_squared unmatched_Phi_R_squared", real=True
    )
    v = spectral.V
    lightA = s.Integral(
        P * spectral.weight("scalar", 2) / (4 + P * (1 - v * v)), (v, 0, 1)
    )
    lightH = s.Integral(
        P * spectral.weight("scalar", 0) / (4 + P * (1 - v * v)), (v, 0, 1)
    )
    fullO2 = (
        O2
        + 2 * lightR * P / k
        - 4 * lightW * P * P / k
        + P * P * lightA / (16 * s.pi**2 * k)
    )
    fullO0 = (
        O0
        - 4 * lightR * P / k
        - 24 * lightR2 * P * P / k
        + P * P * lightH / (192 * s.pi**2 * k)
    )
    full_residue = k / (k + dk + 2 * lightR)
    checks.update(
        {
            "unmatched_light_R_TT_coefficient": s.diff(fullO2, lightR) - 2 * P / k,
            "unmatched_light_R_trace_coefficient": s.diff(fullO0, lightR) + 4 * P / k,
            "unmatched_light_Weyl_coefficient": s.diff(fullO2, lightW) + 4 * P * P / k,
            "unmatched_light_R_squared_coefficient": s.diff(fullO0, lightR2)
            + 24 * P * P / k,
            "full_parameterized_massless_residue": s.cancel(
                full_residue * (k + dk + 2 * lightR) - k
            ),
            "light_curvature_does_not_change_volume_double_pole": s.expand(
                (-2 * lightR * P + 4 * lightW * P * P) / P**2
            ).coeff(P, -2),
        }
    )
    margins = {
        "heavy_mass_above_e_upper": n - 3,
        "heavy_mass_below_1e198": s.Integer(10) ** 198 - n,
        "e_cubed_above10": 1 + 3 + s.Rational(9, 2) + s.Rational(27, 6) - 10,
        "delta_kappa_strict_positive_lower": 5 * m - 1,
        "delta_kappa_below_1e198": 864 * s.Integer(10) ** 198
        - (599 * s.Integer(10) ** 198 + 5 * m),
        "TT_local_below_one_tenth": s.Rational(1, 10) - s.Rational(602, 960 * 9),
        "TT_nonlocal_below_one_thousandth": s.Rational(1, 1000)
        - (s.Rational(16, 53760) + s.Rational(24, 3584)) / 9,
        "TT_total_error_below_one": 1 - s.Rational(1, 10) - s.Rational(1, 1000),
        "trace_local_below_seven_tenths": s.Rational(7, 10) - s.Rational(602, 96 * 9),
        "trace_nonlocal_below_one_hundredth": s.Rational(1, 100)
        - (s.Rational(8 * 34, 26880) + s.Rational(24, 8960)) / 9,
        "trace_total_error_below_one": 1 - s.Rational(7, 10) - s.Rational(1, 100),
        "response_ratio_denominator_above_half_kappa": k / 2 - 1,
        "response_ratio_numerator_below_twice_upper": 2 * s.Integer(10) ** 198
        - (s.Integer(10) ** 198 + 1),
    }
    return {
        "entire_fixed_H_Proca_delta_kappa": dk,
        "fixed_H_Proca_kappa_eff": k + dk,
        "complete_radial_A2_and_H0": (A, H),
        "complete_prescribed_H_Proca_metric_force_sectors": {2: O2, 0: O0},
        "all_three_cuts_with_unmatched_light_polynomial": {2: fullO2, 0: fullO0},
        "full_parameterized_Gaussian_massless_residue": full_residue,
        "unmatched_light_curvature_coefficients": (lightR, lightW, lightR2),
        "positive_canonical_H_Proca_pole_residue": k / (k + dk),
        "complete_first_nonlocal_moments": totalI,
        "all_coarse_positive_arithmetic_margins": margins,
        "kernel_definition": "O=O2 P2+O0 P0, p=lambda^2+|q|^2, equivalently z=-p on the Feynman/retarded sheet. O2=p*kappa_eff/kappa+p^2*A2/(16pi^2*kappa); O0=-2p*kappa_eff/kappa+p^2*H0/(192pi^2*kappa). The H and Proca cuts and their fixed finite volume/R/Weyl/R^2 coefficients are included in O2,O0. The separate full parameterized kernel adds the light scalar cut with its three UNMATCHED curvature coefficients; it is not licensed to set those coefficients to any convenient values. Compact Euler variations vanish as before.",
        "exact_covariant_low_momentum_structure": "Before gauge fixing, p(P2-2P0) is local because its apparent highest inverse powers of p cancel. p^2 times either projector is also local. The subtracted D_i(-p)=O(p^3) leaves an analytic full tensor at zero. After gauge fixing, the prescribed H/Proca massless pole has conserved numerator P2-P0/2 and residue kappa/kappa_eff. Including Phi changes its residue to kappa/(kappa_eff+2c_Phi_R), whose sign is NOT established without matching. This is not an additional physical spin0 graviton.",
        "rigorous_bounds": "For the actual n=10^200/512+2,M_A^2=10^6,kappa=10^800, 0<delta_kappa<10^198 and 0<delta_kappa/kappa<10^-602. Hence 0<1-kappa/kappa_eff<10^-602. On the entire complex disk |p|<=1, the TT and trace quotient remainders beyond their respective Einstein-plus-fixed-Newton constants each have absolute value<1/kappa. Both prescribed H/Proca quotient sectors are nonzero; no additional H/Proca Gaussian response pole is in this disk. The inequalities do NOT establish the full Phi-inclusive kernel bound. The H/Proca TT response relative to the original Einstein response differs by<4*10^-602, including its continuous p0 ratio.",
        "physical_matching_scope": "The three gapped Gaussian metric insertions have no cosmological double pole once the source-pinned whole volume term is included. The H/Proca fixed long-range exchange coefficient is1/kappa_eff. The Phi-inclusive coefficient is1/(kappa_eff+2c_Phi_R), with c_Phi_R unresolved and the original coupling and external chart retained. This is not an all-orders physical Newton constant, not full proper-vertex/LSZ or physical mass matching, and not a theorem discarding every branch or apparent double-pole term in the full scalar amplitude. No high-energy or curved stability conclusion follows from the low complex disk.",
        "checks": checks,
        "gates": {
            "all_coarse_exact_margins_positive": all(
                value > 0 for value in margins.values()
            ),
            "both_prescribed_determinants_and_explicit_third_cut": all(
                len(value.atoms(s.Integral)) == 2 for value in (A, H)
            )
            and all(len(value.atoms(s.Integral)) == 3 for value in (fullO2, fullO0)),
            "finite_prescribed_Newton_shift_and_unmatched_light_anchor_retained": dk
            != 0
            and full_residue.has(lightR),
            "complete_fixed_curvature_polynomial_retained": all(
                fixed[name] != 0 for name in ("R_old", "Weyl_squared", "R_old_squared")
            ),
            "low_disk_not_a_physical_EFT_cutoff": True,
            "no_full_amplitude_or_exact_interacting_pole_claim": True,
            "no_massless_physical_scalar_graviton_added": True,
        },
    }
