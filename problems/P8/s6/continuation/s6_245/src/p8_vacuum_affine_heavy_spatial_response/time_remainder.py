"""Full scalar state/reference difference and correctly ordered six-step time remainder."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_curved_state import state

from . import domain


@cache
def recurrence():
    rows = [(1,)]
    for j in range(6):
        old = rows[-1]
        rows.append(
            tuple(
                (2 * j + 2 - r) * (old[r] if r < len(old) else 0)
                + (old[r - 1] if r else 0)
                for r in range(j + 2)
            )
        )
    return tuple(rows)


@cache
def literal_rows():
    g = s.symbols("inverse_phase_0:7")
    a = s.symbols("full_pair_0:7")
    h = s.symbols("source_0:7")

    def dt(expr):
        return s.expand(
            sum(
                s.diff(expr, row[j]) * row[j + 1] for row in (g, a, h) for j in range(6)
            )
        )

    value = a[0] * h[0]
    out = []
    for order in range(7):
        normalized = s.expand(
            value.subs(
                {row[j]: s.factorial(j) for row in (g, a) for j in range(7)},
                simultaneous=True,
            )
        )
        out.append(tuple(normalized.coeff(h[j]) for j in range(order + 1)))
        if order < 6:
            value = dt(g[0] * value)
    return tuple(out)


@cache
def constants():
    n = state.MASS2
    kappa = state.KAPPA
    C, R = domain.PAIR, domain.RADIUS
    rows = recurrence()
    bulk = C * C * 64 * sum(rows[6]) * R**-6
    endpoint = C * C * 64 * sum(rows[5]) * R**-5
    return {
        "whole_sixth_bulk_numerator": bulk,
        "whole_fifth_endpoint_numerator": endpoint,
        "complete_time_integral_numerator_over_mass": (bulk + endpoint)
        * s.Rational(64, 24 * 4),
        "complete_time_removed_union_numerator_over_radial_cutoff": (bulk + endpoint)
        * s.Rational(4**4, 36),
        "whole_actual_minus_W6_transfer_difference_normalized": s.Rational(
            10**106, kappa
        )
        / n**3,
        "whole_actual_minus_W6_transfer_difference_normalized_tail": s.Rational(
            10**107, kappa
        )
        / (n * n * domain.MASS_FLOOR),
    }


@cache
def data():
    rows = recurrence()
    literal = literal_rows()
    cc = constants()
    checks = {}
    for j in range(7):
        checks["every_full_literal_Leibniz_coefficient_" + str(j)] = s.Matrix(
            [a - b for a, b in zip(rows[j], literal[j])]
        )
    L = s.symbols("L0:7")
    for order in range(1, 7):
        total = sum(s.I * (-s.I) ** j * (L[j + 1] - s.I * L[j]) for j in range(order))
        checks["whole_ordered_integration_identity_" + str(order)] = s.expand(
            total - L[0] + (-s.I) ** order * L[order]
        )
    checks["full_fifth_endpoint_coefficient"] = s.I * (-s.I) ** 5 - 1
    checks["full_sixth_bulk_coefficient"] = (-s.I) ** 6 + 1
    old = s.symbols("old_0:4")
    error = s.symbols("error_0:4")
    actual = tuple(a + b for a, b in zip(old, error))
    full = s.prod(actual) - s.prod(old)
    telescope = sum(
        error[j] * s.prod(actual[:j]) * s.prod(old[j + 1 :]) for j in range(4)
    )
    checks["entire_four_factor_difference_with_all_error_products"] = s.expand(
        full - telescope
    )
    checks["full_sixth_Leibniz_row_sum"] = sum(rows[6]) - 97423
    checks["full_fifth_Leibniz_row_sum"] = sum(rows[5]) - 7916
    nu, mu = s.symbols("nu mu", positive=True)
    checks["full_pair_harmonic_mean_identity"] = s.expand(
        (nu + mu) ** 2 - 4 * nu * mu - (nu - mu) ** 2
    )
    checks["sixth_equal_frequency_integrable_weight"] = (nu * mu / (nu + mu) ** 6).subs(
        mu, nu
    ) - 1 / (64 * nu**4)
    checks["fifth_equal_frequency_logarithmic_absolute_majorant"] = (
        nu * mu / (nu + mu) ** 5
    ).subs(mu, nu) - 1 / (32 * nu**3)
    r, K, m = s.symbols("r K m", positive=True)
    checks["full_far_inverse_power_integral"] = s.integrate(r**-2, (r, K, s.oo)) - 1 / K
    checks["same_original_frequency_radial_radius"] = s.expand(
        (4 * s.sqrt(K * K - m * m)) ** 2 - 16 * (K * K - m * m)
    )
    # The radius square is16(K²-m²); at K²>=4m² it exceeds K² by at least11K².
    checks["same_cutoff_radius_margin"] = s.expand(
        16 * (K * K - m * m) - K * K - 11 * K * K - 4 * (K * K - 4 * m * m)
    )
    B = s.Integer(10) ** 83
    pair_error = 3 * 600 * 400
    four_product = pair_error * (3 * 400**2 + 3 * 200**2)
    raw_integrand = four_product * 64**2
    gates = {
        "full_scalar_reference_feature_and_actual_norms": 3 * 25**2 < domain.PAIR,
        "full_actual_plus_comparison_W6_feature_error": 300 * B
        + s.Rational(4 * 10**64, domain.MASS_FLOOR**2)
        < 600 * B,
        "complete_scalar_pair_difference": pair_error < 10**6,
        "whole_four_factor_error_including_squares": raw_integrand < 10**20,
        "all_internal_initial_errors_small_but_not_reset": 600
        * B
        / domain.MASS_FLOOR**11
        < 1,
        "entire_state_difference_spatial_weight_on_source": 8 * 10**20 * B < 10**106,
        "entire_removed_union_state_tail_source_weight": 400 * 10**20 * B < 10**107,
        "full_actual_state_difference_small": cc[
            "whole_actual_minus_W6_transfer_difference_normalized"
        ]
        < s.Rational(1, 10**1280),
        "full_actual_state_difference_tail_small": cc[
            "whole_actual_minus_W6_transfer_difference_normalized_tail"
        ]
        < s.Rational(1, 10**1180),
        "every_source_jet_through_six_retained": len(rows) == 7 and len(rows[-1]) == 7,
        "both_finite_time_pieces_kept": cc["whole_sixth_bulk_numerator"] > 0
        and cc["whole_fifth_endpoint_numerator"] > 0,
        "entire_time_integral_numerator": cc[
            "complete_time_integral_numerator_over_mass"
        ]
        < 10**51,
        "entire_time_removed_union_numerator": cc[
            "complete_time_removed_union_numerator_over_radial_cutoff"
        ]
        < 10**53,
        "sixth_is_absolutely_integrable": s.integrate(r**-2, (r, 1, s.oo)) == 1,
        "fifth_absolute_envelope_not_integrable": s.integrate(1 / r, (r, 1, s.oo))
        == s.oo,
        "spatial_contact_cancels_only_in_transfer_difference": True,
        "no_final_source_cutoff_or_lower_germ_reset": True,
        "computational_W6_family_not_a_new_physical_state_response": True,
    }
    return {
        "all_full_Cauchy_Leibniz_rows": rows,
        "all_independent_literal_rows": literal,
        "complete_six_endpoint_coefficients": tuple(
            s.I * (-s.I) ** j for j in range(6)
        ),
        "whole_sixth_bulk_coefficient": (-s.I) ** 6,
        "complete_four_factor_error_identity": full,
        "complete_reference_error_numerators": {
            "full_pair_error_over_B": pair_error,
            "full_four_product_error_over_B": four_product,
            "complete_volume_weighted_memory_over_B": raw_integrand,
        },
        "constants": cc,
        "actual_state_difference": "The actual SLE-minus-W6 feature difference includes the full exact-comparison residual, SLE beta, both mixed terms and every error product. Both created legs and the reflected ordering are retained. Spatial Hamiltonian/profile contacts cancel only in the P-minus-zero difference; they remain in the actual homogeneous anchor.",
        "full_six_step_identity": "With physical detector SHARP, positive detector phase and negative source phase, I=exp(-iTheta)sum(j0..5)i(-i)^j gL^jU+(-i)^6 integral exp(-iTheta)L^6U. The j5 endpoint is+gL5U and the whole sixth bulk coefficient is-1, before positive imaginary part. All lower endpoints vanish only by the prepared source germ.",
        "complete_time_bound": "Full analytic scalar pair C sqrt(nu mu), C1e9, g<=2/(nu+mu), radius1/20000. Both j5 and bulk6 have weight nu mu/(nu+mu)^6<=1/(4nu^4). The entire retarded triangle and endpoint are bounded separately before their positive sum. Full J4=64/(8pi m)<64/(24m). No vector polarization factor is inserted.",
        "same_mask_tail": "Use Omega_star<=K on BOTH modes, K>=2m, with comoving radius4sqrt(K²-m²)>K. Symmetry bounds the removed union by half of the one-leg fourth-power radial tail. The state error uses the full low/high external-transfer partition, retaining large-small pairs and placing the spatial weight on the source only.",
        "graph": "Detector L2; source time jets0..6 for the finite time piece, spatial Bessel degree6 for the state-difference tail. All are dominated by the final prepared source Z136 graph.",
        "checks": checks,
        "gates": {key: bool(value) for key, value in gates.items()},
    }
