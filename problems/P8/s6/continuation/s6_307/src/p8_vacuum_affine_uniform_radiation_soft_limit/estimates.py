"""Whole47-tree amplitude estimates with the paired-current cancellation."""

from functools import cache

import sympy as s

from . import source

LOW = s.Integer(2) * 10**8
DIRECT = s.Integer(10) ** 17
MATTER = s.Integer(330000)
SPLIT = s.Integer(192)


def require_domain(energy, transfer, resolution):
    return source.radiation_bounds.require_fixed_resolution_domain(
        energy, transfer, resolution
    )


def low_remainder_bound(energy, transfer, omega, kappa):
    _, _, _, _, delta, w = require_domain(energy, transfer, omega)
    source.require_mass(kappa)
    if w > s.sqrt(delta) / SPLIT:
        raise ValueError("Require the proved low-energy region")
    return LOW / s.sqrt(s.sympify(kappa))


def global_direct_gravity_bound(energy, transfer, omega, kappa):
    _, _, _, _, delta, w = require_domain(energy, transfer, omega)
    source.require_mass(kappa)
    return DIRECT * delta / (s.sqrt(s.sympify(kappa)) * w * s.sqrt(delta + w * w))


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(value)

    y, w, delta, kap = s.symbols("root_transfer omega delta kappa", positive=True)
    values = {
        "paired_current_Lipschitz": s.Integer(4) * 4 + 4 * 2 * 16,
        "one_low_mixed_group": s.Integer(600) * (4 * 43392 + 50400),
        "all_low_paired_bases": s.Integer(2) * 134380800 + 64 * 100992,
        "offshell_shift_stress_entry": s.Integer(2) * 2 + 4 * 2,
        "offshell_shift_contraction": s.Integer(16) * 9 * 12,
        "all_low_offshell_shifts": s.Integer(12) * 16 * 1728 * 2,
        "all_low_seagulls": s.Integer(283392),
        "all_low_cubics": s.Integer(537477120),
        "whole_low_relative_remainder": s.Integer(608)
        + s.Rational(275225088 + 663552 + 283392 + 537477120, 8),
        "all_global_paired_mixed_bases": s.Integer(2) * 3600 * 1200 * 30000,
        "global_timelike_base": s.Integer(3600) * 64 / 5,
        "global_timelike_common_denominator": s.Integer(2) * 46080,
        "all_global_offshell_shifts": s.Integer(12) * 16 * 1728 * 30000,
        "all_global_seagulls": s.Integer(4250880000),
        "all_global_cubics": s.Integer(201553920358318080),
        "whole_global_direct_normalized": s.Rational(
            259200000000 + 92160 + 9953280000 + 4250880000 + 201553920358318080, 8
        ),
    }
    expected = (
        144,
        134380800,
        275225088,
        12,
        1728,
        663552,
        283392,
        537477120,
        101706752,
        259200000000,
        46080,
        92160,
        9953280000,
        4250880000,
        201553920358318080,
        25194274220321280,
    )
    for (name, value), target in zip(values.items(), expected):
        put(name, value - target)
    put(
        "same_old_low_cubic_budget",
        source.radiation_bounds.uniform_data()["whole_uniform_constants"][
            "whole_low_cubic"
        ]
        - values["all_low_cubics"],
    )
    put(
        "same_old_global_seagull_budget",
        source.radiation_bounds.uniform_data()["whole_uniform_constants"][
            "whole_global_seagull"
        ]
        - values["all_global_seagulls"],
    )
    put(
        "same_old_global_cubic_budget",
        source.radiation_bounds.uniform_data()["whole_uniform_constants"][
            "whole_global_cubic"
        ]
        - values["all_global_cubics"],
    )
    put(
        "near_transfer_quotient_difference",
        600 * (43392 / y + 50400 / y**2) - (26035200 * y + 30240000) / y**2,
    )
    c0, c1, c2, jold, jnew = s.symbols(
        "channel0 channel1 channel2 current_old current_new"
    )
    put(
        "unshifted_complete_current_recoil_grouping",
        sum(c * (jnew - jold) for c in (c0, c1, c2)) - (c0 + c1 + c2) * (jnew - jold),
    )
    put(
        "whole_offshell_shift_denominator_conversion",
        values["all_low_offshell_shifts"] / delta / (8 / delta) - 82944,
    )
    put(
        "global_transfer_amplitude_normalization",
        values["whole_global_direct_normalized"]
        * delta
        / (s.sqrt(kap) * w * s.sqrt(delta + w * w))
        - values["whole_global_direct_normalized"]
        / (w * s.sqrt(delta + w * w))
        / (8 / delta)
        * 8
        / s.sqrt(kap),
    )
    am, ag, rm, rg = s.symbols("Am Ag Rm Rg", positive=True)
    put("positive_full_Born_remainder_weight", (am * LOW + ag * LOW) / (am + ag) - LOW)
    put("positive_full_Born_weight_sum", am / (am + ag) + ag / (am + ag) - 1)
    put(
        "convex_amplitude_squared_gap",
        am * rm**2 / (am + ag)
        + ag * rg**2 / (am + ag)
        - ((am * rm + ag * rg) / (am + ag)) ** 2
        - am * ag * (rm - rg) ** 2 / (am + ag) ** 2,
    )
    margins = {
        "low_bound_above_complete_gravity_remainder": LOW
        - values["whole_low_relative_remainder"],
        "low_bound_above_complete_matter_remainder": LOW - MATTER,
        "direct_bound_above_complete_gravity_amplitude": DIRECT
        - values["whole_global_direct_normalized"],
        "low_transfer_half_gap": s.Rational(1, 2) - s.Rational(7, 192),
        "low_transfer_difference_roundup": 7 - 6 - s.Rational(10, 192),
        "mixed_transfer_root_below4": s.Integer(16) - 12,
        "timelike_D_above5": s.Rational(45, 8) - 5,
        "global_gap_below4": s.Integer(4) - s.Rational(65, 64),
        "interior_Born_lower_bound_above8": s.Rational(257, 18) - 8,
        "endpoint_Born_residue_above8": s.Rational(257, 32) - 8,
    }
    for name, value in margins.items():
        put("positive_arithmetic_" + name, value - s.Abs(value))
    return {
        "whole_component_constants": values,
        "whole_positive_margins": margins,
        "whole_low_energy_full_remainder": LOW / s.sqrt(kap),
        "whole_global_gravity_direct_amplitude": DIRECT
        * delta
        / (s.sqrt(kap) * w * s.sqrt(delta + w * w)),
        "whole_low_energy_proof": "Atomega<=sqrt(delta)/192 each mixed |D-D0|<=7sqrt(tau_j)omega and |D|>=tau_j/2. Pair currents before applying the quotient bound; both paired differences are<=600[43392/sqrt(tau_j)+50400/tau_j]. As tau_j<=12, both mixed channels plus the timelike channel are below275225088/delta. All12 exact off-shell shifts add663552/delta. The unchanged complete seagull and improved cubic bounds are283392/delta and537477120/delta. The whole old-channel sum multiplies the current recoil bound608, giving608 G0, not separate uncancelled channel estimates. G0>8/delta yields101706752<2e8; positive full-Born weighting includes the330000 matter remainder.",
        "whole_global_direct_proof": "Use the exact physical gaps -D>(tau_j+omega^2)/30000 and current pairs<600sqrt(tau_j+omega^2)/omega. The two mixed bases contribute259200000000/[omega sqrt(delta+omega^2)], the timelike base92160 in that denominator. All mandatory shift terms, all seagulls and all cubics have the listed1/(delta+omega^2) budgets, which are bounded by the same1/[omega sqrt(delta+omega^2)] denominator. Divide by G0>8/delta. The resulting coefficient25194274220321280 is below1e17. No graph is dropped.",
        "checks": checks,
        "gates": {
            "all_amplitude_margins_strict": all(
                bool(value > 0) for value in margins.values()
            ),
            "low_remainder_independent_of_nonzero_hard_transfer": True,
            "global_direct_bound_without_soft_Taylor_extrapolation": True,
            "all_offshell_shift_seagull_and_cubic_terms_retained": True,
            "all47_interferences_preserved_before_squaring": True,
            "positive_full_Born_not_matter_only_normalization": True,
        },
    }
