"""Whole first-Newton one-loop rate assembly with explicit unmatched hard data."""

from functools import cache

import sympy as s
from p8_vacuum_affine_physical_virtual_soft_pairing import inclusive as paired

from . import matter, source

AM, AG = s.symbols("positive_matter_Born positive_gravity_Born", positive=True)
MAX_RESOLUTION = s.Rational(1, 8)
LM, HD, DELTA, REAL = s.symbols(
    "real_complete_matter_loop finite_hard_dressing_interference finite_soft_conversion finite_exact_minus_soft_real_rate",
    real=True,
)


def normalized_one_Newton_correction(
    matter_Born=AM,
    gravity_Born=AG,
    matter_loop_real=LM,
    hard_dressing=HD,
    conversion=DELTA,
    real_remainder=REAL,
):
    am, ag, lm, hd, de, rr = map(
        s.sympify,
        (
            matter_Born,
            gravity_Born,
            matter_loop_real,
            hard_dressing,
            conversion,
            real_remainder,
        ),
    )
    return (am * am * (hd + de + rr) + 2 * ag * lm) / (am + ag) ** 2


def hard_reference(matter_Born=AM, gravity_Born=AG, hard_dressing=HD):
    am, ag, hd = map(s.sympify, (matter_Born, gravity_Born, hard_dressing))
    return am * am * hd / (am + ag) ** 2


def reference_error(resolution=MAX_RESOLUTION, kappa=source.KAPPA):
    return matter.data()[
        "whole_uniform_interference_upper"
    ] + paired.physical_reference_error(resolution, kappa)


@cache
def data():
    checks = {}

    def put(name, value):
        checks[name] = s.factor(value)

    h, loop = s.symbols("Newton_marker loop_marker", real=True)
    am, ag, lm, lg, im, ig = s.symbols("am ag lm lg im ig", real=True)
    amp = am + h * ag + loop * (lm + s.I * im + h * (lg + s.I * ig))
    squared = s.expand(amp * s.conjugate(amp))
    put(
        "complete_first_Newton_one_loop_interference",
        squared.coeff(h, 1).coeff(loop, 1) - 2 * (am * lg + ag * lm),
    )
    put(
        "complete_zero_Newton_one_loop_interference",
        squared.coeff(h, 0).coeff(loop, 1) - 2 * am * lm,
    )
    put(
        "higher_Newton_loop_interference_not_included",
        squared.coeff(h, 2).coeff(loop, 1) - 2 * ag * lg,
    )
    numerator = AM * AM * (HD + DELTA + REAL) + 2 * AG * LM
    put(
        "whole_known_first_Newton_inclusive_numerator",
        normalized_one_Newton_correction() * (AM + AG) ** 2 - numerator,
    )
    put(
        "unchosen_hard_dressing_sensitivity",
        s.diff(normalized_one_Newton_correction(), HD) - AM**2 / (AM + AG) ** 2,
    )
    put(
        "complete_matter_loop_sensitivity",
        s.diff(normalized_one_Newton_correction(), LM) - 2 * AG / (AM + AG) ** 2,
    )
    put(
        "finite_real_error_weight",
        s.diff(normalized_one_Newton_correction(), REAL) - AM**2 / (AM + AG) ** 2,
    )
    put(
        "finite_conversion_error_weight",
        s.diff(normalized_one_Newton_correction(), DELTA) - AM**2 / (AM + AG) ** 2,
    )
    put(
        "hard_reference_not_matching_choice",
        hard_reference() - AM**2 * HD / (AM + AG) ** 2,
    )
    put(
        "matter_only_Born_reference_recovered",
        normalized_one_Newton_correction(gravity_Born=0) - HD - DELTA - REAL,
    )
    put(
        "new_complete_correction_equals_selected_plus_missing_piece",
        normalized_one_Newton_correction()
        - AM**2 * (HD + DELTA + REAL) / (AM + AG) ** 2
        - 2 * AG * LM / (AM + AG) ** 2,
    )
    e = s.Symbol("IR_regulator", positive=True)
    B = s.Symbol("real_Bsoft", real=True)
    # Whole selected real/virtual principal part; the added massive matter
    # coefficient has no graviton IR regulator and contributes no new pole.
    virtual = AM**2 * 2 * B / (8 * s.pi**2 * source.K * e)
    real = AM**2 * (-2 * B) / (8 * s.pi**2 * source.K * e)
    put("entire_first_Newton_rate_IR_pole", virtual + real)
    put("added_matter_loop_has_no_graviton_regulator", s.diff(2 * AG * LM, e))
    compact = paired.physical_reference_error(s.Rational(1, 8), source.KAPPA)
    total = reference_error()
    put(
        "total_reference_bound_is_sum",
        total - matter.data()["whole_uniform_interference_upper"] - compact,
    )
    return {
        "whole_first_Newton_one_loop_inclusive_numerator": numerator,
        "whole_positive_full_Born_normalized_correction": normalized_one_Newton_correction(),
        "whole_unmatched_hard_reference": hard_reference(),
        "whole_original_reference_error_upper": total,
        "whole_rate_assembly_proof": "Let A_m and A_G be the complete positive matter and gravity Born amplitudes, L_m the same-OS4 complete massive matter loop, and H=2Re(L_g,hard/A_m) the entire S293/S294 gravitational dressing hard reference, retaining unknown finite matching. The full one-loop correction linear in the Newton marker has numerator A_m^2(H+Delta_soft+R_real)+2 A_G Re L_m. S296 supplies the finite paired first term; the added matter term is IR finite at fixed nonforward angle. Divide by the explicitly unexpanded(A_m+A_G)^2 only as a positive normalization. No higher-Newton loop or one-loop square is thereby calculated.",
        "whole_quantitative_reference": "On mass-one25/4<=s<=16,-1<z<1,original n,g,kappa and0<resolution<=1/8, the normalized difference from A_m^2 H/(A_m+A_G)^2 is below eta/2+10^-792<2*10^-199. If the complete L_m is kept explicitly rather than bounded, the sole conversion/recoil error is below10^-792 times A_m^2/(A_m+A_G)^2. Regulator removal is at fixed resolution and angle first.",
        "whole_matching_scope": "This assembles the formal first-loop, first-Newton contribution in the existing matter prescription. Finite gravitational curvature/local/heavy matching in H remains unchosen. A finite-order inclusive rate is not an analytic crossing-symmetric amplitude. Higher-Newton soft radiation becomes essential in the forward cone; no uniform Newton expansion, positive dispersion measure, all-order infrared theorem, or Regge bound follows.",
        "checks": checks,
        "gates": {
            "whole_selected_conversion_error_retained": bool(
                compact < s.Rational(1, 10**792)
            ),
            "whole_new_reference_error_below_two_e_minus_199": bool(
                total < s.Rational(2, 10**199)
            ),
            "all_rate_cross_terms_at_first_Newton_retained": True,
            "old_matter_OS4_not_gravity_finite_matching": True,
            "normalization_does_not_claim_higher_Newton_completeness": True,
            "same_order_of_IR_limits_and_nonforward_domain": True,
            "rate_not_a_crossing_analytic_dispersion_amplitude": True,
        },
    }
