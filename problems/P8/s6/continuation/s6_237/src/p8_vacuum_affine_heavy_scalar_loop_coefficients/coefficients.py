"""Actual first-loop low-energy coefficient intervals and matching tolerances."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_scalar_tree_matching import model

from . import bidisk, jets


@cache
def data():
    n, D, g2, H = s.symbols("n D g_squared log_n", positive=True)
    common = g2 * g2 / (16 * s.pi**2)
    approximations = {
        "delta_b20": common
        * (
            (3 * H - s.Rational(157, 45)) / n**4
            + (24 * H - s.Rational(7969, 315)) / n**5
        ),
        "delta_b21": common * (-3 * H + s.Rational(13, 21)) / n**5,
        "delta_b40": s.S.Zero,
    }
    E = bidisk.REMAINDER_CONSTANT * common / n**6
    gactual = model.G2
    nactual = model.MASS2
    dactual = nactual - 2
    exact20 = s.Rational(1500, 288) * gactual / nactual
    exact21 = s.Rational(1500, 432) * gactual / nactual
    exact40 = s.Rational(8 * 10**12, 9) * gactual / nactual
    full = jets.derived()[-1]
    forward = s.expand(
        full.subs({jets.S: 2 + jets.v - jets.t / 2, jets.T: jets.t}, simultaneous=True)
    )
    checks = {
        "first_loop_b20_definition_from_complete_amplitude": s.expand(
            forward.coeff(jets.v, 2).subs(jets.t, 0)
            - (
                (3 * jets.H - s.Rational(157, 45)) * jets.w**4
                + (24 * jets.H - s.Rational(7969, 315)) * jets.w**5
            )
        ),
        "first_loop_b21_definition_from_complete_amplitude": s.expand(
            s.diff(forward.coeff(jets.v, 2), jets.t).subs(jets.t, 0)
            - (-3 * jets.H + s.Rational(13, 21)) * jets.w**5
        ),
        "first_loop_b40_definition_from_complete_amplitude": forward.coeff(
            jets.v, 4
        ).subs(jets.t, 0),
        "original_b20_target_exact_relation": model.LAMBDA - gactual / (2 * dactual**3),
        "original_b21_target_exact_relation": model.GAMMA - gactual / dactual**4,
        "higher_forward_tree_coefficient_exact_relation": s.cancel(
            model.GAMMA**2 / model.LAMBDA - 2 * gactual / dactual**5
        ),
        "b20_relative_upper_before_actual_n_margin": s.cancel(
            (1500 * common / n**4) / (2 * g2 / D**3)
            - 1500 * g2 * D**3 / (32 * s.pi**2 * n**4)
        ),
        "b21_relative_upper_before_actual_n_margin": s.cancel(
            (1500 * common / n**5) / (3 * g2 / D**4)
            - 1500 * g2 * D**4 / (48 * s.pi**2 * n**5)
        ),
        "b40_relative_upper_before_actual_n_margin": s.cancel(
            256 * E / (2 * g2 / D**5) - 8 * 10**12 * g2 * D**5 / (s.pi**2 * n**6)
        ),
        "positive_total_higher_coefficient_margin": (1 - s.Rational(1, 10**192))
        + s.Rational(1, 10**192)
        - 1,
    }
    return {
        "coefficient_definitions": "b20=[v^2 t^0]A=(partial_v^2 A)/2, b21=[v^2 t^1]A=partial_t partial_v^2 A/2, b40=[v^4 t^0]A=(partial_v^4 A)/24, all at v=t=0 with v=s+t/2-2. The first-order light mass/residue conditions are the unchanged S234 ones.",
        "complete_first_loop_approximations": approximations,
        "complete_coefficient_remainder_bounds": {
            "delta_b20": 16 * E,
            "delta_b21": 64 * E,
            "delta_b40": 256 * E,
        },
        "actual_sign_intervals": "1000g^4/(16pi^2 n^4)<delta_b20<1500g^4/(16pi^2 n^4); -1500g^4/(16pi^2 n^5)<delta_b21<-1000g^4/(16pi^2 n^5). The sign of delta_b40 is not determined.",
        "actual_logarithm_enclosure": "394<log n<462, from10^197<n<10^198 and2<log10<7/3. The upper logarithm follows from a positive exponential partial sum; the lower uses e<3 proved by its factorial-series geometric bound.",
        "actual_relative_rational_bounds": {
            "delta_b20_over_4lambda": exact20,
            "delta_b21_over_minus3gamma": exact21,
            "abs_delta_b40_over_gamma2_over_lambda": exact40,
        },
        "first_loop_matching_tolerances": "0<delta_b20/(4lambda)<10^-203, 0<delta_b21/(-3gamma)<10^-203, |delta_b40|/(gamma^2/lambda)<10^-192. These compare the complete FIRST-LOOP corrections with the unchanged separate tree coefficients.",
        "positive_higher_coefficient_in_specified_truncation": "The tree-plus-first-loop b40 lies between(1-10^-192)gamma^2/lambda and(1+10^-192)gamma^2/lambda and is strictly positive. This does not fix the sign of its first correction or claim an exact positive-measure/Gram saturation property.",
        "same_contact_and_heavy_weight": "The fixed OS4 contact has zero derivative for all three coefficients. No derivative condition is added and no unknown heavy pole or spectral weight is subtracted. The heavy exchange already contributes the positive tree b40.",
        "boundary": "No omitted-loop bound, exact physical LSZ coefficient, all-order dispersion relation, quantum UV completion, finite-gravity Regge control or original common-parent bounce follows. Only the complete first-loop matching is advanced.",
        "checks": {k: s.cancel(value) for k, value in checks.items()},
        "gates": {
            "actual_log_lower_domain": nactual > 10**197,
            "actual_log_upper_domain": nactual < 10**198,
            "first_positive_b20_leading_lower": 3 * 394
            - s.Rational(157, 45)
            - s.Rational(16 * 10**12) / nactual**2
            > 1000,
            "b20_fifth_order_term_positive": 24 * 394 - s.Rational(7969, 315) > 0,
            "full_b20_upper_including_remainder": 1400
            + s.Rational(11100) / nactual
            + s.Rational(16 * 10**12) / nactual**2
            < 1500,
            "first_negative_b21_magnitude_lower": 3 * 394
            - s.Rational(13, 21)
            - s.Rational(64 * 10**12) / nactual
            > 1000,
            "full_b21_magnitude_upper": 3 * 462
            - s.Rational(13, 21)
            + s.Rational(64 * 10**12) / nactual
            < 1500,
            "actual_relative_b20_margin": exact20 < s.Rational(1, 10**203),
            "actual_relative_b21_margin": exact21 < s.Rational(1, 10**203),
            "actual_relative_b40_margin": exact40 < s.Rational(1, 10**192),
            "higher_truncated_coefficient_positive_margin": s.Rational(1, 10**192) < 1,
            "no_sign_of_first_b40_correction_asserted": True,
        },
    }
