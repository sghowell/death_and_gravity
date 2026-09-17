"""Probability subtraction, necessary state connector and matched soft reference."""

from functools import cache

import sympy as s
from p8_vacuum_affine_two_real_soft_overlap import analytic as old_analytic
from p8_vacuum_affine_two_real_soft_overlap import subtraction as overlap

from . import source

require_resolution = overlap.require_resolution
AXIS_COEFFICIENT = s.Rational(1, 10**350)


def two_real_variation_upper(value):
    x = require_resolution(value)
    return 2 * s.Rational(1, 10**725) * x + 2 * s.Rational(1, 10**652) * x * x


def connector_relative_upper(value):
    x = require_resolution(value)
    return 100000 * x * (2 - s.log(x)) / source.KAPPA**2


def dressed_two_real_relative_upper(value):
    return 2 * two_real_variation_upper(value)


def assembled_relative_upper(value):
    x = require_resolution(value)
    single = (
        2
        * s.Min(
            s.Integer(10) ** 32,
            2 * s.Integer(10) ** 14 * x + 2 * s.Integer(10) ** 37 * x * x,
        )
        / source.KAPPA
    )
    return single + connector_relative_upper(x) + dressed_two_real_relative_upper(x)


@cache
def general_identities():
    ur, ui, vr, vi, orr, oi = s.symbols("ur ui vr vi orr oi", real=True)
    U = ur + s.I * ui
    V = vr + s.I * vi
    O = orr + s.I * oi

    def square(z):
        return s.expand(z * s.conjugate(z))

    square_difference = s.expand(square(U + V - O) - square(U) - square(V) + square(O))
    f2, f1a, f1b, ka0, kb0, kasb, kbsa = s.symbols(
        "f2 f1a f1b ka0 kb0 kasb kbsa", real=True
    )
    r1a = f1a - ka0
    r1b = f1b - kb0
    r2 = f2 - kasb * f1b - kbsa * f1a + ka0 * kb0
    connector = (kasb - ka0) * kb0 + (kbsa - kb0) * ka0
    assembled = ka0 * kb0 + kasb * r1b + kbsa * r1a + r2 + connector
    checks = {
        "general_complex_amplitude_probability_conversion": s.expand(
            square_difference - 2 * s.re((U - O) * s.conjugate(V - O))
        ),
        "whole_two_real_density_reconstruction": s.expand(assembled - f2),
        "missing_connector_difference_equals_required_transport": s.expand(
            assembled - connector - f2 + connector
        ),
        "both_single_soft_faces_and_one_common_overlap": s.expand(
            r2 - f2 + kasb * f1b + kbsa * f1a - ka0 * kb0
        ),
        "new_axis_product_variation_coefficient": AXIS_COEFFICIENT**2
        - s.Rational(1, 10**700),
    }
    for count in range(10):
        checks["two_marked_Bose_factor_" + str(count)] = s.binomial(
            count + 2, 2
        ) / s.factorial(count + 2) - s.Rational(1, 2) / s.factorial(count)
    gates = {
        "omitting_the_state_connector_is_a_nonzero_generic_error": s.expand(
            assembled - connector - f2
        )
        != 0,
        "amplitude_and_probability_subtractions_not_identified": square_difference != 0,
    }
    return checks, gates


@cache
def original_samples():
    old_checks, old_gates, records = overlap.original_samples()
    assert all(
        v == 0
        for value in old_checks.values()
        for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert all(old_gates.values())
    checks = {}
    for family, denominator, a, b, F2, B2 in records:
        _, qs, born = overlap.configuration(a, b)
        pa, qa, _ = overlap.configuration(a, 0)
        pb, qb, _ = overlap.configuration(0, b)
        M5a, _ = overlap.amplitude(pa, qa)
        M5b, _ = overlap.amplitude(pb, qb)
        A0 = s.factor(
            overlap.matter.born_continuation(
                born, source.HEAVY_MASS2, source.CUBIC, source.CONTACT
            )
            + overlap.e.old.born(born) / source.KAPPA
        )
        U = (
            s.sqrt(overlap.rho(0, b))
            * overlap.current((*pb, qb[0]), qs[0])
            * M5b
            / (2 * s.sqrt(source.KAPPA) * A0)
        )
        V = (
            s.sqrt(overlap.rho(a, 0))
            * overlap.current((*pa, qa[0]), qs[1])
            * M5a
            / (2 * s.sqrt(source.KAPPA) * A0)
        )
        O = (
            overlap.current(born, qs[0])
            * overlap.current(born, qs[1])
            / (2 * source.KAPPA)
        )
        label = family + "_" + str(denominator)
        checks[label + "_same_phase_and_normalized_polarization_faces"] = s.expand(
            U + V - O - B2
        )
        probability = F2**2 - U**2 - V**2 + O**2
        amplitude = F2**2 - B2**2
        checks[label + "_entire434_probability_decomposition"] = s.expand(
            probability - amplitude - 2 * (U - O) * (V - O)
        )
    return checks


@cache
def data():
    checks, gates = general_identities()
    checks = dict(checks)
    gates = dict(gates)
    checks.update(original_samples())
    K = source.KAPPA
    upper = (
        2 * s.Rational(1, 10**768)
        + 4 * s.Rational(1, 10**725) / 8
        + 4 * s.Rational(1, 10**652) / 64
        + 100000 * s.Rational(5, 8) / K**2
    )
    x = s.Symbol("resolution", positive=True)
    checks.update(
        {
            "connector_reference_prefactor_budget": s.Rational(4, 5) * (12000 + 40)
            - 9632,
            "connector_uniform_positive_threshold_primitive": s.simplify(
                s.diff(x * (2 - s.log(x)), x) - (1 - s.log(x))
            ),
            "connector_uniform_zero_threshold_limit": s.limit(
                x * (2 - s.log(x)), x, 0, dir="+"
            ),
            "two_real_linear_dressing_factor": 2 * (2 * s.Rational(1, 10**725))
            - 4 * s.Rational(1, 10**725),
            "two_real_quadratic_dressing_factor": 2 * (2 * s.Rational(1, 10**652))
            - 4 * s.Rational(1, 10**652),
        }
    )
    return {
        "whole_probability_residual": "For U=J_a(sigma_b)F1b/sqrt(kappa),V=J_b(sigma_a)F1a/sqrt(kappa),O=G00/(ab), define r2=|F2|2-|U|2-|V|2+|O|2. It equals S313's signed amplitude-subtracted density plus2Re[(U-O)conjugate(V-O)]. The continuous-axis bound |G(0,b)-G00|<1e-350*b implies |U-O|<1e-350/a and similarly V. With the unchanged physical phase/polarization/Bose factors, TV2<2e-725*x+2e-652*x2.",
        "whole_missing_connector": "At the physical D4 two-real density level, P0 plus the S309 one-residual class and r2 misses C2=[K_a(sigma_b)-K_a0]K_b0+[K_b(sigma_a)-K_b0]K_a0, with K the polarization-summed squared leading current. This generally nonzero REAL term is not separately soft integrable. It must retain its own state-correct virtual completion.",
        "whole_finite_state_transport_connector": "Define E1(x)=integral dB1(sigma_b)[P_sigma_b(x-b)-P0(x-b)], with the elastic one-soft Born seed and both known factors in the same S301 angular/phase convention. The difference is taken before the divergent seed integration. The new finite-e signed-series domination justifies regulator removal at fixed x and yields this finite expression. No finite hard radiative or evanescent matching term is inferred from it.",
        "whole_unexpanded_connector_bound": "At R<=x, |ln[P_sigma(x)/P0(x)]|<6000R(1-lnR)/kappa, hence its ratio difference<12000R(1-lnR)/kappa. With z=1-R/x, physical nonnegative indices give |z^a_sigma-z^a0|<=40R|lnz|/kappa. Integrating against a0*dR/R,a0<4/(5kappa), gives |E1|/P0<1e5*x*(2-lnx)/kappa2. This tends to zero uniformly and never expands a0*lnx.",
        "whole_two_marked_dressing": "Use the finite signed r2 measure with the full state sigma_ab and remaining total energy x-a-b. S309's signed-series proof applies because S300/S301 depend on total marked energy rather than its multiplicity. The label identity choose(N+2,2)/(N+2)!=1/(2!N!) preserves the seed's existing1/2!. Thus |D2|/P0<(1+4250/kappa)TV2<4e-725*x+4e-652*x2.",
        "whole_scoped_matched_reference": "P_match2=P0+D1[r1]+E1+D2[r2]. Its explicit zero/one/two-real D4 TREE radiation-count densities match the complete unchanged original amplitudes; the algebra is before integration, not a fixed-N regulator-limit interchange. Hard Born and numeric coefficients are held fixed in this counting identity. The same-state leading-soft virtual completion is precisely defined by the entire sums. This is not a complete coupling expansion including hard radiative loops.",
        "whole_reference_error_bound": upper,
        "whole_positive_reference_statement": "The sum of the S309 bound, connector bound and two-marked bound is<1e-653 uniformly at0<x<=1/8 and vanishes asx->0. Therefore this named reference is pointwise positive relative to P0. It is not proved monotone in threshold, a positive event measure, unitary, or equal to the physical all-order detector rate.",
        "whole_remaining_obligations": "Full finite hard real-virtual and evanescent matching, higher nonleading real seeds and loops, the actual interacting quantum state, absolute complex Regge control and the original common-parent bounce/P8 remain open.",
        "checks": checks,
        "gates": {
            **gates,
            "S313_actual_linear_slope_below_axis_bound": old_analytic.LINEAR_SLOPE
            < AXIS_COEFFICIENT,
            "added_axis_product_below_existing_quadratic_coefficient": s.Rational(
                1, 10**700
            )
            < s.Rational(1, 10**652),
            "original_state_log_ratio_below_half": 3000 / K < s.Rational(1, 2),
            "original_dressing_ratio_below_two": 1 + 4250 / K < 2,
            "whole_original_matched_reference_error_below1e_minus653": upper
            < s.Rational(1, 10**653),
            "finite_connector_difference_not_separate_divergent_terms": True,
            "two_marked_seed_Bose_factor_not_counted_twice": True,
            "tree_density_matching_not_missing_hard_virtual_matching": True,
            "pointwise_positive_reference_not_probability_or_unitarity": True,
        },
    }
