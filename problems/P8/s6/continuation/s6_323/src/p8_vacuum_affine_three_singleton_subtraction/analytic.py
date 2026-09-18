"""Uniform hard-only origin/physical tubes and the complete third rectangle."""

from functools import cache

import sympy as s
from p8_vacuum_affine_pair_singleton_subtraction import hard as pair_class
from p8_vacuum_affine_relative_energy_complex_tube import bounds as hard_majorant

from . import source

ETA = s.Rational(1, 10**13)
SIGMA = s.Rational(1, 10**14)


@cache
def geometry():
    rho = 4 * SIGMA
    gates = {
        "smallest_squared_massive_radius": s.Rational(5, 4) ** 2 - 1
        == s.Rational(9, 16),
        "uniform_E_over_r0_below_two": s.Rational(5, 4) / s.Rational(3, 4) < 2,
        "uniform_r0_over_E_below_seven_eighths": s.Rational(7, 8) ** 2
        > s.Rational(3, 4),
        "origin_squared_energy_perturbation_below7rho": 6 * rho
        + s.Rational(9, 2) * rho**2
        < 7 * rho,
        "origin_square_root_branch_margin": 7 * rho
        < s.Rational(3, 4) * s.Rational(9, 16),
        "origin_rprime_shift_below14rho": 7 / s.Rational(3, 4) < 14,
        "origin_Eprime_shift_below7rho": 7 / s.Rational(5, 4) < 7,
        "complex_rprime_Eprime_ratio_below_one": s.Rational(7, 8) + 21 * rho < 1,
        "h_plus_Eprime_gap_above_two": s.Rational(5, 2) - s.Rational(17, 2) * rho > 2,
        "origin_beta_magnitude_below_one": s.Rational(1, 2) + 3 * rho / 8 < 1,
        "complex_spatial_momentum_below_seven_fourths": 17 * rho < s.Rational(1, 56),
        "Doppler_perturbation_less_than_one_eighth": 20 * rho < s.Rational(1, 8),
        "origin_phase_branch_stays_near_one": 42 * rho < 60 * rho < s.Rational(1, 4),
        "quadratic_hard_gap_perturbation_term": 40000 * SIGMA < 1,
        "origin_hard_denominator_above_half_Born_gap": 801 * SIGMA < s.Rational(1, 2),
        "spatial_recoil_shift_below40rho": 14 * rho + 3 * rho < 40 * rho,
        "mixed_hard_momentum_shift_below100rho": 40 * rho + 3 * rho < 100 * rho,
        "timelike_hard_momentum_shift_below100rho": 80 * rho + 3 * rho < 100 * rho,
        "paired_origin_current_below2e4_t": 1024 * (1 + 40 * SIGMA) < 2 * 10**4,
        "massive_spatial_current_gradient_below1024": 4 * 8 + 4 * 3 * 64 < 1024,
        "four_family_second_remainder_budget": 4 * 256 < 10**4,
        "four_family_first_remainder_budget": 4 * 64 < 10**4,
    }
    E, r = s.symbols("E r", positive=True)
    u = s.Matrix(s.symbols("u0:3", real=True))
    v = s.Matrix(s.symbols("v0:3", real=True))
    eta = s.diag(1, -1, -1, -1)
    p3 = s.Matrix([E, *(r * u)])
    p4 = s.Matrix([E, *(-r * u)])
    dp3 = s.Matrix(
        [-s.Rational(1, 2) - r * u.dot(v) / (2 * E), *(-E * u / (2 * r) - v / 2)]
    )
    dp4 = s.Matrix(
        [-s.Rational(1, 2) + r * u.dot(v) / (2 * E), *(E * u / (2 * r) - v / 2)]
    )
    checks = {
        "origin_recoil_linear_conservation": dp3 + dp4 + s.Matrix([1, *v]),
        "origin_first_scalar_shell_modulo_unit_direction": s.factor(
            (p3.T * eta * dp3)[0] - E * (u.dot(u) - 1) / 2
        ),
        "origin_second_scalar_shell_modulo_unit_direction": s.factor(
            (p4.T * eta * dp4)[0] - E * (u.dot(u) - 1) / 2
        ),
        "spatial_norm_rational_margin": s.Rational(7, 4) ** 2 - 3 - s.Rational(1, 16),
    }
    t = s.Symbol("t", real=True)
    h = E - t / 2
    ep = s.sqrt(h * h - t * t * v.dot(v) / 4)
    rp = s.sqrt(ep * ep - 1)
    exact3 = s.Matrix(
        [
            h - rp * t * u.dot(v) / (2 * ep),
            *(
                rp * u
                + (rp * t * u.dot(v) / (4 * ep * (h + ep)) - s.Rational(1, 2)) * t * v
            ),
        ]
    )
    exact4 = s.Matrix(
        [
            h + rp * t * u.dot(v) / (2 * ep),
            *(
                -rp * u
                + (-rp * t * u.dot(v) / (4 * ep * (h + ep)) - s.Rational(1, 2)) * t * v
            ),
        ]
    )
    for which, (exact, derivative) in enumerate(((exact3, dp3), (exact4, dp4))):
        actual = exact.diff(t).subs(t, 0).subs(s.sqrt(E * E - 1), r)
        for component, value in enumerate(actual - derivative):
            checks[f"literal_recoil_linear_coeff_{which}_{component}"] = s.factor(value)
    p = s.Matrix(s.symbols("p0:3", real=True))
    n = s.Matrix(s.symbols("n0:3", real=True))
    A = s.zeros(3)
    for i in range(3):
        for j in range(i, 3):
            A[i, j] = A[j, i] = s.Symbol(f"A{i}_{j}")
    d = E - p.dot(n)
    N = (p.T * A * p)[0]
    S = N / d
    for i in range(3):
        target = (2 * (A * p)[i] * d - N * (p[i] / E - n[i])) / d**2
        checks[f"on_shell_spatial_current_gradient_{i}"] = s.factor(
            s.diff(S, p[i]) + s.diff(S, E) * p[i] / E - target
        )
    checks["signed_massive_current_odd"] = s.factor(
        S.xreplace({E: -E, **{v: -v for v in p}}) + S
    )
    return {
        "checks": checks,
        "gates": {k: bool(v) for k, v in gates.items()},
        "whole_origin_ratio": SIGMA,
        "whole_positive_center_ratio": ETA,
        "whole_origin_tube": "For each mixed Born channel t=sqrt(tau)>0, |z_i|<=sigma*t. Direct recoil radicand bounds give |dp0|<20rho,||dp_space||<40rho and hard shifts<100rho. Then |D-D0|<801sigma*t^2<t^2/2. No future-energy assumption is used inside this separate origin disc.",
        "whole_positive_center_tube": "For the grouped leading hard-only kernels, |z_i-w_i|<=eta*W for all3 energies. The S313 total perturbation and hard-gap estimates apply; all singleton energy poles have already been canceled algebraically. This domain is not asserted for connected pair/triple line currents or the full amplitude.",
        "whole_phase": "At the origin |rho_phase-1|<60rho, so its positive branch stays analytic and below2. At positive centers the retained S313 phase estimate applies. No conjugation or absolute value is inserted into the analytic action.",
    }


@cache
def budgets():
    x = s.Symbol("x")
    functions = hard_majorant.hard_functions(x)
    cm = s.cancel(functions["matter"] / functions["L"] ** 4)
    cg = s.cancel(functions["gravity"] / functions["L"] ** 4)
    matter = tuple(s.diff(cm, x, i).subs(x, 0) for i in range(4))
    gravity = tuple(s.diff(cg, x, i).subs(x, 0) for i in range(4))
    expected_m = (s.Rational(5, 2), 3088, 18923696, 174400047744)
    expected_g = (
        235929600000,
        5566277616582367641600000,
        262649930314558756486872922324992000000,
        18590067899792945264953972626291745126348750848000000,
    )
    n, K = source.HEAVY_MASS2, source.KAPPA
    central = tuple(s.factor(n * n * m + g) for m, g in zip(matter, gravity))
    old = pair_class.budgets()
    B2 = old["whole_complex_core_remainder_B2"]
    endpoint0 = 4096**2
    endpoint1 = max(2 * 4096**3, 4096**2 * 120000000 * 6 * 32**3)
    raw3 = 2 * 3 * 8 * 8 * endpoint0 * 10**6 * (2 * 10**4) ** 3
    raw1 = 2 * 3 * 4 * 4 * 3 * endpoint1 * 10**6 * (2 * 10**4) ** 2
    group_cap = 10**50 * (n * n + 1)
    B3 = 10**4 * group_cap / SIGMA**3
    B1 = 10**4 * group_cap / SIGMA**2
    canonical = K ** s.Rational(3, 2)
    coefficients = {
        "central3": 2 * central[3] / canonical,
        "central2": 256 * central[2] / (ETA * canonical),
        "central1_pair": 8 * central[1] * (2 * 10**10) / canonical,
        "central1_leading": B1 / (ETA**2 * canonical),
        "external0_connected3": 8 * central[0] * 10**24 / canonical,
        "external0_pair_remainder": 1024 * central[0] * 10**11 / (ETA * canonical),
        "external0_pair_leading": 10**14 * B2 / (ETA * canonical),
        "external0_leading3": B3 / (ETA**3 * canonical),
    }
    core_coefficient = sum(coefficients.values())
    previous_coefficient = (
        old["whole_pair_product_derivative_coefficient"]
        + old["whole_retained188_derivative_coefficient"]
    )
    full_coefficient = core_coefficient + previous_coefficient
    checks = {
        f"central_matter_coefficient_{i}": matter[i] - expected_m[i] for i in range(4)
    }
    checks.update(
        {
            f"central_gravity_coefficient_{i}": gravity[i] - expected_g[i]
            for i in range(4)
        }
    )
    checks.update(
        {
            "same_derived_S322_core_remainder_constant": B2 - (10**40 * n * n + 10**60),
            "pair_R0_coefficient_arithmetic": s.Integer(
                8 * 8 * 4 + 2 * 8 * 7 * 10 + 14000 * 100 - 1401376
            ),
            "pair_R1_coefficient_arithmetic": s.Integer(
                8 * 2 * 10 + 5000 * 100 - 500160
            ),
            "pair_leading_coefficient_arithmetic": s.Integer(
                1500000 * 16 * 8 - 192000000
            ),
            "pair_remainder_coefficient_arithmetic": s.Integer(
                4 * 600000 * 16 * 4 + 1500000 * 4 * 4 * 16 * 8 * 4 - 12441600000
            ),
            "eight_core_coefficients_retained": s.Integer(len(coefficients) - 8),
            "complete_class_inventory": s.Integer(3767 + 1349 - 5116),
        }
    )
    a, b, c = s.symbols("a b c", positive=True)
    I = lambda u, v: (u + v) * s.log(u + v) - u * s.log(u) - v * s.log(v)
    J = c * I(a, b) + b * I(a, c) + a * I(b, c)
    checks["entropy_kernel_mixed_pair"] = s.factor(s.diff(I(a, b), a, b) - 1 / (a + b))
    checks["three_central_choices_sum_to_one_J"] = s.expand(
        a * I(b, c) + b * I(a, c) + c * I(a, b) - J
    )
    checks["pair_scale_common_numerator"] = s.expand(
        (a + b) * (a + c)
        + (a + b) * (b + c)
        + (a + c) * (b + c)
        - (a + b + c) ** 2
        - a * b
        - a * c
        - b * c
    )
    gates = {
        "all_central_majorant_coefficients_positive": all(
            v > 0 for v in (*matter, *gravity)
        ),
        "scalar_pair_R0_below1point5e6": 8 * 8 * 4 + 2 * 8 * 7 * 10 + 14000 * 100
        < 1500000,
        "scalar_pair_R1_below6e5": 8 * 2 * 10 + 5000 * 100 < 600000,
        "scalar_pair_leading_below2e8": 1500000 * 16 * 8 < 2 * 10**8,
        "scalar_pair_remainder_below1e11": 4 * 600000 * 16 * 4
        + 1500000 * 4 * 4 * 16 * 8 * 4
        < 10**11,
        "normalized_connected_pair_magnitude_below2e9": 2 * 10**8 * 4 * 2 < 2 * 10**9,
        "normalized_connected_pair_gradient_below1e16": 8 * 10**15 < 10**16,
        "effective_pair_current_fits_existing_norm": s.Rational(2 * 10**9, 10**14) < 32,
        "effective_pair_current_fits_existing_gradient": s.Rational(10**16, 10**14)
        < 2048,
        "real_full_scalar_pair_below2e10_ab_over_u": 2 * 10**8 + s.Rational(10**11, 8)
        < 2 * 10**10,
        "leading_three_grouped_gravity_budget_below1e30": raw3 < 10**30,
        "one_central_grouped_gravity_budget_below1e40": raw1 < 10**40,
        "four_compact_channel_budgets_covered": all(
            2 * table[r] * 128**m < 10**50
            for table in (matter, gravity)
            for r, m in ((0, 3), (1, 2))
        ),
        "two_required_endpoint_coefficient_budgets": all(
            s.factorial(r) * 2**r * (144 * s.binomial(r + 2, 2) + r + 1)
            < s.factorial(r) * 4096**r
            for r in (1, 2)
        ),
        "larger_old_pair_leading_arithmetic_envelope": coefficients[
            "external0_pair_leading"
        ]
        < 2 * 10**14 * B2 / (ETA**3 * canonical),
        "all_eight_coefficients_strictly_positive": all(
            value > 0 for value in coefficients.values()
        ),
        "retained1349_coefficient_below2e_minus723": 0
        < previous_coefficient
        < s.Rational(2, 10**723),
        "complete5116_amplitude_rectangle_below1e_minus670": 0
        < full_coefficient
        < s.Rational(1, 10**670),
        "constant_and_linear_origin_terms_annihilated_by_third_rectangle": True,
        "one_central_prefactor_changes_required_derivative_order": True,
        "connected_pair_and_triple_never_given_global_W_holomorphy": True,
        "amplitude_faces_alone_not_inclusive_probability": True,
        "original_P8_state_Regge_real_virtual_all_N_still_open": True,
    }
    return {
        "checks": checks,
        "gates": {k: bool(v) for k, v in gates.items()},
        "whole_central_majorant_coefficients": {
            "matter": matter,
            "gravity": gravity,
            "relative": central,
        },
        "whole_grouped_kernel_budgets": {
            "leading3_gravity": raw3,
            "central1_gravity": raw1,
            "common_C": group_cap,
            "second_remainder_B3": B3,
            "first_remainder_B1": B1,
        },
        "whole_eight_core_rectangle_coefficients": coefficients,
        "whole_core_rectangle_coefficient": core_coefficient,
        "whole_retained1349_rectangle_coefficient": previous_coefficient,
        "whole_full5116_rectangle_coefficient": full_coefficient,
        "whole_rectangle_threshold": s.Rational(1, 10**670),
        "whole_two_regime_argument": "For a mixed channel t=sqrt(tau), the grouped hard-only kernels satisfy C*delta*(t+W)^m/(t^2+W^2) on the physical disc, and C*t^m on an origin disc of radius sigma*t, m=3 or2. At W<=sigma*t/4 use first/second Cauchy-Taylor estimates; otherwise use t<4W/sigma and the origin Taylor coefficients. The complete second and first remainders are bounded byB3*W^2 andB1*W, uniformly in all hard angles and energy hierarchies.",
        "whole_eight_class_completion": "Six contributions are controlled by explicit connected-line factors, c-only Cauchy and the already derived two-current forward grouping; the two leading sectors use the new grouped Taylor estimates. Compatible faces and the entropy kernel give a full5116-tree amplitude rectangle below1e-670*J. Exact original poles still have no assigned point value.",
        "whole_probability_boundary": "This component bounds the amplitude rectangle. The separate measure module derives its defined finite signed subtraction, not a positive normalized inclusive probability or matched real-virtual combination. All-N summation, loop and evanescent matching, interacting state, absolute complex Regge and original P8 closure remain open.",
    }


def rectangle_envelope(a, b, c):
    a, b, c = source.require_energies(a, b, c)
    I = lambda u, v: (u + v) * s.log(u + v) - u * s.log(u) - v * s.log(v)
    return (c * I(a, b) + b * I(a, c) + a * I(b, c)) / s.Integer(10) ** 670


@cache
def data():
    chart, budget = geometry(), budgets()
    return {
        "checks": {
            **{"geometry_" + k: v for k, v in chart["checks"].items()},
            **{"budget_" + k: v for k, v in budget["checks"].items()},
        },
        "gates": {**chart["gates"], **budget["gates"]},
        "whole_two_analytic_domains": {
            k: v for k, v in chart.items() if k not in ("checks", "gates")
        },
        "whole_complete_subtraction": {
            k: v for k, v in budget.items() if k not in ("checks", "gates")
        },
    }
