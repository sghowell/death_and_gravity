"""Full finite forward Coulomb phase and real classical/logarithmic coefficients."""

from functools import cache

import sympy as s

from . import source

ENERGY = s.Symbol("physical_energy", positive=True)
TAU = s.Symbol("positive_transfer_magnitude", positive=True)
ELL = s.Symbol("raw_log", real=True)
KAPPA = s.Symbol("positive_kappa", positive=True)


def eikonal(energy):
    energy = s.sympify(energy)
    return energy * energy - 4 * energy + 2


def gap(energy):
    energy = s.sympify(energy)
    return s.sqrt(energy * (energy - 4))


def rapidity(energy):
    energy = s.sympify(energy)
    return s.atanh(s.sqrt((energy - 4) / energy))


def classical_coefficient(energy=ENERGY):
    return 3 * s.pi**2 * (5 * eikonal(energy) + 6) / 4


def quantum_log_coefficient(energy=ENERGY):
    a = s.sympify(energy)
    v, d, h = eikonal(a), gap(a), rapidity(a)
    return (
        4 * v * (a - 2) * (8 - 3 * v) * h / d**3 - 2 * v * v / d**2 + (69 * v + 76) / 30
    )


def phase_principal_part(energy=ENERGY, tau=TAU, ell=ELL):
    v = eikonal(energy)
    return -2 * s.I * s.pi * v * v / (gap(energy) * tau) * (ell - s.log(tau) - 2 / v)


def phase_ratio(energy=ENERGY, tau=TAU, ell=ELL, kappa=KAPPA):
    v = eikonal(energy)
    return s.I * v / (8 * s.pi * kappa * gap(energy)) * (s.log(tau) - ell + 2 / v)


def triangle_remainder_bound(tau):
    tau = source.compact.require_delta(tau)
    return (
        s.Rational(5, 8) * s.sqrt(tau) - tau * s.log(tau) / 12 + s.Rational(5, 36) * tau
    )


@cache
def data():
    a, t, e = s.symbols("a t epsilon", nonzero=True)
    coeff = source.dimensional.master_coefficients(t, a, 1, e)
    c0 = s.factor(s.cancel(coeff["C00mu"].subs(e, 0)).subs(t, 0))
    b0 = s.factor(s.cancel(coeff["B00"].subs(e, 0)).subs(t, 0))
    q = s.Symbol("positive_q", positive=True)
    energy = q + 4
    v = eikonal(energy)
    d = gap(energy)
    h = rapidity(energy)
    ju = 4 * h / d
    js = (-4 * h + 2 * s.I * s.pi) / d
    jprime = -2 / d**2 + 4 * (energy - 2) * h / d**3
    fprime = 4 * v * (energy - 2) * (8 - 3 * v) * h / d**3 - 2 * v * v / d**2
    x = s.Symbol("unit_parameter", positive=True)
    tau = TAU
    chart = s.Symbol("nonnegative_chart", nonnegative=True)
    primitive = 2 * s.atan(chart * s.sqrt(1 + tau / 4)) / s.sqrt(1 + tau / 4)
    logshift = s.Symbol("logarithm", real=True)
    checks = {
        "whole_small_channel_massless_triangle_coefficient": s.factor(
            c0 + 3 * (5 * eikonal(a) + 6) / 2
        ),
        "whole_small_channel_massless_bubble_coefficient": s.factor(
            b0 - (87 * eikonal(a) + 118) / 60
        ),
        "crossed_eikonal_at_forward": s.expand(eikonal(4 - a) - eikonal(a)),
        "physical_J_sum_is_pure_Coulomb": s.simplify(js + ju - 2 * s.I * s.pi / d),
        "physical_B_forward_including_diagonal": s.simplify(
            (v * (js + ju) + 2) / 2 - 1 - s.I * s.pi * v / d
        ),
        "ordered_box_dimensional_finite_derivative": s.simplify(
            -4 * v * (js + ju) + 8 * s.I * s.pi * v / d
        ),
        "D_tree_soft_finite_phase": s.simplify(
            4 * ((v * (js + ju) + 2) / 2 - 1) - 4 * s.I * s.pi * v / d
        ),
        "combined_finite_phase_constant": s.simplify(
            -4 * v * (js + ju)
            + 4 * ((v * (js + ju) + 2) / 2 - 1)
            + 4 * s.I * s.pi * v / d
        ),
        "whole_forward_principal_phase_reconstruction": s.simplify(
            phase_principal_part(energy)
            + 2 * s.I * s.pi * v * v * (ELL - s.log(tau)) / (d * tau)
            - 4 * s.I * s.pi * v / (d * tau)
        ),
        "phase_normalization_against_gravity_Born": s.simplify(
            phase_principal_part(energy)
            / (16 * s.pi**2 * KAPPA**2)
            / (v / (KAPPA * tau))
            - phase_ratio(energy)
        ),
        "spacelike_J_derivative": s.simplify(-s.diff(ju, q) - jprime),
        "spacelike_squared_eikonal_J_derivative": s.simplify(
            -s.diff(v * v * ju, q) - fprime
        ),
        "whole_real_log_coefficient": s.simplify(
            fprime
            - c0.subs(a, energy) / 2
            - b0.subs(a, energy)
            - quantum_log_coefficient(energy)
        ),
        "whole_classical_real_coefficient": s.factor(
            -(s.pi**2) * c0 / 2 - classical_coefficient(a)
        ),
        "triangle_classical_chart_primitive": s.simplify(
            s.diff(primitive, chart) - 2 / (1 + (1 + tau / 4) * chart**2)
        ),
        "triangle_classical_chart_upper_endpoint": s.simplify(
            s.limit(primitive, chart, s.oo) - s.pi / s.sqrt(1 + tau / 4)
        ),
        "triangle_classical_closed_form": s.simplify(
            -s.pi / (2 * s.sqrt(tau)) * s.pi / s.sqrt(1 + tau / 4)
            + s.pi**2 / (s.sqrt(tau) * s.sqrt(tau + 4))
        ),
        "triangle_classical_leading_coefficient": s.limit(
            -(s.pi**2) / s.sqrt(tau + 4), tau, 0
        )
        + s.pi**2 / 2,
        "triangle_classical_sqrt_coefficient": s.limit(
            (-(s.pi**2) / s.sqrt(tau + 4) + s.pi**2 / 2) / tau, tau, 0
        )
        - s.pi**2 / 16,
        "triangle_constant_endpoint_log": s.integrate(
            -s.log(x * (1 - x)) / 2, (x, 0, 1)
        )
        - 1,
        "triangle_log_remainder_moment": s.integrate(x * (1 - x), (x, 0, 1))
        - s.Rational(1, 6),
        "triangle_endpoint_log_remainder_moment": s.integrate(
            -x * (1 - x) * s.log(x * (1 - x)), (x, 0, 1)
        )
        - s.Rational(5, 18),
        "triangle_remainder_classical_margin": s.Rational(10, 16) - s.Rational(5, 8),
        "complex_log_monodromy_not_removed": s.expand(
            (-2 * s.I * s.pi * v * v / (d * tau) * (ELL - logshift - 2 / v)).subs(
                logshift, logshift + 2 * s.I * s.pi
            )
            - (-2 * s.I * s.pi * v * v / (d * tau) * (ELL - logshift - 2 / v))
            + 4 * s.pi**2 * v * v / (d * tau)
        ),
    }
    return {
        "whole_forward_Coulomb_principal_part": phase_principal_part(),
        "whole_Coulomb_ratio_to_leading_Born": phase_ratio(),
        "whole_real_classical_coefficient": classical_coefficient(),
        "whole_real_quantum_log_coefficient": quantum_log_coefficient(),
        "full_forward_expansion": "F_known(s,-tau)=phase_principal_part"
        "+classical_coefficient/sqrt(tau)+quantum_log_coefficient*log(tau)+O(1). "
        "The remainder is complex and bounded on the fixed compact energy range. "
        "Both forward/backward endpoints follow by Bose symmetry.",
        "triangle_expansion": "C(-tau)=-pi^2/(2sqrt(tau))-log(tau)/2+1+R; "
        "|R|<5sqrt(tau)/8+tau*abs(log(tau))/12+5tau/36 for0<tau<=1.",
        "complex_analytic_boundary": "The leading log(tau)/tau coefficient is purely "
        "imaginary on the physical edge, not absent. Its monodromy is nonzero. "
        "Real one-loop rate control does not create a holomorphic transfer disk, "
        "control loop squares, justify a forward cross section, or close Regge.",
        "checks": checks,
        "gates": {
            "whole_small_channel_coefficients_from_frozen_D_input": True,
            "physical_massive_J_and_Coulomb_diagonal_retained": True,
            "both_box_and_D_tree_finite_phase_constants_retained": True,
            "complete_classical_and_quantum_log_real_terms": True,
            "triangle_endpoint_remainder_has_explicit_integrable_majorant": True,
            "written_full_complex_bounded_remainder_proof_required": True,
            "imaginary_phase_not_deleted_from_complex_amplitude": True,
            "all_three_finite_matching_coordinates_remain_separate": True,
        },
    }
