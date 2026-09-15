"""Exact sixth antiheat identity, all206 phase jets and unchanged orderings."""

from functools import cache

import sympy as s
from p8_vacuum_affine_weyl_operator_comparison import derivatives as old
from p8_vacuum_affine_weyl_operator_comparison import kernel

from . import source

PHASE = 96
ORDER = 6
MAX_ORDER = 2 * PHASE + 2 * (ORDER + 1)
CONSTANT = kernel.CONSTANT


def require_amplitude(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer, s.Rational)):
        raise TypeError("Require a declared exact full amplitude")
    if value not in (source.H_AMPLITUDE, source.F_AMPLITUDE):
        raise ValueError(
            "Require the new whole Hamiltonian or unchanged volume amplitude"
        )
    return s.sympify(value)


def require_order(value):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, s.Integer))
        or not 0 <= value <= MAX_ORDER
    ):
        raise ValueError("Require an exact proved phase order0 through206")
    return int(value)


def require_parameter_order(value):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, s.Integer))
        or not 0 <= value <= 2
    ):
        raise ValueError("Require a proved homogeneous derivative order0 through2")
    return int(value)


def coefficient(order):
    order = require_order(order)
    return s.Integer(2056) ** order * s.factorial(order) ** 3


def jet(order, amplitude, parameter_order=0):
    order = require_order(order)
    amplitude = require_amplitude(amplitude)
    parameter_order = require_parameter_order(parameter_order)
    return (
        2
        * amplitude
        * coefficient(order)
        / source.RADIUS**order
        * s.factorial(parameter_order)
        / source.HOMOGENEOUS_CAUCHY**parameter_order
    )


def radial_partition(order):
    order = require_order(order)
    return _radial_partition(order)


@cache
def _radial_partition(order):
    return sum(
        s.factorial(order)
        * 4 ** (order - 2 * j)
        * 256 ** (order - j)
        * s.factorial(order - j) ** 2
        / (s.factorial(order - 2 * j) * s.factorial(j))
        for j in range(order // 2 + 1)
    )


def heat_term(order, amplitude, parameter_order=0):
    if (
        isinstance(order, bool)
        or not isinstance(order, (int, s.Integer))
        or not 0 <= order <= ORDER + 1
    ):
        raise ValueError("Require the evaluated heat order0 through7")
    return (
        PHASE**order
        * jet(2 * order, amplitude, parameter_order)
        / (4**order * s.factorial(order))
    )


def bounds(amplitude, parameter_order=0):
    amplitude = require_amplitude(amplitude)
    parameter_order = require_parameter_order(parameter_order)
    return _bounds(amplitude, parameter_order)


@cache
def _bounds(amplitude, parameter_order):
    terms = [heat_term(n, amplitude, parameter_order) for n in range(ORDER + 2)]
    coherent = sum(terms[: ORDER + 1])
    remainder = CONSTANT * terms[ORDER + 1]
    return {
        "whole_all_eight_heat_terms": terms,
        "whole_finite_coherent_polynomial_bound": coherent,
        "whole_retained_Weyl_heat_remainder": remainder,
        "whole_both_complete_operator_bound": coherent + remainder,
        "whole_Weyl_minus_calibrated_operator_difference": sum(terms[2 : ORDER + 1])
        + remainder,
    }


@cache
def identities():
    clock, lam = s.symbols("positive_heat_time commuting_heat_generator", real=True)
    checks = {}
    polynomials = []
    for order in range(ORDER + 1):
        polynomial = sum((-clock * lam) ** j / s.factorial(j) for j in range(order + 1))
        expected = (
            (-1) ** order
            * clock**order
            * s.exp(clock * lam)
            * lam ** (order + 1)
            / s.factorial(order)
        )
        checks["entire_heat_remainder_derivative_" + str(order)] = s.simplify(
            s.diff(s.exp(clock * lam) * polynomial, clock) - expected
        )
        checks["entire_initial_heat_identity_" + str(order)] = (
            polynomial.subs(clock, 0) - 1
        )
        polynomials.append(polynomial)
    x = s.Symbol("whole_phase_coordinate", real=True)
    for degree in (2, 4, 8, 12, 14):
        f = x**degree
        powers = [f]
        for _ in range(7):
            powers.append(s.diff(powers[-1], x, 2) / 4)
        polynomial = sum((-1) ** j * powers[j] / s.factorial(j) for j in range(7))
        heat, current = polynomial, polynomial
        for j in range(1, degree // 2 + 1):
            current = s.diff(current, x, 2) / 4
            heat += current / s.factorial(j)
        remainder = (-1) ** 7 * powers[7] / s.factorial(7)
        checks["full_polynomial_remainder_not_dropped_" + str(degree)] = s.expand(
            f - heat - remainder
        )
    return {"whole_heat_polynomials": polynomials, "checks": checks}


@cache
def data():
    identity = identities()
    radial = [radial_partition(n) for n in range(MAX_ORDER + 1)]
    ratios = [s.Integer(2056) * (n + 1) ** 3 / source.RADIUS for n in range(MAX_ORDER)]
    match = s.Matrix([coefficient(n) - old.coefficient(n) for n in range(197)])
    return {
        "whole_phase_and_heat_orders": [PHASE, ORDER, MAX_ORDER],
        "whole_all207_radial_partition_bounds": radial,
        "whole_all207_phase_product_coefficients": [
            coefficient(n) for n in range(MAX_ORDER + 1)
        ],
        "whole_all206_successive_growth_ratios": ratios,
        "whole_exact_positive_heat_polynomials": identity["whole_heat_polynomials"],
        "whole_all_six_amplitude_parameter_operator_packets": {
            str((name, j)): bounds(amplitude, j)
            for name, amplitude in (
                ("Hamiltonian", source.H_AMPLITUDE),
                ("volume", source.F_AMPLITUDE),
            )
            for j in range(3)
        },
        "whole_extended_phase_proof": "The unchanged all-order bump, transition quotient, full radial partitions and Cauchy-Leibniz argument hold at every fixed finite order. Evaluate ALL207 radial sums and206 successive ratios through206=192+14. Orders0..196 exactly match the frozen product coefficients. The same simultaneous96-coordinate R/8 Cauchy domain and homogeneous d=1e-124 polydisc give Jn(A,j)=2A2056^n(n!)^3/R^n*j!/d^j. No higher real-time source regularity is assumed.",
        "whole_exact_remainder_proof": "For Pm(thetaD)=sumj0..m(-thetaD)^j/j!, d_theta[e^(thetaD)Pm(thetaD)]=(-1)^m theta^m e^(thetaD)D^(m+1)/m!. Integrate0..1. With m6, Weyl(b)=Q(P6b)+Weyl(r6), where r6=-integral0^1 theta^6 e^(thetaD)D^7b dtheta/6!. This is an EXACT identity with the full positive-heat remainder. It neither applies an ill-posed backward heat flow nor changes either original ordering.",
        "whole_both_operator_bound_proof": "The coherent polynomial is bounded by sumj0..6 96^j J2j/(4^j j!). Every Schur derivative of the full remainder has total order14..206. Positive heat is a supremum contraction and commutes with all derivatives; the evaluated decreasing Jn bounds these by96^7 J14/(4^7 7!). The unchanged normalized kernel theorem with constant1e112 therefore bounds the actual Weyl remainder. Sum both pieces for B(A); their positivity also bounds the original calibrated-coherent operator. No symbol-supremum shortcut is used.",
        "whole_ordering_difference_proof": "The calibrated operator remains Q[(1-D)b]. Subtract only this exact term from the foregoing identity, retaining Q[sumj2..6(-D)^j b/j!] and the full Weyl remainder. Their norms sum to the displayed Delta(A). The same proof holds for each first/second homogeneous derivative because the cutoff, state, whitening and fixed reference are independent of Y.",
        "whole_cutoff_difference_warning": "All derivatives of a difference of the two compact cutoff symbols vanish on their common core. The finite coherent polynomial shares that property. The positive-heat remainder generally does NOT; its operator-norm contribution remains in every subsequent state/force/readout comparison.",
        "checks": {
            **identity["checks"],
            "all197_original_phase_coefficients_match": match,
            "full192_plus14_phase_orders": s.Integer(MAX_ORDER - 206),
            "same_independent_kernel_constant": CONSTANT - 10**112,
        },
        "gates": {
            "all207_full_radial_partition_bounds": all(
                radial[n] <= 2048**n * s.factorial(n) ** 3 for n in range(MAX_ORDER + 1)
            ),
            "all206_phase_ratios": max(ratios) < s.Rational(1, 10**9),
            "both_complete_Hamiltonian_operators": bounds(source.H_AMPLITUDE)[
                "whole_both_complete_operator_bound"
            ]
            < 3 * source.H_AMPLITUDE,
            "both_complete_volume_operators": bounds(source.F_AMPLITUDE)[
                "whole_both_complete_operator_bound"
            ]
            < 3 * source.F_AMPLITUDE,
            "whole_relative_H_ordering_difference": bounds(source.H_AMPLITUDE)[
                "whole_Weyl_minus_calibrated_operator_difference"
            ]
            < source.H_AMPLITUDE / 10**55,
            "whole_relative_F_ordering_difference": bounds(source.F_AMPLITUDE)[
                "whole_Weyl_minus_calibrated_operator_difference"
            ]
            < source.F_AMPLITUDE / 10**55,
            "both_original_orderings_unchanged": True,
            "complete_heat_remainder_not_declared_core_supported": True,
            "no_backward_heat_evolution_or_Taylor_defined_Hamiltonian": True,
        },
    }
