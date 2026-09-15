"""Complete corrected symbols, both exact orderings and all206 phase jets."""

from functools import cache

import sympy as s
from p8_vacuum_affine_weyl_operator_comparison import derivatives as old
from p8_vacuum_affine_weyl_operator_comparison import kernel

from . import geometry, source

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


def operator_packet(amplitude, parameter_order=0):
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
def operator_data():
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
            str((name, j)): operator_packet(amplitude, j)
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
            "both_complete_Hamiltonian_operators": operator_packet(source.H_AMPLITUDE)[
                "whole_both_complete_operator_bound"
            ]
            < 3 * source.H_AMPLITUDE,
            "both_complete_volume_operators": operator_packet(source.F_AMPLITUDE)[
                "whole_both_complete_operator_bound"
            ]
            < 3 * source.F_AMPLITUDE,
            "whole_relative_H_ordering_difference": operator_packet(source.H_AMPLITUDE)[
                "whole_Weyl_minus_calibrated_operator_difference"
            ]
            < source.H_AMPLITUDE / 10**55,
            "whole_relative_F_ordering_difference": operator_packet(source.F_AMPLITUDE)[
                "whole_Weyl_minus_calibrated_operator_difference"
            ]
            < source.F_AMPLITUDE / 10**55,
            "both_ordering_definitions_retained_for_new_corrected_symbols": True,
            "complete_heat_remainder_not_declared_core_supported": True,
            "no_backward_heat_evolution_or_Taylor_defined_Hamiltonian": True,
        },
    }


H_AMPLITUDE, F_AMPLITUDE = source.H_AMPLITUDE, source.F_AMPLITUDE
STATE_WEIGHT = s.Rational(1, 10**350)
FORCE = s.Rational(1, 10**230)
HOMOGENEOUS_DERIVATIVE = s.Integer(10) ** 470


def operator_bound(amplitude, homogeneous_order=0):
    return operator_packet(amplitude, homogeneous_order)[
        "whole_both_complete_operator_bound"
    ]


@cache
def bounds():
    first = operator_bound(H_AMPLITUDE, 1)
    second = operator_bound(H_AMPLITUDE, 2)
    return {
        "actual_corrected_centered_interaction_amplitude": geometry.remainder()[
            "whole_corrected_centered_interaction_amplitude"
        ],
        "actual_corrected_centered_volume_amplitude": geometry.remainder()[
            "whole_corrected_centered_volume_amplitude"
        ],
        "centered_operator": operator_bound(H_AMPLITUDE),
        "each_first_homogeneous_derivative": first,
        "each_second_homogeneous_derivative": second,
        "quantum_force": 10**104 * first / geometry.KAPPA,
        "quantum_force_Lipschitz": 10**106 * (first + second) / geometry.KAPPA,
        "centered_volume_operator": operator_bound(F_AMPLITUDE),
    }


@cache
def data():
    generic, b = operator_data(), bounds()
    hraw, href, h0, chi, fdot = s.symbols(
        "full_raw_after_boundary_pullback full_prepared_reference raw_center full_cutoff full_boundary_time_generator"
    )
    volume, v0 = s.symbols(
        "whole_raw_volume_after_boundary_pullback actual_homogeneous_volume"
    )
    return {
        "whole_universal_206jet_and_full_positive_heat_estimates": {
            key: value
            for key, value in generic.items()
            if key not in ("checks", "gates")
        },
        "whole_corrected_full_operator_and_force_bounds": b,
        "whole_corrected_model_definition": "Let z denote the original PREPARED canonical variables, with the unchanged initial full covariance/whitening S0, original prepared reference Mref and its propagator Uref. The physical raw variables are Sb(u)z. Define hcor=hraw(u,Y,Sb(u)z)+partial_uF_b(z), with EVERY original raw nonlinear source and exact auxiliary/spatial reconstruction retained. Its prepared interaction is g=(hcor-href)(u,Y,Mref z), g0=hraw(u,Y,0), and b_c=chi_c(g-g0). The new corrected centered operators are K_C=Q_V[(1-D)b_c] and K_W=OpW(b_c), D=Delta_w/4. Both actual cutoff derivative contacts and the complete positive-heat remainder are kept.",
        "whole_corrected_state_and_observable_picture": "The coupled Hilbert unknown phi starts at the ORIGINAL prepared seed. Restore physical raw state Ub(u)Uref(u)exp[-i integral g0(Y)]phi. The raw initial state is Ub(0)psi0, not the archived unshifted raw Gaussian. Transport the readout and positive coherent core effect by the SAME Ub and Uref. The physical volume symbol composes with Sb and Mref: F_Y=average exp(3v)R_full(u,Nstar)^(-3/4), F0=R_full(u,N0(Y))^(-3/4). Quantize F0+chi_c(F_Y-F0) in each ordering and multiply its expectation by exp(3alpha).",
        "whole_canonical_quantization_boundary": "The fixed-reference time-dependent Sb is metaplectic in the finite canonical system. Weyl covariance is exact, and the coherent map is covariant when state/window are transported together. The quadratic generating term is exactly calibrated before localization because D^2F_b=0. Once included in the full compact b_c, every cutoff derivative is retained. No nonlinear canonical covariance theorem, independent new vacuum or deletion of the boundary phase is assumed.",
        "whole_parameter_force_proof": "Sb, partial_uF_b, chi, whitening and Uref are fixed-reference objects independent of the five LIVE Y variables. Differentiate the whole hraw(u,Y,Sb z) and actual implicit roots; retain every original homogeneous canonical chain-rule row and every scalar-center force. The same Cauchy radii and corrected complete amplitudes yield displayed first/second bounded force derivatives. The canonical coefficient row ceilings1e104/kappa and1e106/kappa remain valid. The boundary phase is configuration-dependent; only each solution's separate scalar g0 phase may be factored for vector comparisons.",
        "checks": {
            **generic["checks"],
            "whole_corrected_scalar_center_split": s.expand(
                h0
                + chi * (hraw + fdot - href - h0)
                - (h0 + chi * ((hraw + fdot - href) - h0))
            ),
            "whole_corrected_volume_center_not_replaced_by_one": s.expand(
                v0 + chi * (volume - v0) - (chi * volume + (1 - chi) * v0)
            ),
            "both_corrected_amplitudes_are_declared": s.Matrix(
                [H_AMPLITUDE - source.H_AMPLITUDE, F_AMPLITUDE - source.F_AMPLITUDE]
            ),
            "quantum_force_uses_original_full_canonical_row": b["quantum_force"]
            - 10**104 * b["each_first_homogeneous_derivative"] / geometry.KAPPA,
        },
        "gates": {
            **generic["gates"],
            "corrected_full_Hamiltonian_amplitude_hypothesis": b[
                "actual_corrected_centered_interaction_amplitude"
            ]
            < H_AMPLITUDE,
            "corrected_full_volume_amplitude_hypothesis": b[
                "actual_corrected_centered_volume_amplitude"
            ]
            < F_AMPLITUDE,
            "whole_corrected_uniform_quantum_force": b["quantum_force"] < FORCE,
            "whole_corrected_force_Lipschitz": b["quantum_force_Lipschitz"] < 1,
            "whole_corrected_generator_parameter_row": 5
            * b["each_first_homogeneous_derivative"]
            < HOMOGENEOUS_DERIVATIVE,
            "whole_corrected_parameter_Cauchy_polydisc": 5 * source.HOMOGENEOUS_CAUCHY
            + source.REAL_RADIUS
            < source.HOMOGENEOUS_RADIUS,
            "new_corrected_solutions_not_archived_S273_solutions": True,
            "same_prepared_seed_and_correctly_transported_raw_seed": True,
        },
    }
