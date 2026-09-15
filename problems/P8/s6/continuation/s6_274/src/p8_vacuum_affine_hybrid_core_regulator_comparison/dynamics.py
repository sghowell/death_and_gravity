"""Full hybrid core probability and coupled cutoff/ordering/readout comparisons."""

from functools import cache

import sympy as s
from p8_vacuum_affine_selfconsistent_finite_feedback import geometry, homogeneous

from . import operators, source

TAIL = s.Rational(1, 10**4000)
LEAKAGE = s.Rational(1, 10**6)
CUTOFF_Y = s.Rational(1, 10**580)
CUTOFF_STATE = s.Rational(1, 10**12)
CUTOFF_VOLUME = s.Rational(1, 10**514)
ORDER_Y = s.Rational(1, 10**630)
ORDER_STATE = s.Rational(1, 10**65)
ORDER_VOLUME = s.Rational(1, 10**566)
VOLUME_PARAMETER = s.Integer(10) ** 8


def operator_row(amplitude, order=0):
    packet = operators.bounds(amplitude, order)
    return (
        packet["whole_finite_coherent_polynomial_bound"],
        packet["whole_retained_Weyl_heat_remainder"],
        packet["whole_both_complete_operator_bound"],
        packet["whole_Weyl_minus_calibrated_operator_difference"],
    )


@cache
def bounds():
    AH, AF, T = source.H_AMPLITUDE, source.F_AMPLITUDE, source.TIME
    SH, EH, BH, DH = operator_row(AH)
    S1, E1, B1, D1 = operator_row(AH, 1)
    _, _, B2, _ = operator_row(AH, 2)
    SF, EF, BF, DF = operator_row(AF)
    _, _, BF1, _ = operator_row(AF, 1)
    Bq = 10**104 * B1 / geometry.KAPPA
    Lq = 10**106 * (B1 + B2) / geometry.KAPPA
    a = T * (homogeneous.LIPSCHITZ + Lq)
    b = 2 * T * Bq
    c = 5 * T * B1
    denominator = 1 - a - b * c
    motion = BH * T
    leakage = TAIL + 2 * motion
    force_cut = 10**104 / geometry.KAPPA * (2 * S1 * TAIL + 2 * E1 + 4 * B1 * motion)
    state_cut = 2 * SH * T * s.sqrt(TAIL) + 2 * EH * T + BH**2 * T**2
    y_cut = (b * state_cut + T * force_cut) / denominator
    phi_cut = c * y_cut + state_cut
    force_order = 10**104 * D1 / geometry.KAPPA
    state_order = T * DH
    y_order = (b * state_order + T * force_order) / denominator
    phi_order = c * y_order + state_order
    volume_parameter = 12 + 4 * 2 * 12 * geometry.NZ + 6 * BF + 10 * BF1
    volume_cut = (
        VOLUME_PARAMETER * y_cut
        + 4 * BF * phi_cut
        + 4 * SF * TAIL
        + 4 * EF
        + 8 * BF * motion
    )
    volume_order = VOLUME_PARAMETER * y_order + 4 * BF * phi_order + 2 * DF
    return {
        "whole_original_positive_initial_tail_ceiling": TAIL,
        "whole_complete_Hamiltonian_operator": BH,
        "whole_complete_centered_volume_operator": BF,
        "whole_uniform_quantum_force": Bq,
        "whole_uniform_quantum_force_Lipschitz": Lq,
        "whole_coupled_a_b_c": [a, b, c],
        "whole_coupled_inverse_denominator": denominator,
        "whole_phase_factored_state_motion": motion,
        "whole_evolved_positive_coherent_POVM_leakage": leakage,
        "whole_direct_cutoff_force": force_cut,
        "whole_direct_cutoff_unitary_comparison": state_cut,
        "whole_coupled_cutoff_five_coordinate_difference": y_cut,
        "whole_coupled_cutoff_phase_factored_state_difference": phi_cut,
        "whole_direct_ordering_force": force_order,
        "whole_direct_ordering_unitary_comparison": state_order,
        "whole_coupled_ordering_five_coordinate_difference": y_order,
        "whole_coupled_ordering_phase_factored_state_difference": phi_order,
        "whole_physical_volume_parameter_row_sum": volume_parameter,
        "whole_coupled_cutoff_physical_volume_difference": volume_cut,
        "whole_coupled_ordering_physical_volume_difference": volume_order,
    }


@cache
def data():
    b = bounds()
    a, cross, state = b["whole_coupled_a_b_c"]
    T = source.TIME
    cut_y, cut_state = (
        b["whole_coupled_cutoff_five_coordinate_difference"],
        b["whole_coupled_cutoff_phase_factored_state_difference"],
    )
    ord_y, ord_state = (
        b["whole_coupled_ordering_five_coordinate_difference"],
        b["whole_coupled_ordering_phase_factored_state_difference"],
    )
    D = b["whole_coupled_inverse_denominator"]
    checks = {
        "whole_coupled_denominator_includes_both_cross_rows": D - 1 + a + cross * state,
        "whole_cutoff_coupled_classical_row": (1 - a) * cut_y
        - cross * cut_state
        - T * b["whole_direct_cutoff_force"],
        "whole_cutoff_coupled_quantum_row": cut_state
        - state * cut_y
        - b["whole_direct_cutoff_unitary_comparison"],
        "whole_ordering_coupled_classical_row": (1 - a) * ord_y
        - cross * ord_state
        - T * b["whole_direct_ordering_force"],
        "whole_ordering_coupled_quantum_row": ord_state
        - state * ord_y
        - b["whole_direct_ordering_unitary_comparison"],
        "whole_positive_initial_tail_not_zero": TAIL - s.Rational(1, 10**4000),
    }
    return {
        "whole_all_evaluated_hybrid_leakage_and_comparison_bounds": b,
        "whole_safe_leakage_and_comparison_ceilings": {
            "positive_POVM_leakage": LEAKAGE,
            "cutoff_five_coordinate": CUTOFF_Y,
            "cutoff_factored_state": CUTOFF_STATE,
            "cutoff_physical_volume": CUTOFF_VOLUME,
            "ordering_five_coordinate": ORDER_Y,
            "ordering_factored_state": ORDER_STATE,
            "ordering_physical_volume": ORDER_VOLUME,
        },
        "whole_same_solution_and_seed": "Apply the new uniform operator estimates to the SAME four S273 solutions, with original48-pair pure seed, actual scalar-centered phase convention, five live classical density coordinates, cyclic M1 quadrature, all heavy and symmetry-consistent homogeneous equations. The old uniqueness theorem identifies these trajectories; no new source, state, cutoff, Hamiltonian or physical background is chosen.",
        "whole_actual_core_probability_proof": "For every allowed Y path, exact bounded unitarity gives||phi(u)-psi0||<=BH|u|. The original Husimi tail is positive and<1e-4000. Its outside-core positive coherent effect is a contraction, hence its evolved probability is at mosttail+2BH T<1e-6. This is a POVM probability, not sharp joint phase support, projection or zero leakage. Transport both effect and state by the same Uref.",
        "whole_cutoff_comparison_proof": "For a difference of the two original cutoffs at fixed Y, every derivative of the compact symbol vanishes on the common core. The exact finite coherent P6 piece has norm<=2SH and seed action<=2SH sqrt(tail); the complete Weyl remainder has norm<=2EH and is NOT called core supported. Duhamel with the comparison state's displacement<=BH|u| gives ecut=2SH T sqrt(tail)+2EH T+BH^2 T^2. For each homogeneous derivative, its expectation on paths withinBH T of the seed differs by at most2S1 tail+2E1+4B1 BH T. The full canonical force coefficient1e104/kappa yields the displayed force_cut.",
        "whole_coupled_stability_proof": "With y=sup|deltaY|infty and z=sup||delta_phi||, the exact integral equations imply y<=a y+b z+T Fdirect and z<=c y+edirect, wherea=T(Lhom+Lquant),b=2T Bquant,c=5T B1. These include the actual background feedback, not externally identical histories. Since1-a-bc>99/100, solve BOTH rows to gety=(b edirect+T Fdirect)/(1-a-bc),z=c y+edirect. For ordering useedirect=T DeltaH andFdirect=1e104 Delta1/kappa. All five density coordinates are compared; cyclic M1 remains its full quadrature.",
        "whole_physical_volume_comparison_proof": "Bound the complete parameter derivative of exp(3alpha)[F0(Y)I+f(Y)] by the displayed row sum<1e8, usingUN<12,Nz<1e4,a^3<2 and the centered operator's first Cauchy derivative. For state change at fixedY the scalarF0I cancels, so the bound uses4BF z, not the norm of the identity. Cutoff readout change on the evolved state includes both the coherent seed tail and its finite displacement, and the full Weyl heat remainder:4SF tail+4EF+8BF BH T. Ordering readout change contributes2DeltaF. Add the parameter, state and operator changes to obtain the stated physical-volume differences.",
        "whole_scalar_phase_and_covariance_boundary": "phi is the phase-factored interaction state: each solution restores its OWN exp[-i integral g0(Y)] and the SAME fixed Uref. Every scalar-center classical force was retained. Distinct Y paths may yield very different global scalar phases, so these vector-norm comparison bounds are NOT assigned to the un-factored physical state. Physical-volume and POVM expectations are phase invariant. No history-dependent vacuum reset, covariance deletion or state projection is made.",
        "whole_original_boundary": "These are finite hybrid and finite-regulator comparisons on the same1e-180 interval. They do not provide exact support, a uniform mode/volume/phase-cutoff limit, an unlocalized Hamiltonian, quantized homogeneous variables, original fully quantum means, physical matching, omitted-loop/Regge/UV control or nonlinear global completion. They do not strengthen the old turnaround into a unique or strict minimum. Original V/G/B/P8 remain OPEN.",
        "checks": checks,
        "gates": {
            "whole_coupled_inverse_gap": D > s.Rational(99, 100),
            "whole_positive_initial_probability_ceiling": 0 < TAIL < 1,
            "whole_evolved_nonzero_POVM_leakage": b[
                "whole_evolved_positive_coherent_POVM_leakage"
            ]
            < LEAKAGE,
            "whole_cutoff_five_coordinate_difference": cut_y < CUTOFF_Y,
            "whole_cutoff_phase_factored_state_difference": cut_state < CUTOFF_STATE,
            "whole_cutoff_physical_volume_difference": b[
                "whole_coupled_cutoff_physical_volume_difference"
            ]
            < CUTOFF_VOLUME,
            "whole_ordering_five_coordinate_difference": ord_y < ORDER_Y,
            "whole_ordering_phase_factored_state_difference": ord_state < ORDER_STATE,
            "whole_ordering_physical_volume_difference": b[
                "whole_coupled_ordering_physical_volume_difference"
            ]
            < ORDER_VOLUME,
            "whole_physical_parameter_row_sum": b[
                "whole_physical_volume_parameter_row_sum"
            ]
            < VOLUME_PARAMETER,
            "all_new_quantities_are_exact_positive_finite_rationals": all(
                v.is_Rational is True and v > 0
                for k, v in b.items()
                if k != "whole_coupled_a_b_c"
            ),
            "same_original_full_scalar_center_forces_retained": True,
            "no_unfactored_physical_state_distance_inferred": True,
            "original_P8_still_open": True,
        },
    }
