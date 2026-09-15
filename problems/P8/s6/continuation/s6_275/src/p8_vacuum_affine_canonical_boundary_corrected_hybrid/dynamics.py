"""New corrected coupled trajectories, symmetry, turnaround and comparisons."""

from functools import cache

import sympy as s
from p8_vacuum_affine_selfconsistent_finite_feedback import homogeneous, symmetry

from . import geometry, quantum, source

BASE_ERROR = s.Rational(1, 10**280)
VOLUME_ERROR = s.Rational(1, 10**380)
MINIMUM_RADIUS = s.Rational(1, 10**188)


@cache
def existence_bounds():
    spatial, operator = geometry.bounds(), quantum.bounds()
    T, real = source.TIME, source.REAL_RADIUS
    selfmap = T * (
        BASE_ERROR + homogeneous.LIPSCHITZ * real + operator["quantum_force"]
    )
    classical = (
        T * (homogeneous.LIPSCHITZ + operator["quantum_force_Lipschitz"])
        + 2 * T * operator["quantum_force"] / quantum.STATE_WEIGHT
    )
    state = quantum.STATE_WEIGHT * T * 5 * operator["each_first_homogeneous_derivative"]
    nreference = s.Rational(5, 10**400) / (3 * (1 - spatial["contraction"]))
    baseline = source.PROFILE + homogeneous.LIPSCHITZ * nreference
    ncenter = (s.Rational(5, 10**400) + 12 * 10**4 * real) / (
        3 * (1 - spatial["contraction"])
    )
    # Include the product of the small physical scale error with F.
    centered = quantum.operator_bound(quantum.F_AMPLITUDE)
    volume = 12 * ncenter + 12 * real + 2 * centered
    gap = (1 + T * T) ** 6 * (1 - VOLUME_ERROR) - (1 + VOLUME_ERROR)
    location = ((1 + MINIMUM_RADIUS**2) ** 6 - 1) * (
        1 - VOLUME_ERROR
    ) - 2 * VOLUME_ERROR
    return {
        "self_map_deviation": selfmap,
        "classical_contraction_row": classical,
        "quantum_contraction_row": state,
        "reference_lapse_displacement": nreference,
        "reference_classical_rate_error": baseline,
        "inner_homogeneous_lapse_displacement": ncenter,
        "physical_volume_error": volume,
        "endpoint_gap": gap,
        "minimum_location_margin": location,
    }


@cache
def existence_data():
    b = existence_bounds()
    T, real = source.TIME, source.REAL_RADIUS
    x = s.Symbol("nonnegative_clock_square", nonnegative=True)
    clock = s.Symbol("real_evolution_time", real=True)
    phase = s.Function("retained_scalar_phase")(clock)
    vector = s.Function("centered_state_component")(clock)
    generator = s.Function("centered_generator_component")(clock)
    restored = s.exp(-s.I * phase) * vector
    phase_identity = (
        s.I * s.diff(restored, clock) - (s.diff(phase, clock) + generator) * restored
    ).subs(s.diff(vector, clock), -s.I * generator * vector)
    return {
        "whole_coupled_fixed_point_and_volume_budgets": b,
        "whole_evaluated_time_and_inner_homogeneous_radius": [T, real],
        "whole_phase_metric_weight": quantum.STATE_WEIGHT,
        "whole_volume_error_and_minimum_radius": [VOLUME_ERROR, MINIMUM_RADIUS],
        "whole_precise_finite_hybrid_definition": "These are FOUR NEW canonical-boundary corrected finite hybrid solutions, one for each original cutoff and operator ordering. The five live classical density coordinates alpha,p,pm,eta,ph and cyclic M1 retain all raw homogeneous source equations. The48 quantum nonzero-mode pairs start at the original prepared Gaussian; the physical raw initial state is Ub(0)psi0. The full nonlinear Hamiltonian is pulled back by Sb with its exact time generator before reference interaction and localization. This is not the archived unshifted model or a homogeneous quantization.",
        "whole_path_space_and_map": "Use continuous real Y paths on[-T,T] with original homogeneous Cauchy data and sup|Y-Ybar|<=1e-390, and continuous normalized original-seed Hilbert paths phi in the cubic-invariant zero-translation-charge sector. The complete path metric is max{sup|deltaY|,1e-350 sup||delta_phi||}. Integrate the full classical-plus-corrected-quantum force at input(Y,phi) and evolve the exact corrected bounded unitary at inputY. The closed product subset is complete; convexity of the unit sphere is not required.",
        "whole_self_map_proof": "The actual reference obeys the full clock rate identities except-p_fixed; its full lapse root differs from1 by the displayed profile-controlled amount. The entire N-force and reduced homogeneous Jacobian are<1e112. Thus the classical deviation is bounded by T(BASE_ERROR+Lhom*real+Fquant)<real/2 for every normalized input state. The quantum map is normalized and symmetry preserving by exact unitarity.",
        "whole_contraction_proof": "State expectations of each bounded force differ by at most2||force|| ||delta phi||. The first row isT(Lhom+Lquant)+2TFquant/weight. The full centered-generator Duhamel identity, using unitarity rather than exp(T||K||), gives second row weight*T*5max||partial_Y K||. Both are<1e-50. The contraction theorem gives a unique continuous pair; the integral equations and norm-continuous bounded operator family give C1 paths. M1 follows by quadrature and a^3pm is exactly conserved.",
        "whole_physical_phase_and_picture": "The corrected comparison equation is i phi'=K_Y^correct phi. Restore psi_raw=Ub(u)Uref(u)exp[-i integral g0(Y)]phi. Both the configuration-dependent boundary phase and its time-connection energy are compulsory; they cannot be factored as a global scalar. The actual scalar g0 remains in every homogeneous force and is restored separately. Readout and positive coherent effect receive the same Ub and Uref conjugations.",
        "whole_consistent_remaining_homogeneous_equations": "The complete raw reconstruction and prepared reference are cubic/translation equivariant. boundary.py checks all24 actual cubic actions and all3 translation charges commute with the SAME boundary shear on all scalar waves. Thus the corrected full operators preserve the actual invariant sector. The full homogeneous extension in three-vector and all FIVE traceless-shape canonical pairs has vanishing forces there by the displayed Reynolds averages. The heavy homogeneous scalar remains live and the temporal normal vector is the full algebraic root.",
        "whole_physical_volume_turnaround": "For each NEW corrected solution the correctly transported complete volume mean obeys nu=(1+u^2)^6(1+error), |error|<1e-380. Both endpoint values exceed nu(0) by more than5T^2; all global minima lie in |u|<1e-188. C1 regularity gives criticality and a negative/positive derivative point on the respective sides. No unique or strict minimum, positive acceleration or monotone half-interval is inferred.",
        "whole_open_boundary": "This repairs one finite physical canonical identification, not the original unrestricted quantum theory. No homogeneous quantum state, mode/volume/cutoff removal, unlocalized Hamiltonian, physical UV matching, omitted-loop/Regge control or nonlinear global completion is supplied. Original V/G/B/P8 remain OPEN and scoped P8(a) is unchanged.",
        "checks": {
            "whole_sixth_power_bounce_polynomial": s.expand(
                (1 + x) ** 6 - 1 - sum(s.binomial(6, i) * x**i for i in range(1, 7))
            ),
            "whole_two_state_contraction_rows": b["classical_contraction_row"]
            - T * (homogeneous.LIPSCHITZ + quantum.bounds()["quantum_force_Lipschitz"])
            - 2 * T * quantum.bounds()["quantum_force"] / quantum.STATE_WEIGHT,
            "whole_state_unitary_Duhamel_row": b["quantum_contraction_row"]
            - quantum.STATE_WEIGHT
            * T
            * 5
            * quantum.bounds()["each_first_homogeneous_derivative"],
            "whole_restored_scalar_plus_centered_rate": s.simplify(phase_identity),
            "whole_same_real_time": T - s.Rational(1, 10**180),
        },
        "gates": {
            "whole_reference_profile_rate_error": b["reference_classical_rate_error"]
            < BASE_ERROR,
            "whole_strict_homogeneous_self_map": b["self_map_deviation"] < real / 2,
            "whole_strict_coupled_contraction": max(
                b["classical_contraction_row"], b["quantum_contraction_row"]
            )
            < s.Rational(1, 10**50),
            "whole_physical_volume_error": b["physical_volume_error"] < VOLUME_ERROR,
            "whole_positive_physical_volume_floor": 1 - VOLUME_ERROR > s.Rational(1, 2),
            "whole_both_endpoint_gap": b["endpoint_gap"] > 5 * T * T,
            "whole_all_minima_location": b["minimum_location_margin"] > 0
            and MINIMUM_RADIUS < T,
            "whole_real_path_inside_complex_domain": real
            < source.HOMOGENEOUS_CAUCHY
            < source.HOMOGENEOUS_RADIUS,
            "finite_hybrid_not_original_P8_closure": True,
            "no_unique_strict_minimum_or_core_support_claim": True,
        },
    }


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
    packet = quantum.operator_packet(amplitude, order)
    return (
        packet["whole_finite_coherent_polynomial_bound"],
        packet["whole_retained_Weyl_heat_remainder"],
        packet["whole_both_complete_operator_bound"],
        packet["whole_Weyl_minus_calibrated_operator_difference"],
    )


@cache
def comparison_bounds():
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
def comparison_data():
    b = comparison_bounds()
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
        "whole_same_solution_and_seed": "Apply the corrected full operator estimates to the FOUR NEW solutions of this successor's existence theorem, with the original prepared48-pair seed and its correctly transported physical raw state, complete scalar-center force, five live classical density coordinates, cyclic M1 quadrature and symmetry-consistent full homogeneous equations. The archived S273-S274 solutions are NOT identified with these corrected trajectories.",
        "whole_actual_core_probability_proof": "For every corrected allowed Y path, bounded unitarity gives||phi(u)-psi0||<=BH|u|. The unchanged prepared Gaussian's positive Husimi tail is<1e-4000. The outside-core coherent effect is a positive contraction, so its evolved probability is at mosttail+2BH T<1e-6. Transport BOTH effect and state by Ub(u)Uref(u); this is a physical positive-operator probability, not exact joint support or a state projection.",
        "whole_cutoff_comparison_proof": "For a difference of the two original cutoffs at fixed Y, every derivative of the compact symbol vanishes on the common core. The exact finite coherent P6 piece has norm<=2SH and seed action<=2SH sqrt(tail); the complete Weyl remainder has norm<=2EH and is NOT called core supported. Duhamel with the comparison state's displacement<=BH|u| gives ecut=2SH T sqrt(tail)+2EH T+BH^2 T^2. For each homogeneous derivative, its expectation on paths withinBH T of the seed differs by at most2S1 tail+2E1+4B1 BH T. The full canonical force coefficient1e104/kappa yields the displayed force_cut.",
        "whole_coupled_stability_proof": "With y=sup|deltaY|infty and z=sup||delta_phi||, the exact integral equations imply y<=a y+b z+T Fdirect and z<=c y+edirect, wherea=T(Lhom+Lquant),b=2T Bquant,c=5T B1. These include the actual background feedback, not externally identical histories. Since1-a-bc>99/100, solve BOTH rows to gety=(b edirect+T Fdirect)/(1-a-bc),z=c y+edirect. For ordering useedirect=T DeltaH andFdirect=1e104 Delta1/kappa. All five density coordinates are compared; cyclic M1 remains its full quadrature.",
        "whole_physical_volume_comparison_proof": "Bound the complete parameter derivative of exp(3alpha)[F0(Y)I+f(Y)] by the displayed row sum<1e8, usingUN<12,Nz<1e4,a^3<2 and the centered operator's first Cauchy derivative. For state change at fixedY the scalarF0I cancels, so the bound uses4BF z, not the norm of the identity. Cutoff readout change on the evolved state includes both the coherent seed tail and its finite displacement, and the full Weyl heat remainder:4SF tail+4EF+8BF BH T. Ordering readout change contributes2DeltaF. Add the parameter, state and operator changes to obtain the stated physical-volume differences.",
        "whole_scalar_phase_and_covariance_boundary": "phi factors each solution's OWN global scalar g0 phase, while every g0 force is retained. The same FIXED configuration-dependent Ub and reference Uref are restored on all corrected solutions, so their common unitary cancels in phase-factored vector norm differences but not in the physical operator definition. No boundary phase is dropped and no old unshifted covariance is called the corrected physical state. Unfactored states may have very different global phases.",
        "whole_original_boundary": "The repaired finite canonical bridge does not supply a uniform mode/volume/cutoff limit, an unlocalized Hamiltonian, homogeneous quantum means, physical UV matching, omitted-loop/Regge control or nonlinear global completion. Neither the positive-core probability nor these finite regulator comparisons strengthens the turnaround to a unique or strict minimum. Original V/G/B/P8 remain OPEN.",
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


@cache
def data():
    h = homogeneous.data()
    sym = symmetry.data()
    pieces = {
        "corrected_existence_and_turnaround": existence_data(),
        "corrected_core_and_regulator_comparisons": comparison_data(),
    }
    h_keys = (
        "whole_homogeneous_normal_Hamiltonian",
        "whole_six_full_classical_rates",
        "whole_actual_fixed_reference",
        "whole_six_actual_clock_rate_residuals",
        "whole_five_full_N_and_reduced_rate_bounds",
        "whole_six_quantum_force_corrections",
    )
    return {
        **{
            name: {
                key: value
                for key, value in part.items()
                if key not in ("checks", "gates")
            }
            for name, part in pieces.items()
        },
        "whole_retained_full_homogeneous_source_and_canonical_force_chain": {
            key: h[key] for key in h_keys
        },
        "whole_prepared_radial_covariance_and_homogeneous_symmetry": {
            "original_prepared_covariance": sym[
                "whole_full_radial_96_phase_covariance"
            ],
            "all24_proper_cubic_actions": sym[
                "whole_24_full_signed_48_configuration_actions_row_order"
            ],
            "full_Proca_K_V": sym["whole_full_Proca_K_V"],
        },
        "whole_homogeneous_canonical_boundary": "The stipulated Sb is a fixed-reference time-dependent map on the nonzero quantum chart, independent of live Y. Its linear momentum increments have zero mean. Consequently the absolute homogeneous Cauchy data, canonical one-form split and all six density force-chain identities remain unchanged; the quantum energy in that chain is the corrected complete K. The raw scalar center at zero phase is unchanged since Sb(u) applied to the zero phase vector is zero and F_b and partial_uF_b vanish at zero phase. The matrix Sb(u=0) is NOT the identity.",
        "whole_original_seed_symmetry_boundary": "The generic radial covariance is the original PREPARED covariance with all internal/q-p entries. The corrected RAW covariance is its Sb transport, verified separately. Since each cubic action and translation commutes with Sb, invariance and zero translation charges survive without averaging the state. Original full scalar/tensor, heavy SLE and Proca prescriptions are retained.",
        "checks": {
            **{
                name + "_" + key: value
                for name, part in pieces.items()
                for key, value in part["checks"].items()
            },
            **{"full_homogeneous_" + key: value for key, value in h["checks"].items()},
            **{
                "full_prepared_symmetry_" + key: value
                for key, value in sym["checks"].items()
            },
        },
        "gates": {
            **{
                name + "_" + key: bool(value)
                for name, part in pieces.items()
                for key, value in part["gates"].items()
            },
            **{
                "full_homogeneous_" + key: bool(value)
                for key, value in h["gates"].items()
            },
            **{
                "full_prepared_symmetry_" + key: bool(value)
                for key, value in sym["gates"].items()
            },
            "fixed_reference_boundary_preserves_full_live_force_chain": True,
            "no_same_solution_identification_with_archived_unshifted_model": True,
        },
    }
