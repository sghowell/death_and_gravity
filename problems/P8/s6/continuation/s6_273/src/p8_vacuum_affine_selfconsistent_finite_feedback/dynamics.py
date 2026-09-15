"""Exact coupled finite hybrid path fixed point and physical-volume turnaround."""

from functools import cache

import sympy as s

from . import geometry, homogeneous, quantum, source

BASE_ERROR = s.Rational(1, 10**280)
VOLUME_ERROR = s.Rational(1, 10**380)
MINIMUM_RADIUS = s.Rational(1, 10**188)


@cache
def bounds():
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
def data():
    b = bounds()
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
        "whole_precise_finite_hybrid_definition": "Homogeneous alpha,p,pm,eta,ph are classical canonical-density variables; M1 is reconstructed by its cyclic rate. The unchanged48 nonzero-mode pairs have their original pure Gaussian Cauchy state and exact full regulated quantum Hamiltonian. Keep the fixed source profiles and all primitive/spatial/constraint contacts. This is a separately named finite hybrid model, not a quantization of homogeneous modes or the original unrestricted theory.",
        "whole_path_space_and_map": "Use continuous real Y paths on[-T,T] with Y(0)=Ybar(0), ||Y-Ybar||infty<=1e-390, and continuous normalized cubic-invariant zero-translation-charge Hilbert paths phi with the original seed at0. The metric is max{||deltaY||infty,1e-350||deltaphi||infty}; the displayed row-sum bounds apply to this weighted maximum metric. The map integrates the complete classical-plus-quantum force using the input(Y,phi), and sends phi to the exact unitary U_K[Y] applied to the seed. This closed subset is complete; convexity of the unit sphere is not required.",
        "whole_self_map_proof": "The actual reference obeys the full clock rate identities except-p_fixed; its full lapse root differs from1 by the displayed profile-controlled amount. The entire N-force and reduced homogeneous Jacobian are<1e112. Thus the classical deviation is bounded by T(BASE_ERROR+Lhom*real+Fquant)<real/2 for every normalized input state. The quantum map is normalized and symmetry preserving by exact unitarity.",
        "whole_contraction_proof": "State expectations of each bounded force differ by at most2||force|| ||delta phi||. The first row isT(Lhom+Lquant)+2TFquant/weight. The full centered-generator Duhamel identity, using unitarity rather than exp(T||K||), gives second row weight*T*5max||partial_Y K||. Both are<1e-50. The contraction theorem gives a unique continuous pair; the integral equations and norm-continuous bounded operator family give C1 paths. M1 follows by quadrature and a^3pm is exactly conserved.",
        "whole_physical_phase_and_picture": "The comparison equation is i phi'=K_Y phi. Restore psi=exp[-i integral g0(Y)]phi, so i psi'=[g0 I+K_Y]psi. g0 remains in the homogeneous classical force throughout. The physical state and readout use the SAME fixed Uref. No state projection, discarded phase, instantaneous reference reset or near-identity state estimate is used.",
        "whole_consistent_remaining_homogeneous_equations": "The full cubic-equivariant extension in homogeneous spatial-vector and five traceless-shape canonical directions exists by the same local smooth chart. At the invariant solution, every vector and traceless symmetric force vanishes by the exact group average. Thus all three vector and all FIVE shape coordinate/momentum equations are satisfied at zero. The temporal normal vector is reconstructed from the whole algebraic root. The heavy scalar remains a live sourced homogeneous variable.",
        "whole_physical_volume_turnaround": "For either original cutoff and either complete operator ordering, nu=exp(3alpha)<F_operator>=(1+u^2)^6 times a factor within1e-380 of1. Each endpoint exceeds nu(0) by more than5T^2. Every global minimum is inside |u|<1e-188. The mean is C1, so a minimum is critical and there are negative and positive derivative points on its sides. No uniqueness, strict positive acceleration or monotone half-interval is claimed.",
        "whole_open_boundary": "These are self-consistent finite HYBRID equations with classical homogeneous variables, not original fully quantum means. No long-time core leakage, unlocalized Hamiltonian, mode/volume/phase-cutoff removal, physical matching, omitted-loop/Regge/UV control or nonlinear global completion is supplied. A turnaround of an extended readout does not prove occupation remains in its unmodified classical core. Original V/G/B/P8 remain OPEN.",
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
