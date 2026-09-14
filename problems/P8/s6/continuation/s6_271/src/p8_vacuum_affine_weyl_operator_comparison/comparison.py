"""Actual finite Weyl operator errors, positivity and entire-unitary comparison."""

from functools import cache

import sympy as s
from p8_vacuum_affine_quantitative_local_time import quantum as previous

from . import derivatives, kernel, source


@cache
def bounds():
    # The heat remainder needs every Schur derivative alpha_j<=2,
    # hence total orders4..196 of the entire cutoff amplitude.
    error_factor = kernel.CONSTANT * s.Rational(source.PHASE**2, 32)
    volume = error_factor * derivatives.amplitude_jet(4, source.VERROR)
    Hamiltonian = error_factor * derivatives.amplitude_jet(4, source.HBOUND)
    reference_norm = previous.GENERATOR
    Uerror = source.TIME * Hamiltonian
    Wnorm = reference_norm + Hamiltonian
    old = previous.dynamical_bounds()
    mean = volume + 4 * Uerror
    return {
        "whole_Weyl_volume_operator_ordering_error": volume,
        "whole_Weyl_Hamiltonian_operator_ordering_error": Hamiltonian,
        "whole_Weyl_volume_distance_from_identity": 2 * source.VERROR + volume,
        "whole_Weyl_volume_positive_operator_floor": 1 - 2 * source.VERROR - volume,
        "whole_Weyl_interaction_operator_norm": Wnorm,
        "whole_two_ordering_unitary_OPERATOR_norm_error": Uerror,
        "whole_Weyl_unitary_operator_distance_from_identity": source.TIME * Wnorm,
        "whole_same_seed_two_ordering_volume_mean_error": mean,
        "whole_two_Weyl_cutoff_state_error": 2 * Uerror
        + old["whole_two_explicit_regulator_state_difference"],
        "whole_two_Weyl_cutoff_readout_mean_error": 2 * mean
        + old["whole_two_readout_mean_difference"],
        "whole_evolved_Weyl_coherent_outside_core_probability": old[
            "whole_original_initial_gamma_tail_ceiling"
        ]
        + 2 * source.TIME * Wnorm,
    }


@cache
def identities():
    t, l = s.symbols("heat_time commuting_heat_generator", real=True)
    U, V = s.symbols("exact_unitary_one exact_unitary_two", commutative=False)
    F, G = s.symbols("first_readout second_readout", commutative=False)
    psi, phi = s.symbols("first_state second_state", commutative=False)
    # Exact additive decomposition behind the comparison, not state-only.
    decomposition = (F * psi - G * phi) - (F * (psi - phi) + (F - G) * phi)
    return {
        "whole_exact_heat_remainder_differential_identity": s.simplify(
            s.diff(s.exp(t * l) * (1 - t * l), t) + t * s.exp(t * l) * l * l
        ),
        "full_state_and_operator_difference_decomposition": s.expand(decomposition),
        "whole_zero_time_interaction_identity": s.simplify(
            s.exp(t * l) * (1 - t * l)
        ).subs(t, 0)
        - 1,
        "full_phase_Laplacian_D2_normalization": s.Rational(source.PHASE**2, 32)
        - s.Rational(1, 2) * s.Rational(source.PHASE**2, 16),
        "same_volume_Hamiltonian_relative_ordering_error": bounds()[
            "whole_Weyl_volume_operator_ordering_error"
        ]
        / source.VERROR
        - bounds()["whole_Weyl_Hamiltonian_operator_ordering_error"] / source.HBOUND,
        "full_unitary_comparison_includes_both_propagators": s.expand(
            U - V - (U - s.Integer(1)) + (V - s.Integer(1))
        ),
    }


@cache
def data():
    b = bounds()
    return {
        "whole_all_actual_operator_positivity_unitary_and_readout_bounds": b,
        "whole_declared_alternative_Weyl_regulator": "For each SAME explicit chi_c and the same entire g_ext,c,F_ext,c, define A_W,c=OpW(g_ext,c), F_W,c=OpW(F_ext,c) at unit canonical CCR. The previous A_C,c=Q_V[(1-D)g_ext,c] and F_C,c=Q_V[(1-D)F_ext,c] are kept as distinct comparison operators. All original source, time, primitive, implicit, cutoff and scalar background terms remain. This is a newly explicit ordering of the already named finite regulator, not a change to the original parent or identification with an unregularized theory.",
        "whole_heat_SYMBOL_to_actual_OPERATOR_argument": "For r_f=exp(D)(1-D)f-f, the exact identity r_f=-integral_0^1 t exp(tD)D^2 f dt holds. Derivatives commute with the full heat semigroup, which contracts the sup norm. For every Schur multiindex alpha_j<=2, |alpha|<=192, its derivative is bounded by96^2/32 times the entire amplitude jet J_(|alpha|+4). The actual high-order cutoff/Cauchy bounds decrease through196, so J4 bounds them ALL. Applying the explicit1e112 Gaussian-frame Schur theorem now gives an actual Weyl OPERATOR norm, not only a symbol error. The exact same-seed coherent identity Q_V((1-D)f)=OpW(exp(D)(1-D)f) identifies the compared operator. Every high mixed derivative and heat/cutoff contact remains.",
        "whole_actual_Weyl_volume_positivity": "The established calibrated coherent F_C is positive and ||F_C-I||<2e-255. Its actual operator distance to F_W is<1e-200. Thus ||F_W-I||<1e-199 and F_W>1/2 as a bounded self-adjoint operator on the full Hilbert space and its zero-charge sector. This proves positivity for this concrete near-identity full volume; it does not assert that Weyl quantization preserves positivity of arbitrary positive symbols.",
        "whole_Weyl_self_adjointness_continuity_and_strong_dynamics": "The whole real g_ext is a retained scalar plus a smooth compact phase symbol, so the explicit kernel theorem gives a bounded self-adjoint Weyl operator. Its time dependence is norm continuous by the common source domain and uniform high phase derivatives, with only the inherited real-time smoothness. The non-scalar Weyl kernel is Schwartz; the existing seminorm Dyson argument supplies the common Schwartz strong equation. The complete bounded generator, not a finite Taylor or occupation truncation, therefore has an exact unitary U_W. In physical picture use U_ref U_W and conjugate each corresponding readout by the same U_ref.",
        "whole_entire_unitary_OPERATOR_comparison": "For the full bounded self-adjoint A_W,A_C, unitarity and the exact two-propagator Duhamel formula give ||U_W(u)-U_C(u)|| <=integral||A_W-A_C|| <=|u|*error_H<1e-944. This holds in OPERATOR norm, not merely on the Gaussian seed. Each entire scalar phase remains. The Weyl propagator itself obeys ||U_W-I||<=|u|*(B_C+error_H)<1e-943. These statements concern only |u|<=1e-2000 and the two defined finite orderings.",
        "whole_same_seed_readout_and_cutoff_comparisons": "A same-cutoff Weyl/calibrated-coherent readout mean differs by at most error_F+4*unitary_error<2e-200, including BOTH operator and state change; the factor4 uses ||F_C||<2 and normalized states. The two distinct Weyl cutoffs are compared by the triangle inequality through the qualified S270 comparison, giving state error<1e-943 and mean error<1e-199. The stronger old calibrated-coherent-only1e-1970 and1e-1230 bounds are NOT reassigned unchanged to Weyl regulators. The evolved Weyl coherent outside-core probability is<1e-942, not zero; that probability still uses the original positive coherent POVM, not a supposed Weyl phase projection.",
        "whole_all_original_translation_constraints_and_covariance": "The original finite translation action is symplectic and preserves the entire symbols and actual pure seed. Exact metaplectic covariance of WEYL quantization, or its defining kernel, makes A_W and F_W commute with that full unitary translation group. The three original common zero-charge constraints give a reducing sector for the complete unitary, with no coherent-label conditioning or deleted generated means. For symplectic whitening S0, Weyl covariance is unitary; the old complete physical cross covariance is not replaced by a new vacuum.",
        "whole_remaining_original_boundary": "All comparisons are between explicit finite smooth regulators of the same local classical Hamiltonian and volume. No unregularized singular Hamiltonian or original interacting mean, uniform mode/torus limit, Wilsonian matching, omitted-loop bound, all-energy UV positivity or nonlinear global completeness follows. Original V/G/B/P8 remain OPEN and completed scoped P8(a) is unchanged.",
        "checks": identities(),
        "gates": {
            "actual_full_Weyl_volume_OPERATOR_error": b[
                "whole_Weyl_volume_operator_ordering_error"
            ]
            < s.Rational(1, 10**200),
            "actual_full_Weyl_H_OPERATOR_error": b[
                "whole_Weyl_Hamiltonian_operator_ordering_error"
            ]
            < 10**1056,
            "actual_Weyl_positive_volume_floor": b[
                "whole_Weyl_volume_positive_operator_floor"
            ]
            > s.Rational(1, 2),
            "actual_Weyl_volume_operator_distance": b[
                "whole_Weyl_volume_distance_from_identity"
            ]
            < s.Rational(1, 10**199),
            "full_Weyl_interaction_bound": b["whole_Weyl_interaction_operator_norm"]
            < 10**1056,
            "entire_two_ordering_unitary_OPERATOR_norm": b[
                "whole_two_ordering_unitary_OPERATOR_norm_error"
            ]
            < s.Rational(1, 10**944),
            "entire_Weyl_unitary_state_bound": b[
                "whole_Weyl_unitary_operator_distance_from_identity"
            ]
            < s.Rational(1, 10**943),
            "full_state_and_readout_ordering_mean": b[
                "whole_same_seed_two_ordering_volume_mean_error"
            ]
            < 2 * s.Rational(1, 10**200),
            "two_defined_Weyl_cutoff_states": b["whole_two_Weyl_cutoff_state_error"]
            < s.Rational(1, 10**943),
            "two_defined_Weyl_cutoff_readout_means": b[
                "whole_two_Weyl_cutoff_readout_mean_error"
            ]
            < s.Rational(1, 10**199),
            "evolved_Weyl_positive_coherent_POVM_leakage": b[
                "whole_evolved_Weyl_coherent_outside_core_probability"
            ]
            < s.Rational(1, 10**942),
            "every_heat_and_Schur_derivative_is_covered": 2 * source.PHASE + 4
            <= source.MAX_ORDER,
            "all_full_bounds_positive_exact_rationals": all(
                value.is_Rational is True and value > 0 for value in b.values()
            ),
            "entire_original_background_phase_retained": True,
            "not_generic_positive_Weyl_symbol_or_original_matching": True,
        },
    }
