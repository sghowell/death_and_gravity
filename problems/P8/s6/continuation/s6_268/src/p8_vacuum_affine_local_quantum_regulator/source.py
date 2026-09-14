"""Whole original source, time-dependent canonical chart and free reference."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import phase
from p8_vacuum_affine_nonlinear_auxiliary_measure import canonical
from p8_vacuum_affine_nonlinear_lapse_branch import response
from p8_vacuum_affine_nonlinear_lapse_branch import source as auxiliary
from p8_vacuum_affine_nonlinear_spatial_reduction import translations


@cache
def moving_chart():
    u, q, p, dq = s.symbols(
        "clock log_coordinate canonical_momentum log_velocity", real=True
    )
    a = s.Function("positive_scale", positive=True)(u)
    b = s.Function("moving_momentum_mean", real=True)(u)
    old_q = a * s.exp(q)
    old_p = (p + b) / old_q
    pullback = s.expand(old_p * (s.diff(old_q, q) * dq + s.diff(old_q, u)))
    beta = b * q
    rate = s.diff(a, u) / a
    connection = -(p + b) * rate + s.diff(b, u) * q
    H = s.Symbol("entire_old_Hamiltonian", real=True)
    v, Pv, Hb = s.symbols("hat_log_scale Pi_v hat_Hubble", real=True)
    return {
        "whole_general_chart_rule": "If j_u^*Theta at fixed u equals P dQ+d_z beta_u, then h_red=H_parent(j_u)-Theta_parent(partial_u j_u)+partial_u beta_u. The chosen S267 cotangent chart and constant physical kappa normalization have beta=0. Any subsequent canonical boundary shear requires its entire beta; the original primitive boundary is already in the parent.",
        "whole_nonlinear_point_chart": [old_q, old_p],
        "whole_fixture_exact_boundary": beta,
        "whole_fixture_time_connection": connection,
        "whole_fixture_reduced_Hamiltonian": H + connection,
        "whole_metric_scale_connection": -Hb * Pv,
        "checks": {
            "nonlinear_point_chart_full_symplectic_pairing": s.diff(old_q, q)
            * s.diff(old_p, p)
            - 1,
            "whole_moving_chart_boundary_and_time_contact": s.expand(
                pullback - p * dq - s.diff(beta, q) * dq - s.diff(beta, u) + connection
            ),
            "deleting_boundary_time_derivative_is_nonzero": s.diff(beta, u)
            - q * s.diff(b, u),
            "whole_metric_scale_connection_from_Pi_v": -Hb
            * (2 * Pv * s.exp(-2 * v) / 2)
            * s.exp(2 * v)
            + Hb * Pv,
        },
        "gates": {
            "moving_mean_time_contact_not_deleted": connection.has(s.diff(b, u)),
            "scale_time_contact_not_deleted": connection.has(rate),
            "no_extra_boundary_phase_silently_removed": True,
        },
    }


@cache
def data():
    full = canonical.full_trace()
    bindings = canonical.current_data()
    ref = translations.data()
    branch = auxiliary.data()
    contact = response.physical()
    moving = moving_chart()
    checks = dict(moving["checks"])
    checks["literal_original_bounce_slice_of_whole_time_parent"] = branch["checks"][
        "literal_entire_parent_Hamiltonian_slice"
    ]
    checks["original_complete_temporal_elimination"] = full["checks"][
        "entire_reduced_Hamiltonian"
    ]
    checks["original_heavy_mass_binding"] = bindings["checks"][
        "actual_fixed_heavy_mass"
    ]
    checks["physical_kappa_and_density_CCR"] = (
        phase.C * phase.OMEGA * phase.C.T - phase.kappa * phase.a**3 * phase.OMEGA
    )
    checks["complete_physical_volume_second_chain"] = contact["checks"][
        "whole_complete_physical_second_chain"
    ]
    return {
        "whole_original_time_dependent_normal_Hamiltonian": full[
            "whole_reduced_Hamiltonian"
        ],
        "whole_original_temporal_solution": full["whole_temporal_solution"],
        "whole_actual_source_function_bindings": bindings[
            "whole_abstract_to_current_function_bindings"
        ],
        "whole_original_R_F_and_fixed_profiles": bindings[
            "whole_actual_R_F_fixed_profile_and_constant_bindings"
        ],
        "whole_fixed_profile_bindings": bindings["whole_fixed_profile_bindings"],
        "whole_primitive_boundary": bindings["primitive_boundary"],
        "whole_original_spatial_generator": full[
            "whole_spatial_shift_generator_on_primary_surface"
        ],
        "whole_original_reference_maps": ref["whole_original_scalar_canonical_maps"],
        "whole_original_reference_prescription": "Evaluate the SAME S251 full scalar/tensor, original S240 H SLE and S190 three-mode Proca prescriptions at finitely many nonzero torus modes, propagated to the fixed Cauchy slice u=0. Keep both complete canonical shears, central map, physical kappa and Fourier normalization. Do not replace H by the W6 comparison state or drop internal cross covariances. No homogeneous quantum state or literal R3 identification is made.",
        "whole_physical_unit_CCR_map": phase.C,
        "whole_all_implicit_lapse_first_second_contacts": {
            "first": response.jets()["whole_all_twelve_first_implicit_derivatives"],
            "second": response.jets()[
                "whole_all_144_mixed_second_implicit_derivatives"
            ],
            "canonical_pullback": response.jets()["whole_canonical_pullback_rule"],
        },
        "whole_physical_volume_first_second_variations": [
            contact["whole_complete_physical_volume_first_variation"],
            contact["whole_complete_physical_volume_second_variation"],
        ],
        "whole_canonical_time_connection": {
            k: v for k, v in moving.items() if k not in ("checks", "gates")
        },
        "whole_finite_local_classical_domain": "Fix finite L and a finite symmetric nonzero Fourier set Lambda. Each real spatial basis function has eight physical configuration channels, so d=8|Lambda| when Lambda counts both signs. Homogeneous reduced variables are external reference parameters, not quantized degrees of freedom. Reconstruct all nonlinear spatial harmonics with the S267 full convolution and exact adjoint ghost inverse. Choose a closed phase ball mapping strictly inside S266's twelve-invariant box and S267's chart at u=0. Smooth finite-parameter reconstruction, compactness, algebraic lapse/temporal constraints and their strict pivots give an unevaluated common positive time interval for the WHOLE u-dependent parent. This is a finite variational regulator, not an invariant Fourier submanifold of the nonlinear PDE.",
        "whole_interaction_picture_definition": "Let h_ref be the complete quadratic Hamiltonian generating the original fixed free preparation flow S_u, with S_0=I and its metaplectic propagator U_ref. Define g(u,z)=(h_red-h_ref)(u,S_u z). Retain the entire difference, including scalar, linear, mixed-quadratic and higher terms; the reference need not be the full nonlinear Hessian. The moving-chart connection occurs in h_red once, and the reference-flow connection subtracts h_ref once.",
        "whole_localized_symbols": "On a common interaction-picture phase ball choose time-independent smooth invariant chi with 0<=chi<=1, equal1 on a smaller core and compact support in the domain. Extend g_ext=g(u,0)+chi(z)[g(u,z)-g(u,0)]. Let F_red be the spatial average exp(3v) R_full(u,Nstar)^(-3/4), and F(u,z)=F_red(u,S_u z) its ENTIRE interaction-picture pullback. Define F_ext=F(u,0)+chi(z)[F(u,z)-F(u,0)]. The background factors are scalar identities, not discarded vacuum phases. F_ext>0 by convexity; neither extension is the original function outside the local domain.",
        "checks": checks,
        "gates": {
            **moving["gates"],
            "whole_original_all_spatial_matter_vector_channels_retained": all(
                full["whole_reduced_Hamiltonian"].has(z)
                for z in (
                    canonical.G,
                    canonical.shear,
                    canonical.electric,
                    canonical.magnetic,
                    canonical.wmass,
                    canonical.gm,
                    canonical.gh,
                    canonical.curvature,
                    canonical.j,
                )
            ),
            "full_time_source_not_bounce_germ_off_slice": True,
            "same_free_reference_without_reminimization_or_conditioning": True,
            "finite_time_and_canonical_radius_not_evaluated": True,
            "entire_canonical_invariant_second_contacts_retained": True,
            "finite_variational_regulator_not_full_original_PDE": True,
        },
    }
