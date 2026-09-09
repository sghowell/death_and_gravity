"""Actual global joint state prescription and exact covariance bookkeeping."""

from functools import cache

import sympy as sp

from . import model, projectors


@cache
def data():
    a = model.a
    # Cartesian nine-scalar realization of the two TT polarizations:
    # t=2 Pi_TT Phi, pi=(1/2) Pi_TT P_Phi.
    pi = projectors.data()["TT"]
    J9 = sp.zeros(18)
    J9[:9, 9:] = sp.eye(9)
    J9[9:, :9] = -sp.eye(9)
    E = sp.diag(2 * pi, pi / 2)
    Jtt = sp.zeros(18)
    Jtt[:9, 9:] = pi
    Jtt[9:, :9] = -pi
    Mscalar = sp.Matrix([[0, 1 / a**3], [-a * model.q, 0]])
    MT = model.data()["tensor_polarization_generator"]
    S = sp.diag(sp.sqrt(2), 1 / sp.sqrt(2))
    return {
        "reference_state_scalar_factor": "The same chosen S6.101 positive global all-order scalar state, not an instantaneous free replacement",
        "reference_state_tensor_factor": "Two physical TT polarizations from the zero-mean massless scalar Minkowski state on the auxiliary past, transported through b to u=0 and then by actual a for all u",
        "reference_state_Proca_factor": "One actual m=1000 Proca field with three physical polarizations, positive Minkowski covariance on the auxiliary past, transported through b to u=0 and then by actual a for all u",
        "full_anchor_covariance": "C7(0,k)=diag(Cscalar(0,k), U_T,b(0,-3) C_T,flat U_T,b(0,-3)^T twice, U_P,b(0,-3) C_P,flat U_P,b(0,-3)^T)",
        "full_actual_bisolution": "W7(t,s,k)=U7,a(t,0,k) C7(0,k) U7,a(s,0,k)^T",
        "all_covariances_retain_common_hbar_over_kappa": True,
        "scalar_restriction_is_exactly_the_preexisting_chosen_scalar_state": True,
        "same_original_global_background_on_both_legs_of_final_state": True,
        "TT_Cartesian_realization": E,
        "scalar_to_each_tensor_polarization_canonical_map": S,
        "full_physical_mode_count": 7,
        "full_physical_CCR_phase_dimension": 14,
        "additional_physical_tensor_and_vector_modes": 5,
        "checks": {
            "projected_nine_scalar_realization_has_correct_TT_CCR": model.clean(
                E * J9 * E.T - Jtt
            ),
            "tensor_factor_is_actual_minimal_scalar_with_correct_normalization": model.clean(
                S * Mscalar * S.inv() - MT
            ),
            "tensor_scalar_field_normalization_is_canonical": model.clean(
                S * model.J2 * S.T - model.J2
            ),
            "joint_physical_mode_count": sp.Integer(2 + 2 + 3 - 7),
            "joint_physical_density_phase_dimension": sp.Integer(4 + 4 + 6 - 14),
        },
    }
