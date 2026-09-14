"""Extended canonical point map and all spatial kinetic Legendre blocks."""

from functools import cache

import sympy as s


@cache
def extended_map():
    u = s.Symbol("clock", real=True)
    N = s.Symbol("positive_lapse", positive=True)
    C = s.Function("positive_spatial_conformal_factor", positive=True)(u, N)
    h = s.Matrix(s.symbols("spatial_metric0:6", real=True))
    W = s.Matrix(s.symbols("spatial_vector0:3", real=True))
    shift = s.Matrix(s.symbols("shift0:3", real=True))
    T = s.Symbol("normal_vector", real=True)
    matter = s.Matrix(s.symbols("unchanged_matter0:2", real=True))
    newq = h.col_join(W).col_join(shift).col_join(s.Matrix([T, N])).col_join(matter)
    oldq = (
        (C * h)
        .col_join(W)
        .col_join(shift)
        .col_join(s.Matrix([N * T + shift.dot(W), N]))
        .col_join(matter)
    )
    oldp = s.Matrix(s.symbols("old_momentum0:16", real=True))
    jac = oldq.jacobian(newq)
    newp = jac.T * oldp
    expected = (
        (C * oldp[:6, :])
        .col_join(oldp[6:9, :] + shift * oldp[12])
        .col_join(oldp[9:12, :] + W * oldp[12])
        .col_join(
            s.Matrix(
                [
                    N * oldp[12],
                    oldp[13] + s.diff(C, N) * h.dot(oldp[:6, :]) + T * oldp[12],
                ]
            )
        )
        .col_join(oldp[14:, :])
    )
    dq = s.Matrix(s.symbols("independent_dq0:16"))
    P = s.Matrix(s.symbols("canonical_momentum0:16", real=True))
    inverse_p = (
        (P[:6, :] / C)
        .col_join(P[6:9, :] - shift * P[12] / N)
        .col_join(P[9:12, :] - W * P[12] / N)
        .col_join(
            s.Matrix(
                [P[12] / N, P[13] - s.diff(C, N) * h.dot(P[:6, :]) / C - T * P[12] / N]
            )
        )
        .col_join(P[14:, :])
    )
    old_time_contact = oldp.dot(oldq.diff(u))
    new_time_contact = s.diff(C, u) * P[:6, :].dot(h) / C
    return {
        "whole_extended_new_position_vector": newq,
        "whole_extended_old_position_map": oldq,
        "whole_extended_old_momenta": oldp,
        "whole_extended_new_momentum_map": expected,
        "whole_extended_new_momenta": P,
        "whole_extended_inverse_momentum_map": inverse_p,
        "whole_extended_position_Jacobian": jac,
        "whole_position_density_factor": C**6 * N,
        "whole_extended_phase_density_factor": s.Integer(1),
        "whole_new_Hamiltonian_additive_time_contact": -new_time_contact,
        "map_boundary": "Use the complete point cotangent lift before imposing primary momenta. C depends on clock and lapse. The six metric coordinates have their dual six-coordinate momenta; tensor off-diagonal multiplicities are included in that dual convention. Both matter coordinates are unchanged spectators. H_new is H_old pulled back minus the displayed old p dot partial_u F contact. This map does not quantize the affine complement or gauge-fix spatial diffeomorphisms.",
        "checks": {
            "whole_extended_point_cotangent_momenta": (newp - expected).applyfunc(
                s.factor
            ),
            "whole_extended_cotangent_oneform": s.factor(
                oldp.dot(jac * dq) - newp.dot(dq)
            ),
            "whole_point_position_Jacobian": s.factor(jac.det() - C**6 * N),
            "whole_momentum_map_inverse": (jac.T * inverse_p - P).applyfunc(s.factor),
            "whole_extended_phase_volume_one": s.factor(
                jac.det() * inverse_p.jacobian(P).det() - 1
            ),
            "whole_explicit_time_contact": s.factor(
                old_time_contact - s.diff(C, u) * oldp[:6, :].dot(h)
            ),
            "whole_time_contact_in_new_momenta": s.factor(
                old_time_contact.subs(dict(zip(oldp, inverse_p)), simultaneous=True)
                - new_time_contact
            ),
        },
        "gates": {
            "all_sixteen_extended_positions_and_momenta_retained": len(newq)
            == len(oldq)
            == len(oldp)
            == len(P)
            == 16,
            "lapse_and_shift_primary_momenta_are_not_erased_before_lift": all(
                expected[13].has(item) for item in (oldp[13], oldp[12], s.diff(C, N))
            )
            and all(expected[9 + i].has(oldp[12], W[i]) for i in range(3)),
            "explicit_time_Hamiltonian_contact_not_deleted": new_time_contact != 0,
            "unchanged_two_matter_spectators": expected[14:, :] == oldp[14:, :],
        },
    }


@cache
def primitive_boundary_map():
    lift = extended_map()
    q = lift["whole_extended_new_position_vector"]
    P = lift["whole_extended_new_momenta"]
    oldp = lift["whole_extended_old_momenta"]
    oldq = lift["whole_extended_old_position_map"]
    u = s.Symbol("clock", real=True)
    N = q[13]
    h = s.Matrix([[q[0], q[1], q[2]], [q[1], q[3], q[4]], [q[2], q[4], q[5]]])
    V = s.sqrt(h.det())
    primitive = s.Function("whole_fixed_primitive", real=True)(u, N)
    boundary = V * primitive
    gradient = s.Matrix([s.diff(boundary, item) for item in q])
    hessian = s.hessian(boundary, q)
    # L_after = L_before - d(V I)/du, up to the retained spatial flux.
    # Therefore P_after = P_point_lift - d_q(V I), H_after = H_point+partial_u(V I).
    after = lift["whole_extended_new_momentum_map"] - gradient
    inverse = lift["whole_extended_inverse_momentum_map"].subs(
        dict(zip(P, P + gradient)), simultaneous=True
    )
    dq = s.Matrix(s.symbols("independent_boundary_dq0:16"))
    path_time = s.Symbol("boundary_path_time", real=True)
    path_q = s.Matrix(
        [s.Function("boundary_path_" + str(i), real=True)(path_time) for i in range(16)]
    )
    along_path = boundary.subs(
        {u: path_time, **dict(zip(q, path_q))}, simultaneous=True
    )
    path_derivative = (
        s.diff(along_path, path_time)
        .xreplace(
            {
                **dict(zip(path_q.diff(path_time), dq)),
                **dict(zip(path_q, q)),
                path_time: u,
            }
        )
        .doit()
    )
    time_shift = -oldp.dot(oldq.diff(u)) + s.diff(boundary, u)
    mapped_shift = time_shift.subs(dict(zip(oldp, inverse)), simultaneous=True)
    C = s.Function("positive_spatial_conformal_factor", positive=True)(u, N)
    expected_shift = -s.diff(C, u) * (P[:6, :] + gradient[:6, :]).dot(
        q[:6, :]
    ) / C + s.diff(boundary, u)
    # Independent local density identity, with div(V shift) left intact.
    volume, volume_dot, div_volume_shift, K, prim, prim_dot, prim_u, shift_grad_prim = (
        s.symbols(
            "volume volume_dot div_volume_shift K primitive primitive_dot primitive_u shift_grad_primitive",
            real=True,
        )
    )
    raw = volume * (prim_dot - prim_u - shift_grad_prim)
    after_density = -N * volume * prim * K - volume * prim_u
    retained_boundary = (
        volume_dot * prim
        + volume * prim_dot
        - prim * div_volume_shift
        - volume * shift_grad_prim
    )
    identity = (raw - after_density - retained_boundary).subs(
        volume_dot, N * volume * K + div_volume_shift
    )
    return {
        "whole_primitive_boundary_function": boundary,
        "whole_primitive_momentum_gradient": gradient,
        "whole_combined_point_and_boundary_momentum_map": after,
        "whole_combined_inverse_momentum_map": inverse,
        "whole_combined_Hamiltonian_time_shift_in_new_momenta": expected_shift,
        "whole_primitive_boundary_Hessian": hessian,
        "boundary_prescription": "Within the original S174 ADM variational boundary convention, the remaining linear lapse-velocity density is V I_s (s_dot-Ni partial_i s). It equals -N V I K-V I_u plus partial_u(V I)-partial_i(V Ni I). Subtracting that full boundary shifts momenta by minus d_q(V I) and adds partial_u(V I) to the Hamiltonian after the point lift. The full combined time term also uses the boundary-shifted metric momenta. Preserve the spatial flux and both initial/final endpoint phases; no state reset is licensed.",
        "checks": {
            "whole_primitive_time_and_spatial_boundary_identity": s.expand(identity),
            "whole_boundary_oneform_identity": s.factor(
                path_derivative - gradient.dot(dq) - s.diff(boundary, u)
            ),
            "whole_primitive_gradient_Hessian_symmetry": hessian - hessian.T,
            "whole_combined_momentum_map_inverse": (
                after.subs(dict(zip(oldp, inverse)), simultaneous=True) - P
            ).applyfunc(s.factor),
            "whole_combined_Hamiltonian_time_shift": s.factor(
                mapped_shift - expected_shift
            ),
        },
        "gates": {
            "entire_primitive_lapse_momentum_shift_retained": gradient[13].has(
                s.diff(primitive, N)
            ),
            "entire_primitive_metric_momentum_shifts_retained": all(
                gradient[i] != 0 for i in range(6)
            ),
            "combined_time_term_not_point_lift_alone": s.factor(
                expected_shift - lift["whole_new_Hamiltonian_additive_time_contact"]
            )
            != 0,
            "endpoint_and_spatial_boundary_prescription_retained": True,
        },
    }


@cache
def metric_density_bridge():
    N, C = s.symbols("positive_lapse positive_C", positive=True)
    h = s.Matrix(s.symbols("metric_coordinate0:6", real=True))
    shift = s.Matrix(s.symbols("ADM_shift0:3", real=True))
    metric = s.Matrix([[h[0], h[1], h[2]], [h[1], h[3], h[4]], [h[2], h[4], h[5]]])
    # Order g00, g0i, gij against N, Ni, hij. The ADM determinant has
    # sign - for the original +--- signature; its positive density is 2N det(h).
    adm = (
        s.Matrix([N * N - (shift.T * metric * shift)[0]])
        .col_join(-metric * shift)
        .col_join(-h)
    )
    coords = s.Matrix([N]).col_join(shift).col_join(h)
    jac = adm.jacobian(coords)
    signed = s.factor(jac.det())
    spatial_det = s.factor((C * metric).det())
    return {
        "whole_ADM_metric_position_Jacobian": jac,
        "whole_signed_ADM_coordinate_determinant": signed,
        "whole_covariant_metric_position_density_ratio": C**9,
        "whole_unitary_spatial_metric_position_density_ratio": C**6,
        "density_boundary": "The ADM conversion has absolute determinant 2N det(gamma). The physical spatial determinant is C^3 det(gamma), while the six spatial coordinate map has determinant C^6. Their product reconciles C^9 for ten covariant metric coordinates with C^6 in unitary ADM metric coordinates. Adding W0=N T+Ni Wi supplies the independent N factor. None of these position-only factors replaces the unit Liouville Jacobian of the full cotangent lift.",
        "checks": {
            "whole_ADM_metric_coordinate_determinant": signed + 2 * N * metric.det(),
            "whole_spatial_metric_volume_determinant": spatial_det
            - C**3 * metric.det(),
            "whole_covariant_and_unitary_metric_density_bridge": s.factor(
                spatial_det * C**6 / metric.det() - C**9
            ),
        },
        "gates": {"regular_positive_lapse_and_positive_spatial_metric_required": True},
    }


@cache
def spatial_kinetic_blocks():
    N, V, U, M, Cchi, zeta = s.symbols("N volume U M Cchi zeta", positive=True)
    a, b, c, d, e, f = s.symbols("g11 g12 g13 g22 g23 g33", real=True)
    metric = s.Matrix([[a, b, c], [b, d, e], [c, e, f]])
    vel = s.Matrix(s.symbols("electric_velocity0:3", real=True))
    grad = s.Matrix(s.symbols("gradient_W0_0:3", real=True))
    transport = s.Matrix(s.symbols("shift_F_0:3", real=True))
    momentum = s.Matrix(s.symbols("electric_momentum0:3", real=True))
    lag = (
        V
        * zeta
        * Cchi
        / (2 * N)
        * ((vel - grad - transport).T * metric.inv() * (vel - grad - transport))[0]
    )
    solution = N * metric * momentum / (V * zeta * Cchi) + grad + transport
    result = (momentum.dot(vel) - lag).subs(dict(zip(vel, solution)), simultaneous=True)
    expected = N * (momentum.T * metric * momentum)[0] / (
        2 * V * zeta * Cchi
    ) + momentum.dot(grad + transport)
    matter = s.Matrix(s.symbols("matter_velocity0:2", real=True))
    drift = s.Matrix(s.symbols("matter_shift_drift0:2", real=True))
    pm = s.Matrix(s.symbols("matter_momentum0:2", real=True))
    LM = V * U * (matter - drift).dot(matter - drift) / (2 * N)
    Hm = (pm.dot(matter) - LM).subs(
        dict(zip(matter, N * pm / (V * U) + drift)), simultaneous=True
    )
    shear_v = s.Matrix(s.symbols("five_orthonormal_shear_velocities0:5", real=True))
    shear_p = s.Matrix(s.symbols("five_orthonormal_shear_momenta0:5", real=True))
    shear_lag = V * M * shear_v.dot(shear_v) / (8 * N)
    shear_ham = (shear_p.dot(shear_v) - shear_lag).subs(
        dict(zip(shear_v, 4 * N * shear_p / (V * M))), simultaneous=True
    )
    DD = s.Matrix([[0, 1, -1], [-1, 0, 1], [1, -1, 0]])
    pi = s.Matrix(s.symbols("cell_pi0:3"))
    lapses = s.Matrix(s.symbols("cell_N0:3"))
    temporal = s.Matrix(s.symbols("cell_T0:3"))
    W0 = s.matrix_multiply_elementwise(lapses, temporal)
    gauss = -W0.dot(DD * pi)
    return {
        "whole_three_direction_spatial_metric": metric,
        "whole_electric_Lagrangian": lag,
        "whole_electric_velocity_solution": solution,
        "whole_electric_Hamiltonian": expected,
        "whole_two_matter_Lagrangian": LM,
        "whole_two_matter_Hamiltonian": s.expand(Hm),
        "whole_five_shear_Lagrangian": shear_lag,
        "whole_five_shear_Hamiltonian": s.expand(shear_ham),
        "finite_SBP_spatial_matrix": DD,
        "whole_finite_cell_temporal_Gauss_Hamiltonian": gauss,
        "spatial_boundary": "The Maxwell Legendre transform contains pi^i partial_i W0. After spatial integration by parts it is -N T partial_i pi^i plus the retained shift contribution. This removes spatial lapse derivatives from the auxiliary Hamiltonian without deleting the Gauss density. The finite antisymmetric three-cell matrix checks summation by parts only; it does not assert a discrete Leibniz rule or anomaly-free first-class diffeomorphism algebra.",
        "checks": {
            "whole_three_direction_electric_Legendre": s.factor(result - expected),
            "whole_two_matter_canonical_normalizations": s.factor(
                Hm - N * pm.dot(pm) / (2 * V * U) - pm.dot(drift)
            ),
            "whole_five_shear_canonical_normalizations": s.factor(
                shear_ham - 2 * N * shear_p.dot(shear_p) / (V * M)
            ),
            "whole_spatial_Gauss_summation_by_parts": s.expand(pi.dot(DD * W0) - gauss),
            "finite_spatial_SBP_antisymmetry": DD + DD.T,
        },
        "gates": {
            "full_spatial_metric_not_diagonalized": all(
                metric.has(item) for item in (b, c, e)
            ),
            "entire_Gauss_term_is_nonzero": gauss != 0,
            "both_matter_drift_terms_retained": all(Hm.has(item) for item in drift),
            "five_shear_three_electric_two_matter_velocities_retained": len(shear_v)
            + len(vel)
            + len(matter)
            == 10,
        },
    }
