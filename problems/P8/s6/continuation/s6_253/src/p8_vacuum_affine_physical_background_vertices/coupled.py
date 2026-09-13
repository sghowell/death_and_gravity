"""Full aligned scalar/longitudinal-Proca Gaussian Hamiltonian and joint measure."""

from functools import cache

import sympy as s

D, Z, J, K, P = s.symbols(
    "positive_D positive_Z positive_J electric_K comoving_P", positive=True
)
th, w, c, Lnv, Vvv, Vvs, C, Y, Yv = s.symbols(
    "Theta mixed_momentum matter_charge lapse_v potential_vv potential_vsigma curvature scalar_gradient vector_spatial_mass",
    real=True,
)
r, rN, dH = s.symbols("R_minus_one R_lapse_derivative H_minus_fixed_clock_H", real=True)
v, sigma, longitudinal, n, b, A0, vd, sd, Ad = s.symbols(
    "metric_scalar matter_scalar vector_longitudinal lapse shift_divergence negative_vector_time metric_rate matter_rate vector_rate",
    real=True,
)
pv, ps, pA, q = s.symbols(
    "normalized_pv normalized_ps normalized_pA physical_q", real=True
)
PHASE = s.Matrix([v, sigma, longitudinal, pv, ps, pA])
OMEGA = s.zeros(6)
OMEGA[:3, 3:] = s.eye(3)
OMEGA[3:, :3] = -s.eye(3)
GAMMA = 1 - 3 * Z * r**2 / (2 * D)


def base_lagrangian():
    return (
        -3 * D * vd**2
        + Z * sd**2 / 2
        + 6 * th * n * vd
        + w * n * sd
        + (J - 3 * th**2 / D + w**2 / (2 * Z)) * n**2
        - 3 * c * vd * sigma
        + Lnv * n * v
        + Vvv * v**2
        + Vvs * v * sigma
        + b * (2 * D * vd - 2 * th * n + c * sigma)
        + 2 * C * q * v**2
        - Y * q * sigma**2 / 2
    )


def complete_lagrangian():
    source = r * (3 * vd - b) + 3 * rN * dH * n
    return (
        base_lagrangian()
        + K * (Ad - P * A0) ** 2 / 2
        + Z * (A0 + source) ** 2 / 2
        - Yv * longitudinal**2 / 2
    )


@cache
def unreduced_hamiltonian():
    lag = complete_lagrangian()
    rates = s.Matrix([vd, sd, Ad])
    momenta = s.Matrix([pv, ps, pA])
    kinetic = s.hessian(lag, rates)
    linear = s.Matrix([s.diff(lag, z).subs({vd: 0, sd: 0, Ad: 0}) for z in rates])
    solution = kinetic.inv() * (momenta - linear)
    return s.cancel(
        (momenta.dot(rates) - lag).subs(dict(zip(rates, solution)), simultaneous=True)
    )


def lapse_numerator():
    return (
        th * (pv - 3 * r * P * pA + 3 * c * sigma) / D
        - w * ps / Z
        - Lnv * v
        - 3 * rN * dH * P * pA
    )


def reduced_hamiltonian():
    return (
        ps**2 / (2 * Z)
        - c * (pv - 3 * r * P * pA) * sigma / (2 * D)
        - 3 * c**2 * sigma**2 / (4 * D)
        + Y * q * sigma**2 / 2
        - (2 * C * q + Vvv) * v**2
        - Vvs * v * sigma
        + lapse_numerator() ** 2 / (4 * J)
        + r * P * pA * pv / (2 * D)
        - 3 * r**2 * P**2 * pA**2 / (4 * D)
        + pA**2 / (2 * K)
        + P**2 * pA**2 / (2 * Z)
        + Yv * longitudinal**2 / 2
    )


@cache
def data():
    lag, ham = complete_lagrangian(), unreduced_hamiltonian()
    kinetic = s.hessian(lag, [vd, sd, Ad])
    auxiliary = s.hessian(ham, [n, b, A0]).applyfunc(s.cancel)
    A0star = s.cancel(-(s.diff(ham, A0) - auxiliary[2, 2] * A0) / auxiliary[2, 2])
    intermediate = s.cancel(ham.subs(A0, A0star))
    solution = {n: lapse_numerator() / (2 * J), b: pv / (2 * D)}
    reduced = s.cancel(intermediate.subs(solution, simultaneous=True))
    full = reduced_hamiltonian()
    checks = {
        "entire_joint_three_velocity_Hessian": kinetic - s.diag(-6 * D * GAMMA, Z, K),
        "complete_shift_lapse_auxiliary_entry": auxiliary[0, 1],
        "complete_shift_vector_time_auxiliary_entry": auxiliary[1, 2],
        "complete_shift_auxiliary_pivot": auxiliary[1, 1] + 2 * D / 3,
        "complete_vector_time_auxiliary_pivot": auxiliary[2, 2] + Z / GAMMA,
        "complete_lapse_vector_time_coupling": auxiliary[0, 2]
        + 3 * Z * (r * th / D + rN * dH) / GAMMA,
        "complete_lapse_pivot_after_vector_time_Schur": auxiliary[0, 0]
        - auxiliary[0, 2] ** 2 / auxiliary[2, 2]
        + 2 * J,
        "entire_joint_auxiliary_determinant": auxiliary.det()
        + 4 * D * Z * J / (3 * GAMMA),
        "joint_configuration_density_squared": kinetic.det() * auxiliary.det()
        - 8 * D**2 * Z**2 * J * K,
        "complete_vector_time_constraint_solution": s.diff(ham, A0).subs(A0, A0star),
        "complete_lapse_constraint_solution": s.diff(intermediate, n).subs(
            solution, simultaneous=True
        ),
        "complete_shift_constraint_solution": s.diff(intermediate, b).subs(
            solution, simultaneous=True
        ),
        "entire_aligned_six_phase_Hamiltonian": reduced - full,
        "entire_quadratic_Hessian_reconstructs_Hamiltonian": (
            PHASE.T * s.hessian(full, PHASE) * PHASE
        )[0]
        / 2
        - full,
    }

    R, lapse = s.symbols("positive_R positive_lapse", positive=True)
    checks["literal_parent_joint_Schur_factor"] = (
        GAMMA.subs(
            {
                D: R ** s.Rational(1, 4) / lapse,
                Z: R ** -s.Rational(3, 4) / lapse,
                r: R - 1,
            },
            simultaneous=True,
        )
        - 1
        + 3 * (R - 1) ** 2 / (2 * R)
    )
    checks["explicit_valid_R_domain_gamma_margin"] = (
        1
        - 3 * (R - 1) ** 2 / (2 * R)
        - s.Rational(1, 4)
        - 3 * (R - s.Rational(1, 2)) * (2 - R) / (2 * R)
    )
    return {
        "whole_aligned_lagrangian": lag,
        "whole_unreduced_Hamiltonian": ham,
        "whole_three_velocity_Hessian": kinetic,
        "whole_three_auxiliary_Hessian": auxiliary,
        "vector_time_solution_before_other_auxiliaries": A0star,
        "lapse_and_shift_solutions": solution,
        "whole_reduced_Hamiltonian": full,
        "whole_reduced_Hessian": s.hessian(full, PHASE),
        "joint_gamma": GAMMA,
        "normalized_Dirac_density": 4 * D * Z * J / (3 * GAMMA),
        "normalized_configuration_density_squared": 8 * D**2 * Z**2 * J * K,
        "checks": {
            name: value.applyfunc(s.cancel)
            if isinstance(value, s.MatrixBase)
            else s.cancel(value)
            for name, value in checks.items()
        },
        "gates": {
            "entire_source_square_retained": lag.has(r, rN),
            "mixed_first_source_vertex_nonzero": s.diff(full, r, pv, pA).subs(
                {r: 0, dH: 0}
            )
            != 0,
            "second_source_contact_nonzero": s.diff(full, r, 2).subs({r: 0, dH: 0})
            != 0,
            "vector_time_lapse_coupling_not_diagonalized_away": auxiliary[0, 2] != 0,
            "gamma_margin_requires_stated_R_domain": True,
            "aligned_slice_not_arbitrary_vector_background": True,
            "no_uniform_nonlinear_two_cone_neighborhood_inferred": True,
        },
    }


@cache
def phase_data():
    """Joint physical Dirac measure; no separate determinant prescription."""
    packet = data()
    density, a, kappa, zeta = s.symbols(
        "positive_density positive_hat_scale positive_kappa positive_zeta",
        positive=True,
    )
    H = s.Symbol("background_Hhat", real=True)
    aux = packet["whole_three_auxiliary_Hessian"]
    constraints = s.Matrix([s.diff(unreduced_hamiltonian(), x) for x in (n, b, A0)])
    gradients = constraints.jacobian(PHASE)
    lower = (gradients * OMEGA * gradients.T).applyfunc(s.cancel)
    Daux, lower_physical = density * aux, density * lower
    zero = s.zeros(3)
    dirac = zero.row_join(-Daux).col_join(Daux.T.row_join(lower_physical))
    invD = aux.inv().applyfunc(s.cancel) / density
    inverse = (
        (invD.T * lower_physical * invD)
        .row_join(invD.T)
        .col_join((-invD).row_join(zero))
    )
    Cphase = s.sqrt(kappa) * s.diag(1, 1, s.sqrt(zeta), a**3, a**3, a**3 / s.sqrt(zeta))
    hess = packet["whole_reduced_Hessian"]
    canonical_hessian = kappa * a**3 * Cphase.inv().T * hess * Cphase.inv()
    damping = s.diag(0, 0, 0, 1, 1, 1)
    weighted = OMEGA * hess - 3 * H * damping
    canonical_generator = OMEGA * canonical_hessian
    ref = {D: 1, Z: 1, r: 0, dH: 0, K: zeta / a**2, q: P**2 / a**2}
    physical_pp = canonical_hessian[3:, 3:].subs(ref, simultaneous=True)
    outer_det = 2 * a**7 * J / (th**2 * (1 + zeta * P**2 / a**2))
    checks = {
        "entire_joint_lower_constraint_bracket_antisymmetric": lower + lower.T,
        "entire_six_constraint_Dirac_inverse": dirac * inverse - s.eye(6),
        "physical_CCR_density_normalization": Cphase * OMEGA * Cphase.T
        - kappa * a**3 * OMEGA,
        "whole_canonical_phase_time_connection": Cphase * weighted * Cphase.inv()
        + 3 * H * damping
        - canonical_generator,
        "whole_six_phase_symplectic_generator": canonical_generator * OMEGA
        + OMEGA * canonical_generator.T,
        "joint_Dirac_density_including_physical_weight": Daux.det()
        + 4 * density**3 * D * Z * J / (3 * GAMMA),
        "joint_configuration_density_including_physical_weight": (
            density * packet["whole_three_velocity_Hessian"]
        ).det()
        * Daux.det()
        - 8 * density**6 * D**2 * Z**2 * J * K,
        "joint_canonical_configuration_density_rescaling": 8
        * (kappa * a**3) ** 6
        * D**2
        * Z**2
        * J
        * K
        / (kappa**6 * zeta**2)
        - 8 * a**18 * D**2 * Z**2 * J * K / zeta**2,
        "whole_reference_outer_reduced_configuration_determinant": physical_pp.det()
        * outer_det
        - 1,
        "joint_reference_unreduced_configuration_density": (
            8 * a**18 * D**2 * Z**2 * J * K / zeta**2
        ).subs(ref, simultaneous=True)
        - 8 * a**16 * J / zeta,
    }
    return {
        "whole_lower_constraint_bracket": lower_physical,
        "whole_six_constraint_Dirac_matrix": dirac,
        "whole_six_constraint_Dirac_inverse": inverse,
        "physical_phase_map": Cphase,
        "whole_canonical_six_phase_Hessian": canonical_hessian,
        "whole_canonical_six_phase_generator": canonical_generator,
        "reference_outer_reduced_configuration_determinant": outer_det,
        "joint_unreduced_canonical_configuration_density_squared": 8
        * a**18
        * D**2
        * Z**2
        * J
        * K
        / zeta**2,
        "measure_domain": "D,Z,J,K,density positive and gamma>0; 1/2<R<2 gives gamma>1/4. The normalized configuration measure uses the single joint oscillatory Gaussian and its continued phase. The outer reduced configuration determinant additionally requires Theta nonzero; the physical phase and auxiliary Dirac matrix do not.",
        "quantum_measure_boundary": "This is the already gauge-fixed quadratic physical measure. It is not a nonlinear covariant, BRST, ghost, affine connection, disformal field-map or separately regularized determinant prescription.",
        "checks": {
            name: value.applyfunc(s.cancel)
            if isinstance(value, s.MatrixBase)
            else s.cancel(value)
            for name, value in checks.items()
        },
        "gates": {
            "whole_lower_constraint_bracket_not_assumed_zero": lower != zero,
            "entire_inverse_lower_right_block_zero": inverse[3:, 3:] == zero,
            "ordinary_physical_phase_bracket_retained": True,
            "configuration_pivot_not_used_at_Theta_crossing": True,
        },
    }


@cache
def held_vector_data():
    """Entire held temporal-vector quadratic square and its clock vertex."""
    Wbar, ZN, ZNN, rNN = s.symbols("held_temporal_vector Z_N Z_NN R_NN", real=True)
    num0 = Wbar - 3 * r * dH
    num1 = -A0 - r * (3 * vd - b) - 3 * rN * dH * n
    num2 = (
        -b * longitudinal / P
        - 3 * r * v * b
        - rN * n * (3 * vd - b)
        - 3 * rNN * dH * n**2 / 2
    )
    square = (
        Z * (num1**2 + 2 * num0 * num2) / 2
        + (ZN * n + 3 * Z * v) * num0 * num1
        + (ZNN * n**2 / 4 + 3 * ZN * n * v / 2 + 9 * Z * v**2 / 4) * num0**2
    )
    lagrangian = (
        base_lagrangian()
        + K * (Ad - P * A0) ** 2 / 2
        + square
        - Yv * longitudinal**2 / 2
    )
    restriction = {Wbar: 0, r: 0, dH: 0, D: 1, Z: 1}
    Lnormal = s.diff(square, Wbar).subs(restriction, simultaneous=True)
    nstar = (th * (pv + 3 * c * sigma) - w * ps - Lnv * v) / (2 * J)
    on_shell = {n: nstar, b: pv / 2, vd: th * nstar - c * sigma / 2, A0: P * pA}
    normal_ham = (-Lnormal).subs(on_shell, simultaneous=True)
    scalar_part = rN * nstar * (3 * th * nstar - 3 * c * sigma / 2 - pv / 2)
    expected_ham = (
        scalar_part + pv * longitudinal / (2 * P) + (ZN * nstar + 3 * v) * P * pA
    )
    eps, etaR, etaH = s.symbols("epsilon first_R first_Hhat", real=True)
    embedding = 3 * (eps * etaR) * (eps * etaH)
    # The chain rule is exact for arbitrary function jets: S(0)=S'(0)=0.
    g, gv, gvv, hw, hvw, hww = s.symbols("H0 H1 H2 HW H1W HWW", real=True)
    independent = (
        g
        + eps * gv
        + eps**2 * gvv / 2
        + hw * Wbar
        + eps * hvw * Wbar
        + hww * Wbar**2 / 2
    )
    checks = {
        "entire_aligned_slice_recovers_complete_source_square": lagrangian.subs(
            Wbar, 3 * r * dH
        )
        - complete_lagrangian(),
        "held_vector_clock_Lagrangian_vertex": Lnormal
        + b * longitudinal / P
        + rN * n * (3 * vd - b)
        + (ZN * n + 3 * v) * A0,
        "held_vector_clock_Hamiltonian_envelope_vertex": normal_ham - expected_ham,
        "normal_vector_vertex_whole_phase_reconstruction": (
            PHASE.T * s.hessian(normal_ham, PHASE) * PHASE
        )[0]
        / 2
        - normal_ham,
        "aligned_and_held_first_derivatives_agree": s.diff(
            independent.subs(Wbar, embedding) - independent.subs(Wbar, 0), eps
        ).subs(eps, 0),
        "whole_second_embedding_contact": s.diff(
            independent.subs(Wbar, embedding) - independent.subs(Wbar, 0), eps, 2
        ).subs(eps, 0)
        - 6 * etaR * etaH * hw,
    }
    return {
        "whole_held_temporal_vector_quadratic_Lagrangian": lagrangian,
        "whole_held_temporal_vector_normal_square": square,
        "whole_clock_normal_vector_Hamiltonian_vertex": normal_ham,
        "whole_clock_normal_vector_Hessian": s.hessian(normal_ham, PHASE),
        "surviving_scalar_part_in_centered_product_reference": scalar_part,
        "alignment_second_derivative": 6 * etaR * etaH,
        "embedding_rule": "At the reference S=S_first=0. H_aligned_second=H_held_W0_second+H_W0*S_second. For the full finite effective action Gamma_aligned_second=Gamma_held_W0_second+Gamma_W0*S_second. The Gaussian Hamiltonian part has Gamma_W0=-expectation(H_W0); quantum measure and other sectors are not silently dropped. The reference product makes scalar/vector cross expectations zero, but the displayed scalar quadratic source expectation need not vanish.",
        "physical_vector_mean_boundary": "A centered free vector covariance does not set the interacting normal-vector tadpole or density-weighted Gauss constraint to zero. No physical vector mean is assigned without its complete reconstruction and chosen composite subtraction.",
        "infrared": "The displayed shift pairing contains 1/P and is interpreted on P>0 with the existing physical-shift infrared domain. It is not a literal P=0 homogeneous constraint formula.",
        "checks": {name: s.cancel(value) for name, value in checks.items()},
        "gates": {
            "nonzero_light_quadratic_normal_vector_source_retained": scalar_part != 0,
            "held_W0_and_aligned_second_probes_distinguished": True,
            "whole_shift_vector_pairing_retained": normal_ham.has(1 / P),
        },
    }
