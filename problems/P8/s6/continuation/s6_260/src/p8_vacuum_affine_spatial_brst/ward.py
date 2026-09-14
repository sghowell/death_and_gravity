"""Conditional fixed-source Ward transport and exact finite orbit diagnostics."""

from functools import cache

import sympy as s


@cache
def orbit():
    v, t = s.symbols("positive_variance positive_gauge_parameter", positive=True)
    Jq, Jy = s.symbols("source_q source_y", real=True)
    phi, mean_y = s.symbols("mean_q mean_y", real=True)
    D = 1 - 2 * v * t * Jy
    W = -s.log(D) / 2 + v * Jq**2 / (2 * D)
    qmean = s.diff(W, Jq)
    qsecond = v / D + qmean**2
    ymean = s.diff(W, Jy)
    Legendre_D = v / (mean_y / t - phi**2)
    sources = {Jq: Legendre_D * phi / v, Jy: (1 - Legendre_D) / (2 * v * t)}
    Gamma = mean_y / (2 * v * t) - s.Rational(1, 2) + s.log(Legendre_D) / 2
    coordinates = s.Matrix([phi, mean_y])
    K = s.Matrix([0, mean_y / t])
    H = s.hessian(Gamma, coordinates)
    root = {phi: 0, mean_y: t * v}
    Kjac = K.jacobian(coordinates)
    moving_H = s.diff(H, t) + sum(
        (s.diff(H, coordinates[a]) * K[a] for a in range(2)), s.zeros(2)
    )
    generated_source = (
        s.Matrix([s.diff(Gamma, phi), s.diff(Gamma, mean_y)])
        - s.Matrix([sources[Jq], sources[Jy]])
    ).applyfunc(s.factor)
    q = s.Symbol("integration_q", real=True)
    moments = [
        s.integrate(
            q**n * s.exp(-q * q / (2 * v)) / s.sqrt(2 * s.pi * v), (q, -s.oo, s.oo)
        )
        for n in (0, 2, 4, 6)
    ]
    eta = s.Symbol("positive_boundary_weight", positive=True)
    defect_second = -eta * (moments[3] - moments[1] * moments[2])
    return {
        "physical_orbit_action": "q^2/(2v), with gauge translations of y; gauge chi=y-t q^2. The gauge determinant is one. Sources remain Jq q+Jy y before gauge reduction.",
        "variance": v,
        "gauge_parameter": t,
        "sources": s.Matrix([Jq, Jy]),
        "whole_connected_generator": W,
        "whole_source_domain": D,
        "whole_retained_first_and_second_q_moments": s.Matrix([qmean, qsecond]),
        "whole_gauge_coordinate_mean": ymean,
        "whole_two_source_Legendre_effective_action": Gamma,
        "whole_mean_coordinate_domain": mean_y / t - phi**2,
        "whole_Nielsen_vector": K,
        "whole_moving_stationary_mean": s.Matrix([0, t * v]),
        "whole_stationary_Hessian": H.subs(root).applyfunc(s.factor),
        "whole_boundary_weight_counterexample_second_variation": defect_second,
        "boundary_counterexample": "Multiplying by exp(-eta y^2/2) before the same gauge reduction changes the physical q distribution to exp[-q^2/(2v)-eta t^2 q^4/2]. Its normalized q^2 mean has second t derivative -12 eta v^3 at t=0. A noninvariant boundary weight cannot be dropped from a Ward identity.",
        "checks": {
            "whole_gauge_parameter_source_Ward_identity": s.factor(
                s.diff(W, t) - Jy * qsecond
            ),
            "whole_composite_gauge_coordinate_mean": s.factor(ymean - t * qsecond),
            "whole_two_source_Legendre_gradient": generated_source,
            "whole_two_source_mean_inverse": s.Matrix(
                [
                    s.factor(qmean.subs(sources, simultaneous=True) - phi),
                    s.factor(ymean.subs(sources, simultaneous=True) - mean_y),
                ]
            ),
            "whole_effective_action_Nielsen_identity": s.factor(
                s.diff(Gamma, t)
                + (s.Matrix([s.diff(Gamma, z) for z in coordinates]).T * K)[0]
            ),
            "whole_stationary_mean_equations": s.Matrix(
                [s.diff(Gamma, z).subs(root) for z in coordinates]
            ),
            "whole_on_shell_Hessian_congruence_transport": (
                moving_H + Kjac.T * H + H * Kjac
            )
            .subs(root)
            .applyfunc(s.factor),
            "whole_actual_Gaussian_integrals": s.Matrix(
                [
                    moments[0] - 1,
                    moments[1] - v,
                    moments[2] - 3 * v**2,
                    moments[3] - 15 * v**3,
                ]
            ),
            "whole_noninvariant_boundary_weight_defect": s.factor(
                defect_second + 12 * eta * v**3
            ),
        },
        "gates": {
            "finite_convergent_source_domain_D_positive_is_required": True,
            "Legendre_domain_t_positive_and_y_over_t_minus_phi_squared_positive": True,
            "gauge_coordinate_mean_moves_without_physical_state_reset": s.diff(t * v, t)
            != 0,
            "composite_mean_not_evaluated_at_zero_q_mean": v != 0,
            "boundary_weight_defect_is_nonzero": defect_second != 0,
            "finite_orbit_diagnostic_not_the_P8_state": True,
        },
    }


@cache
def differentiated_identity():
    x, y, t = s.symbols("mean_coordinate1 mean_coordinate2 gauge_parameter", real=True)
    z = (x, y)
    Gamma = s.Function("full_effective_action")(*z, t)
    K = s.Matrix([s.Function("full_Nielsen_vector" + str(i))(*z, t) for i in range(2)])
    E = s.Function("full_Ward_defect")(*z, t)
    first = s.Matrix([s.diff(Gamma, a) for a in z])
    H = s.hessian(Gamma, z)
    residual = s.diff(Gamma, t) + (first.T * K)[0] - E
    Kjac = K.jacobian(z)
    convective = s.diff(H, t) + sum(
        (s.diff(H, a) * K[i] for i, a in enumerate(z)), s.zeros(2)
    )
    contact = sum((s.hessian(K[a], z) * first[a] for a in range(2)), s.zeros(2))
    expected = convective + Kjac.T * H + H * Kjac + contact - s.hessian(E, z)
    T = s.Matrix([[s.exp(t), t], [0, s.exp(-t)]])
    reference = s.diag(2, -3)
    moving = T.inv().T * reference * T.inv()
    B = s.diff(T, t) * T.inv()
    return {
        "conditional_Ward_identity": "For a left odd differential and S+sPsi+J.F, define K^a=(i/hbar)<deltaPsi sF^a>. A measure/state/regulator/endpoint defect E remains explicit: deltaGamma+Gamma_,a K^a=E. No E=0 claim for the complete P8 quantum functional is supplied.",
        "full_formal_residual": residual,
        "whole_second_field_derivative_identity": expected,
        "whole_off_shell_onepoint_contact": contact,
        "mean_rule": "Only if E=0 and the appropriate stationary Hessian complement is invertible does the stationary mean move by delta phi=K. On that path, dH=-K_,phi^T H-H K_,phi. Off shell the displayed Gamma_,a K^a_,ij term remains; with defects its derivatives remain as well.",
        "finite_regular_transport_matrix": T,
        "finite_transport_Hessian": moving,
        "checks": {
            "whole_second_derivative_Nielsen_chain": (
                s.hessian(residual, z) - expected
            ).applyfunc(s.expand),
            "whole_regular_finite_Hessian_transport": (
                s.diff(moving, t) + B.T * moving + moving * B
            ).applyfunc(s.simplify),
            "whole_regular_transport_inertia_determinant": s.simplify(moving.det() + 6),
        },
        "gates": {
            "off_shell_onepoint_contact_nonzero": contact != s.zeros(2),
            "state_measure_and_endpoint_defect_not_asserted_zero": True,
            "regular_on_shell_gauge_transport_does_not_repair_inertia": True,
            "actual_P8_Nielsen_vector_and_quantum_mean_not_computed": True,
        },
    }
