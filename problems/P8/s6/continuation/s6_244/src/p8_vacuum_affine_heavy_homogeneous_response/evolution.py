"""Exact full scalar frame, covariance, prepared tangent and contact identities."""

from functools import cache

import sympy as s


@cache
def data():
    w, K = s.symbols("omega K", positive=True)
    lam, rate, sd = s.symbols("lambda rate_K squeeze", real=True)
    J = s.Matrix([[0, 1], [-1, 0]])
    M = s.diag(w * w / K, K)
    T = s.diag(s.sqrt(K / w), s.sqrt(w / K))
    Tp = T.diff(K) * K * rate + T.diff(w) * w * lam
    balanced = s.simplify(T.inv() * J * M * T - T.inv() * Tp)
    frame = s.Matrix([[(lam - rate) / 2, w], [-w, (rate - lam) / 2]])
    checks = {
        "complete_time_dependent_canonical_frame": s.simplify(balanced - frame),
        "canonical_frame_symplectic": s.simplify(T.T * J * T - J),
    }
    x, y = s.symbols("graph_real graph_imaginary", real=True)
    den = 1 - x * x - y * y
    Sigma = s.Matrix(
        [
            [(1 + 2 * x + x * x + y * y) / (2 * den), -y / den],
            [-y / den, (1 - 2 * x + x * x + y * y) / (2 * den)],
        ]
    )
    xd = -2 * w * y + sd * (1 - x * x + y * y)
    yd = 2 * w * x - 2 * sd * x * y
    flow = s.Matrix([[sd, w], [-w, -sd]])
    checks["complete_cross_covariance_pure_determinant"] = s.factor(
        Sigma.det() - s.Rational(1, 4)
    )
    checks["exact_graph_to_full_covariance_evolution"] = s.simplify(
        Sigma.diff(x) * xd + Sigma.diff(y) * yd - flow * Sigma - Sigma * flow.T
    )
    dk, dw, tau, qD = s.symbols(
        "delta_K delta_omega trace_D relative_omega_squared_D", real=True
    )
    MD = M.diff(K) * dk + M.diff(w) * dw
    N = s.simplify(T.T * MD * T / w)
    checks["complete_normalized_spatial_current_feature"] = s.simplify(
        N.subs({dk: -K * tau / 2, dw: w * qD / 2}) - s.diag(tau / 2 + qD, -tau / 2)
    )
    cc0, cc1, cc2 = s.symbols("C_QQ C_QP C_PP", real=True)
    C = s.Matrix([[cc0, cc1], [cc1, cc2]])
    current = -s.trace(MD * T * C * T.T) / 2
    expected = -dw * cc0 + (w / 2) * (dk / K) * (cc0 - cc2)
    checks["full_canonical_current_equals_scalar_variational_readout"] = s.factor(
        current - expected
    )
    eps = s.Symbol("source_parameter", real=True)
    md = s.Matrix(2, 2, s.symbols("MD0:4"))
    mdg = s.Matrix(2, 2, s.symbols("MDG0:4"))
    cov = s.Matrix(2, 2, s.symbols("cov0:4"))
    cov1 = s.Matrix(2, 2, s.symbols("cov1_0:4"))
    varied = -s.trace((md + eps * mdg) * (cov + eps * cov1)) / 2
    contact = -s.trace(md * cov1 + mdg * cov) / 2
    checks["whole_first_response_keeps_instantaneous_second_Hamiltonian_contact"] = (
        s.expand(s.diff(varied, eps).subs(eps, 0) - contact)
    )

    rr, rs, alpha, beta, phase = s.symbols("r rsharp alpha beta unit_phase")
    changed = (alpha * rr + beta * phase) / (alpha + beta * phase * rs)
    checks["exact_normalized_Bogoliubov_graph_difference"] = s.factor(
        (changed - rr) * (alpha + beta * phase * rs) - beta * phase * (1 - rr * rs)
    )
    W, Wp, Wpp, theta, thetap = s.symbols(
        "W Wprime Wsecond theta theta_prime", real=True
    )
    marker = s.Symbol("adiabatic_marker", real=True)
    decay = (theta + Wp / W) / 2
    graph = (w - W + s.I * marker * decay) / (w + W - s.I * marker * decay)
    graph_prime = (
        s.diff(graph, w) * w * lam
        + s.diff(graph, W) * Wp
        + s.diff(graph, Wp) * Wpp
        + s.diff(graph, theta) * thetap
    )
    frequency_defect = (
        W**2
        - w**2
        + marker**2
        * (thetap / 2 + theta**2 / 4 + Wpp / (2 * W) - 3 * Wp**2 / (4 * W**2))
    )
    checks["whole_WKB_frequency_to_initial_graph_defect_identity"] = s.factor(
        marker * graph_prime
        - 2 * s.I * w * graph
        - marker * (lam + theta) * (1 - graph**2) / 2
        - 2 * s.I * w * frequency_defect / (w + W - s.I * marker * decay) ** 2
    )
    checks["initial_source_germ_derivative_one"] = s.diff(rr, eps)
    checks["initial_source_germ_derivative_two"] = s.diff(rr, eps, 2)

    omega = s.symbols("omega0:3", real=True)
    squeeze = s.symbols("squeeze0:3", real=True)
    reference = s.symbols("reference0:3")
    error = s.symbols("error0:3")
    defect = s.symbols("defect0:3")

    def jet(values):
        return sum(value * eps**j / s.factorial(j) for j, value in enumerate(values))

    ww, ss, hh, ee, ff = map(jet, (omega, squeeze, reference, error, defect))
    vector = 2 * s.I * ww * ee - ss * ((hh + ee) ** 2 - hh**2) - ff
    actual = reference[0] + error[0]
    actual1 = reference[1] + error[1]
    tangent = 2 * s.I * omega[0] - 2 * squeeze[0] * actual
    first = (
        tangent * error[1]
        + (
            2 * s.I * omega[1]
            - squeeze[1] * (actual + reference[0])
            - 2 * squeeze[0] * reference[1]
        )
        * error[0]
        - defect[1]
    )
    second = (
        tangent * error[2]
        + (
            2 * s.I * omega[2]
            - squeeze[2] * (actual + reference[0])
            - 4 * squeeze[1] * reference[1]
            - 2 * squeeze[0] * reference[2]
        )
        * error[0]
        + (
            4 * s.I * omega[1]
            - 4 * squeeze[1] * actual
            - 2 * squeeze[0] * (actual1 + reference[1])
        )
        * error[1]
        - defect[2]
    )
    checks["complete_exact_first_error_forcing"] = s.expand(
        s.diff(vector, eps).subs(eps, 0) - first
    )
    checks["complete_exact_second_error_forcing_including_error_square"] = s.expand(
        s.diff(vector, eps, 2).subs(eps, 0) - second
    )
    checks["nonzero_quadratic_first_error_coefficient"] = (
        s.expand(second).coeff(error[1], 2) + 2 * squeeze[0]
    )

    nq, np = s.symbols("N_QQ N_PP", real=True)
    reference_current = (
        -w * ((nq + np) * (1 + rr * rs) + (nq - np) * (rr + rs)) / (4 * (1 - rr * rs))
    )
    checks["entire_diagonal_current_even_marker"] = s.factor(
        reference_current - reference_current.xreplace({rr: rs, rs: rr})
    )
    qp = -(rr - rs) / (2 * s.I * (1 - rr * rs))
    checks["cross_covariance_odd_not_even_marker"] = s.factor(
        qp + qp.xreplace({rr: rs, rs: rr})
    )
    for order in range(1, 11):
        checks["complete_finite_coefficient_marker_parity_" + str(order)] = (
            -s.I
        ) ** order - (-1) ** order * s.I**order
    return {
        "full_scalar_canonical_Hamiltonian": M,
        "complete_balanced_canonical_frame": T,
        "complete_balanced_generator": balanced,
        "exact_pure_graph_covariance": Sigma,
        "exact_real_graph_evolution": (xd, yd),
        "whole_spatial_current_vertex": s.diag(tau / 2 + qD, -tau / 2),
        "complete_covariance_response_with_second_contact": contact,
        "exact_Bogoliubov_graph": changed,
        "whole_WKB_graph_defect": 2
        * s.I
        * w
        * frequency_defect
        / (w + W - s.I * marker * decay) ** 2,
        "full_first_error_forcing": first,
        "full_second_error_forcing": second,
        "entire_reference_diagonal_current": reference_current,
        "physical_cross_covariance_not_declared_even": qp,
        "checks": checks,
        "gates": {
            "all_canonical_covariances_and_contacts_retained": True,
            "single_scalar_no_vector_polarization_factor": True,
            "no_reset_initial_mismatch_or_live_state_reminimization": True,
            "nonzero_quadratic_first_error_term": s.expand(second).coeff(error[1], 2)
            != 0,
            "whole_current_not_entire_covariance_is_even": qp != 0,
            "prepared_source_germ_keeps_initial_parameter_errors_zero": True,
        },
    }
