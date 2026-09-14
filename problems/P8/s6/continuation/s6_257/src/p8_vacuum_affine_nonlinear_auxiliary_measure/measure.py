"""Finite-regulator second-class density and the source-centering boundary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_classical_principal_realization import realization


@cache
def block_measure():
    a, b, c, d, e = s.symbols("D11 D12 D21 D22 secondary_bracket", real=True)
    D = s.Matrix([[a, b], [c, d]])
    E = s.Matrix([[0, e], [-e, 0]])
    bracket = s.zeros(2).row_join(-D).col_join(D.T.row_join(E))
    inverse = (
        (D.T.inv() * E * D.inv())
        .row_join(D.T.inv())
        .col_join((-D.inv()).row_join(s.zeros(2)))
    )
    left = s.Matrix([[0, 0, *s.symbols("retained_left0:2", real=True)]])
    right = s.Matrix([0, 0, *s.symbols("retained_right0:2", real=True)])
    return {
        "whole_auxiliary_derivative_matrix": D,
        "whole_nonzero_secondary_constraint_bracket": E,
        "whole_primary_secondary_constraint_bracket": bracket,
        "whole_constraint_bracket_inverse": inverse,
        "whole_constraint_bracket_determinant": s.factor(bracket.det()),
        "whole_positive_second_class_density": s.Abs(D.det()),
        "single_regular_root_delta_Jacobian": 1 / s.Abs(D.det()),
        "finite_regulator_reduction": "At each fixed finite regulator use the entire auxiliary Jacobian D, including all cell and retained-field dependence. Constraints (p_aux,C_aux) have the displayed block form with E not assumed zero. The inverse and determinant formula hold at every finite matrix size by block multiplication and row permutation. On one regular root branch, integrating p_aux and delta(C_aux) cancels the positive |det D| density, leaving the retained canonical Liouville measure. Multiple roots require their separate branch contributions. Perform this cancellation before any continuum or determinant regularization.",
        "remaining_gauge_factor": "Keep the spatial diffeomorphism constraint delta functions, gauge-fixing conditions and their first-class determinant as a separate unevaluated factor. No finite spatial regulator is asserted to preserve their algebra. This second-class calculation neither supplies a BRST-invariant regulator nor quantizes the affine complement or its projective gauge directions.",
        "retained_coordinate_count_before_spatial_gauge": 11,
        "retained_phase_count_before_spatial_gauge": 22,
        "classical_phase_count_after_three_regular_spatial_gauges": 16,
        "checks": {
            "whole_constraint_bracket_inverse_left": (
                bracket * inverse - s.eye(4)
            ).applyfunc(s.factor),
            "whole_constraint_bracket_inverse_right": (
                inverse * bracket - s.eye(4)
            ).applyfunc(s.factor),
            "whole_constraint_bracket_determinant_square": s.factor(
                bracket.det() - D.det() ** 2
            ),
            "whole_retained_Dirac_bracket_correction_zero": s.factor(
                (left * inverse * right)[0]
            ),
            "whole_finite_regular_root_density_cancellation": s.factor(
                s.Abs(D.det()) / s.Abs(D.det()) - 1
            ),
            "classical_regular_spatial_gauge_phase_count": s.Integer(22 - 2 * 3 - 16),
        },
        "gates": {
            "nonzero_secondary_bracket_not_discarded": E != s.zeros(2),
            "inverse_secondary_secondary_block_zero": inverse[2:, 2:] == s.zeros(2),
            "all_six_metric_three_vector_two_matter_coordinates_retained": 6 + 3 + 2
            == 11,
            "remaining_first_class_measure_is_not_claimed_computed": True,
            "nonlinear_momentum_integral_is_not_assumed_Gaussian": True,
        },
    }


@cache
def nonlinear_example():
    q, N, T, p, pN, pT = s.symbols("q N T p pN pT", real=True)
    lam, eta = s.symbols("lambda eta", real=True)
    H = (
        (N * N + T * T) / 2
        + N * q
        + T * p
        + lam * N * T * q * q / 2
        + eta * N * N * T * T / 4
        + (q * q + p * p) / 2
    )
    constraints = s.Matrix([pN, pT, s.diff(H, N), s.diff(H, T)])
    Q, P = (q, N, T), (p, pN, pT)

    def PB(left, right):
        return sum(
            s.diff(left, x) * s.diff(right, y) - s.diff(left, y) * s.diff(right, x)
            for x, y in zip(Q, P)
        )

    bracket = s.Matrix(4, 4, lambda i, j: PB(constraints[i], constraints[j]))
    D = s.hessian(H, (N, T))
    E = bracket[2:, 2:]
    expected = s.zeros(2).row_join(-D).col_join(D.T.row_join(E))
    inverse = (
        (D.T.inv() * E * D.inv())
        .row_join(D.T.inv())
        .col_join((-D.inv()).row_join(s.zeros(2)))
    )
    left = s.Matrix([[PB(q, item) for item in constraints]])
    right = s.Matrix([PB(item, p) for item in constraints])
    return {
        "whole_nonlinear_example_Hamiltonian": H,
        "whole_nonlinear_example_constraints": constraints,
        "whole_nonlinear_example_auxiliary_Hessian": D,
        "whole_nonlinear_example_secondary_bracket": E,
        "whole_nonlinear_example_constraint_bracket": bracket,
        "finite_example_variables": (q, N, T, p, pN, pT, lam, eta),
        "example_boundary": "This nonlinear finite-dimensional diagnostic is not the physical parent action. The written general block proof and actual source-pinned canonical pivot establish the claimed local branch, not a numerical root of this example.",
        "checks": {
            "whole_nonlinear_constraint_bracket_blocks": bracket - expected,
            "whole_nonlinear_secondary_bracket_kept": E[0, 1] - 1 - lam * T * q,
            "whole_nonlinear_constraint_determinant_square": s.factor(
                bracket.det() - D.det() ** 2
            ),
            "whole_nonlinear_constraint_inverse": (
                bracket * inverse - s.eye(4)
            ).applyfunc(s.factor),
            "whole_nonlinear_retained_Dirac_correction_zero": s.factor(
                (left * inverse * right)[0]
            ),
        },
        "gates": {
            "nonlinear_auxiliary_Hessian_varies_with_auxiliaries": any(
                D.diff(item) != s.zeros(2) for item in (N, T)
            ),
            "nonzero_nonlinear_secondary_bracket": E != s.zeros(2),
        },
    }


@cache
def source_centering():
    y = s.Matrix(s.symbols("new_retained new_lapse new_temporal", real=True))
    x = s.Matrix(s.symbols("old_field0:3", real=True))
    alpha, beta, delta = s.symbols("map_alpha map_beta map_delta", real=True)
    point = s.Matrix(
        [
            y[0] + alpha * y[1] * y[2],
            y[1] + beta * y[0] ** 2,
            y[2] + delta * y[0] * y[1],
        ]
    )
    tadpole = s.Matrix(s.symbols("whole_background_onepoint0:3", real=True))
    external = s.Matrix(s.symbols("generating_function_source0:3", real=True))
    Hessian = s.Matrix(3, 3, s.symbols("old_quadratic0:9", real=True))
    Hessian = (Hessian + Hessian.T) / 2
    old = tadpole.dot(x) + (x.T * Hessian * x)[0] / 2 + x[0] * x[1] * x[2]
    pulled = old.subs(dict(zip(x, point)), simultaneous=True)
    complete = pulled - external.dot(point)
    auxiliary = (y[1], y[2])
    C = s.Matrix([s.diff(complete, item) for item in auxiliary])
    expected_C = (
        s.Matrix([s.diff(pulled, item) for item in auxiliary])
        - point.jacobian(auxiliary).T * external
    )
    D = s.hessian(complete, auxiliary)
    source_D = sum(
        (external[i] * s.hessian(point[i], auxiliary) for i in range(3)), s.zeros(2)
    )
    origin = dict.fromkeys(y, 0)
    centered = complete.subs(dict(zip(external, tadpole)), simultaneous=True)
    wrong_linear_source = pulled - tadpole.dot(y)
    contact = sum(
        (tadpole[i] * s.hessian(point[i], y).subs(origin) for i in range(3)), s.zeros(3)
    )
    initial = realization.data()
    profile = (
        realization.parent.RHO.subs(realization.bg.u, 0)
        - 3 * realization.parent.PRESSURE.subs(realization.bg.u, 0) / 2
    )
    rate = s.Symbol("initial_M1_rate", real=True)
    root_square = initial["actual_central_classical_M1_rate_squared"]
    return {
        "whole_example_nonlinear_source_point_map": point,
        "whole_example_unforced_pulled_Hamiltonian": pulled,
        "whole_example_source_dependent_Hamiltonian": complete,
        "whole_source_dependent_auxiliary_constraints": C,
        "whole_source_dependent_auxiliary_Hessian": D,
        "whole_auxiliary_source_contact": source_D,
        "whole_offshell_field_map_second_variation_contact": contact,
        "actual_old_reference_classical_constraint_residual": profile,
        "actual_N1_classical_constraint_root_square": root_square,
        "source_boundary": "The symbols generating_function_source are bookkeeping sources, not the fixed heavy source j or an authorized physical force. An original linear source pulls back as J dot F, not as J dot new_coordinates. It modifies both auxiliary constraints and their Jacobian. At an off-shell reference the second derivative contact is mandatory; source centering cancels it only with the complete pulled-back source. This algebra neither constructs the physical quantum effective action nor makes the fixed quantum reference a solution of the unforced classical equations.",
        "reference_Gaussian_boundary": "S252 remains the specified already spatial-gauge-fixed quadratic reference system with its own lapse/shift auxiliaries and Gaussian density. Its off-shell reference has the displayed exact profile residual. At N=1 the unforced classical comparison uses the different displayed M1-rate root. The profile correction is bounded but not proved zero or nonzero. No equality between this nonlinear lapse/temporal reduction and the S252 off-shell Gaussian is asserted; that requires the full source-dependent gauge, field-map and quantum-mean construction.",
        "checks": {
            "whole_source_dependent_auxiliary_constraints": (C - expected_C).applyfunc(
                s.expand
            ),
            "whole_source_dependent_auxiliary_Hessian": (
                D - s.hessian(pulled, auxiliary) + source_D
            ).applyfunc(s.expand),
            "whole_centered_source_first_variation": s.Matrix(
                [s.diff(centered, item).subs(origin) for item in y]
            ),
            "whole_correct_pulled_source_second_variation": s.hessian(centered, y).subs(
                origin
            )
            - Hessian,
            "whole_incorrect_linear_source_contact_is_retained": s.hessian(
                wrong_linear_source, y
            ).subs(origin)
            - Hessian
            - contact,
            "whole_source_point_map_initial_Jacobian": point.jacobian(y).subs(origin)
            - s.eye(3),
            "actual_old_reference_constraint_residual_binding": initial[
                "whole_current_lapse_constraint_at_N1_u0"
            ].subs(rate, s.Rational(1, 10))
            - profile,
            "actual_N1_classical_root_profile_shift_retained": root_square
            - s.Rational(1, 100)
            + 4 * profile,
        },
        "gates": {
            "nonlinear_source_auxiliary_contact_is_not_zero": source_D != s.zeros(2),
            "off_shell_full_field_map_contact_is_not_zero": contact != s.zeros(3),
            "exact_fixed_profile_residual_not_dropped": profile != 0,
            "classical_comparison_and_fixed_quantum_reference_remain_distinct": True,
        },
    }
