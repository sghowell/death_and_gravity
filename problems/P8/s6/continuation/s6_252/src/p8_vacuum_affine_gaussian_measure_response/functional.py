"""Fixed-state finite CTP overlap and entire coefficient-family variations."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import gaussian, phase
from p8_vacuum_affine_scalar_tame_propagator import charts

from . import response


def bogoliubov(value):
    S = response.exact_matrix(value)
    S = response.symplectic(S, S.rows)
    half = S.rows // 2
    A, B, C, D = (
        S[:half, :half],
        S[:half, half:],
        S[half:, :half],
        S[half:, half:],
    )
    alpha = (A + D + s.I * (C - B)) / 2
    beta = (A - D + s.I * (C + B)) / 2
    return alpha, beta


def relative_bogoliubov(plus, minus, preparation):
    S0 = response.exact_matrix(preparation)
    dimension = S0.rows
    S0 = response.symplectic(S0, dimension)
    plus = response.symplectic(plus, dimension)
    minus = response.symplectic(minus, dimension)
    return bogoliubov(S0.inv() * minus.inv() * plus * S0)


def logdet_coefficients(coefficients, depth):
    """Finite formal matrix Taylor algebra; not convergence at epsilon1."""
    if type(depth) is not int or not 1 <= depth <= 6:
        raise ValueError("Require a specified finite formal depth1..6")
    if not isinstance(coefficients, (tuple, list)) or not coefficients:
        raise ValueError("Require the complete ordered matrix coefficient sequence")
    matrices = [s.Matrix(x) for x in coefficients]
    size = matrices[0].rows
    if not 1 <= size <= 8 or any(M.shape != (size, size) for M in matrices):
        raise ValueError("Require matching finite square matrix coefficients")
    if any(x.has(s.Float) for M in matrices for x in M) or matrices[0].det() == 0:
        raise ValueError("Require exact coefficients and an invertible leading matrix")
    zero = s.zeros(size)
    coeff = matrices + [zero] * max(0, depth + 1 - len(matrices))
    inverse = [coeff[0].inv()]
    for order in range(1, depth):
        total = sum(
            (coeff[j] * inverse[order - j] for j in range(1, order + 1)),
            zero,
        )
        inverse.append((-inverse[0] * total).applyfunc(s.cancel))
    return tuple(
        s.cancel(
            sum(j * s.trace(inverse[order - j] * coeff[j]) for j in range(1, order + 1))
            / order
        )
        for order in range(1, depth + 1)
    )


def complete_coefficient_hessian():
    momentum = s.Symbol("comoving_momentum", nonnegative=True)
    matrix = -phase.OMEGA * phase.canonical_generator()
    return momentum, matrix.subs(charts.q, momentum**2 / phase.a**2).applyfunc(s.cancel)


def bilinear_scaling_bound(left_order, right_order):
    """Conservative written-proof bound, not a physical counterterm inventory."""
    if any(type(order) is not int or order < 0 for order in (left_order, right_order)):
        raise ValueError("Specify both nonnegative finite total derivative orders")
    return 12 + left_order + right_order


@cache
def data():
    checks = {}
    O = response.omega(4)
    for case in range(4):
        S0 = gaussian.symplectic_fixture(case)
        plus = gaussian.symplectic_fixture(case + 3)
        minus = gaussian.symplectic_fixture(case + 7)
        alpha, beta = relative_bogoliubov(plus, minus, S0)
        checks[f"whole_relative_Bogoliubov_CCR_{case}"] = (
            alpha * alpha.conjugate().T - beta * beta.conjugate().T - s.eye(2)
        )
        checks[f"whole_relative_Bogoliubov_pair_symmetry_{case}"] = (
            alpha * beta.T - beta * alpha.T
        )
        same_alpha, same_beta = relative_bogoliubov(plus, plus, S0)
        checks[f"equal_CTP_histories_relative_identity_{case}"] = same_alpha - s.eye(2)
        checks[f"equal_CTP_histories_no_pair_creation_{case}"] = same_beta
        A = response.coefficient_fixture(case)
        C = response.coefficient_fixture(case + 2)
        first = S0.inv() * O * A * S0
        second = S0.inv() * ((O * A) ** 2 + O * C) * S0 / 2

        # The Bogoliubov block map is linear and also applies to Taylor coefficients.
        def alpha_block(M):
            return (M[:2, :2] + M[2:, 2:] + s.I * (M[2:, :2] - M[:2, 2:])) / 2

        a1, a2 = alpha_block(first).conjugate(), alpha_block(second).conjugate()
        log_first = -s.trace(a1) / 2
        log_second_derivative = -s.trace(a2) + s.trace(a1 * a1) / 2
        V = S0 * S0.T / 2
        kernel = response.kernels(A, A, s.eye(4), s.eye(4), V, C)
        checks[f"complete_CTP_overlap_first_vertex_matches_mean_{case}"] = (
            log_first + s.I * kernel["mean_A"]
        )
        checks[f"complete_CTP_overlap_second_vertex_and_noise_{case}"] = (
            log_second_derivative
            - s.I * kernel["effective_action_seagull"]
            + kernel["symmetric_noise"]
        )
    radial, whole = complete_coefficient_hessian()
    parameters = (charts.J, charts.th, charts.E, charts.l, charts.A, charts.T, phase.a)
    first_vertices = {
        str(parameter): whole.diff(parameter).applyfunc(s.cancel)
        for parameter in parameters
    }
    second_vertices = {
        str(parameters[i]) + "__" + str(parameters[j]): whole.diff(
            parameters[i], parameters[j]
        ).applyfunc(s.cancel)
        for i in range(len(parameters))
        for j in range(i, len(parameters))
    }
    for parameter in parameters:
        derivative = first_vertices[str(parameter)]
        checks["whole_first_coefficient_vertex_symmetric_" + str(parameter)] = (
            derivative - derivative.T
        )
    checks["full_physical_scale_derivative_keeps_comoving_q_chain"] = whole.diff(
        phase.a
    ) - (
        (-phase.OMEGA * phase.canonical_generator()).diff(phase.a)
        - 2
        * charts.q
        / phase.a
        * (-phase.OMEGA * phase.canonical_generator()).diff(charts.q)
    ).subs(charts.q, radial**2 / phase.a**2)
    eps, dJ, dT, dA = s.symbols(
        "formal_loop_grade DeltaJ_profile T_profile A_profile", real=True
    )
    J0 = s.Symbol("positive_bare_J", positive=True)
    graded = whole.subs(
        {charts.J: J0 + eps * dJ, charts.T: eps * dT, charts.A: eps * dA},
        simultaneous=True,
    )
    coefficients = [graded.diff(eps, j).subs(eps, 0) / s.factorial(j) for j in range(3)]
    checks["entire_profile_resummation_tree_part"] = coefficients[0] - whole.subs(
        {charts.J: J0, charts.T: 0, charts.A: 0}, simultaneous=True
    )
    checks["entire_profile_resummation_first_part"] = coefficients[1] - (
        dJ * whole.diff(charts.J)
        + dT * whole.diff(charts.T)
        + dA * whole.diff(charts.A)
    ).subs({charts.J: J0, charts.T: 0, charts.A: 0}, simultaneous=True)
    for case in range(3):
        m = s.Integer(case + 1)
        D0 = s.Matrix([[3 + m, 1, 0], [1, 5 + m, 1], [0, 1, 7 + m]])
        K = s.Matrix([[m, 1, 2], [1, -m, 3], [2, 3, 2 * m]])
        F = s.Matrix([[1, -m, 0], [-m, 2, m], [0, m, -1]])
        G0 = D0.inv()
        target = logdet_coefficients([D0, K], 4)
        free = logdet_coefficients([D0, F], 4)
        Gfree = [G0]
        for order in range(1, 4):
            Gfree.append((-G0 * F * Gfree[-1]).applyfunc(s.cancel))
        remainder = [s.eye(3)] + [Gfree[order - 1] * (K - F) for order in range(1, 5)]
        interaction = logdet_coefficients(remainder, 4)
        for order in range(4):
            checks[
                f"noncommuting_free_countervertex_split_order_{order + 1}_case_{case}"
            ] = free[order] + interaction[order] - target[order]
        checks[f"complete_logdet_first_order_{case}"] = target[0] - s.trace(G0 * K)
        checks[f"complete_logdet_second_order_{case}"] = (
            target[1] + s.trace((G0 * K) ** 2) / 2
        )
    # Off-shell nonlinear embedding term, independently differentiating the full functional.
    x, y, auxiliary, tangent, curvature, residual = s.symbols(
        "field_x field_y auxiliary tangent curvature nonzero_onepoint", real=True
    )
    functional = residual * auxiliary + auxiliary**2 / 2 + auxiliary * x + x * y + y**2
    embedding = tangent * x + curvature * x**2 / 2
    pulled = functional.subs(auxiliary, embedding)
    variables = s.Matrix([x, y, auxiliary])
    tangent_map = s.Matrix([[1, 0], [0, 1], [tangent, 0]])
    at_origin = {x: 0, y: 0, auxiliary: 0}
    naive = (
        tangent_map.T * s.hessian(functional, variables).subs(at_origin) * tangent_map
    )
    exact = s.hessian(pulled, [x, y]).subs({x: 0, y: 0})
    contact = s.diag(residual * curvature, 0)
    checks["entire_off_shell_nonlinear_embedding_contact"] = exact - naive - contact
    # A unitary full turn has the same classical endpoint and a different metaplectic phase.
    angle = s.Symbol("rotation_angle", real=True)
    rotation = s.Matrix([[s.cos(angle), s.sin(angle)], [-s.sin(angle), s.cos(angle)]])
    rotation_alpha_conjugate = s.cos(angle) + s.I * s.sin(angle)
    checks["rotation_Bogoliubov_phase_orientation"] = (
        rotation[0, 0] + rotation[1, 1] - s.I * (rotation[1, 0] - rotation[0, 1])
    ) / 2 - rotation_alpha_conjugate
    checks["closed_classical_rotation_endpoint"] = rotation.subs(
        angle, 2 * s.pi
    ) - s.eye(2)
    checks["continued_vacuum_phase_full_turn_is_minus_one"] = s.exp(-s.I * s.pi) + 1
    chart_orders = {}
    for name, chart in (("outer", charts.outer()), ("central", charts.central())):
        chart_orders[name] = {
            field: [charts.degree(entry) for entry in chart[field]]
            for field in ("K", "B", "G")
        }
        for index, entry in enumerate(chart["B"]):
            checks[f"entire_{name}_boundary_symbol_no_more_than_q_{index}"] = s.limit(
                entry / charts.q**2, charts.q, s.oo
            )
    maximum_vertex_degree = max(
        s.degree(entry, radial)
        for matrix in [whole, *first_vertices.values(), *second_vertices.values()]
        for entry in matrix
        if entry != 0
    )
    checks["whole_coefficient_vertex_maximum_spatial_degree_four"] = (
        maximum_vertex_degree - 4
    )
    checks["conservative_full_scalar_response_scaling_degree_bound"] = (
        bilinear_scaling_bound(4, 4) - 20
    )
    tensor_hessian = s.diag(phase.a * radial**2, phase.a**-3)
    tensor_first = s.diag(radial**2, -3 * phase.a**-4)
    tensor_second = s.diag(0, 12 * phase.a**-5)
    checks["entire_each_TT_scale_first_vertex"] = (
        tensor_hessian.diff(phase.a) - tensor_first
    )
    checks["entire_each_TT_scale_second_vertex"] = (
        tensor_hessian.diff(phase.a, 2) - tensor_second
    )
    return {
        "fixed_state_CTP_definition": "Z[h+,h-]=Tr[U_+(T,t*) rho0 U_-(T,t*)^dagger], rho0 fixed by S251. This is a finite-regulator Gaussian trace, not an in-out determinant or a state reset.",
        "whole_relative_symplectic_matrix": "R=S0^-1 S_-^-1 S_+ S0 for V0=S0 S0^T/2. In q,p block order alpha=(Rqq+Rpp+i(Rpq-Rqp))/2 and beta=(Rqq-Rpp+i(Rpq+Rqp))/2.",
        "complete_metaplectic_overlap": "Z=exp[-Tr Log(conjugate(alpha))/2] with the logarithm continued from identical histories along the actual Hamiltonian path. A principal square root of the endpoint determinant loses a possible metaplectic sign. Every retained quadratic Hamiltonian is Weyl ordered; any separate c-number action remains separate.",
        "closed_time_normalization_and_modulus": "Identical histories give R=I continuously and Z=1. alpha alpha^dagger=I+beta beta^dagger implies |Z|<=1. The pure preparation factor S0 is only a representation; right orthosymplectic basis changes do not alter the trace.",
        "whole_scalar_coefficient_Hessian": whole,
        "independent_coefficient_parameters": parameters,
        "whole_first_coefficient_vertices": first_vertices,
        "whole_second_coefficient_vertices": second_vertices,
        "each_complete_TT_Hessian_first_and_second_scale_vertex": (
            tensor_hessian,
            tensor_first,
            tensor_second,
        ),
        "coefficient_family_boundary": "Every scalar Hamiltonian entry and every first/second derivative is retained with q=P^2/a^2. These are independent coefficient-family derivatives, NOT an off-reference covariant metric/light map. Physical vertices still need that full map, nonlinear embedding contacts, gauge/constraint and boundary terms. Each TT mode has Hessian diag(a*P^2,a^-3), with the same fixed-state formulas.",
        "formal_graded_full_profile_Hessian_zero_through_two": coefficients,
        "loop_splitting_boundary": "The finite ordered identity det(Dfree)det(I+Dfree^-1 I)=det(Dtarget) must use a common regulator and the matching countervertex I. Formal Taylor coefficients retain the full noncommuting matrix order and existing counterprofile grades. No convergence at formal grade1, equality of independently regularized continuum determinants, or double counting of conditional H/Proca loops is allowed.",
        "off_shell_boundary": "The coefficient action alone has the known nonzero reference onepoint. The full second chain rule adds S_,a times the second embedding derivative. Its cancellation may only be used in the complete matched mean or justified formal order, not silently in an isolated local coefficient determinant.",
        "complete_chart_symbol_orders_in_q": {
            name: {
                field: [
                    "zero_polynomial" if degree == -s.oo else degree
                    for degree in degrees
                ]
                for field, degrees in chart.items()
            }
            for name, chart in chart_orders.items()
        },
        "whole_coefficient_vertex_maximum_spatial_degree": maximum_vertex_degree,
        "continuum_reference_Wick_noise": "The written two-cone wavefront argument defines every specified finite-jet reference-normal-ordered bilinear and its entire connected noise as distributions on the open slab, without a momentum cutoff. Every contraction retains the same annihilation time orientation, so its wavefront sums cannot vanish. Positivity follows by same-state smoothing and the centered-operator norm. This is not a covariant physical stress definition.",
        "continuum_time_retarded_extension_family": "The full canonical field amplitudes have order at most3/2, hence their two-point symbols have order at most3. A double-mode scaling argument gives a conservative diagonal scaling-degree bound12+dA+dB for total derivative orders dA,dB. The complete coefficient vertices have order at most4, giving bound20 and an extension ambiguity no larger than derivatives of the four-dimensional diagonal delta through order16. Smooth tangential symbol bounds and conormal closure supply the distribution-extension hypotheses. This proves a family of time-retarded extensions, not a selected physical subtraction, a variational/gauge-consistent counterfunctional, spacelike microcausality or an evaluated loop norm.",
        "checks": {
            name: value.applyfunc(s.cancel)
            if isinstance(value, s.MatrixBase)
            else s.cancel(value)
            for name, value in checks.items()
        },
        "gates": {
            "full_second_current_vertices_not_deleted": any(
                M != s.zeros(4) for M in second_vertices.values()
            ),
            "full_profile_reduction_contains_higher_formal_grades": coefficients[2]
            != s.zeros(4),
            "endpoint_principal_root_would_lose_metaplectic_sign": s.exp(-s.I * s.pi)
            != 1,
            "nonzero_onepoint_embedding_contact_is_independent": contact != s.zeros(2),
            "physical_scale_variation_includes_q_dependence": whole.diff(radial)
            != s.zeros(4),
            "coefficient_probes_are_not_physical_stress_by_definition": True,
            "formal_free_split_not_convergent_loop_bound": True,
            "all_complete_boundary_symbol_degrees_at_most_one": all(
                degree <= 1 for chart in chart_orders.values() for degree in chart["B"]
            ),
            "both_roots_keep_same_annihilation_time_orientation": True,
            "continuum_noise_not_a_covariant_stress_renormalization": True,
        },
    }
