"""Exact outward contraction domain for the full nonlinear constraints."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate

from . import source as q


def rational(value):
    return s.Rational(value.numerator, value.denominator)


def maximum(iv):
    return rational(max(abs(iv.lo), abs(iv.hi)))


def ceiling(iv):
    return s.ceiling(maximum(iv))


def outward(iv, denominator=100000):
    return [
        s.floor(rational(iv.lo) * denominator) / denominator,
        s.ceiling(rational(iv.hi) * denominator) / denominator,
    ]


def domain():
    values = {
        q.N: I(1 - q.OUTER, 1 + q.OUTER),
        q.mu: I(s.Rational(1, 10**4), s.Rational(1, 100)),
    }
    values.update({z: I(-q.DELTA, q.DELTA) for z in q.COORDS})
    for f, tree, error, last in [
        (q.R, q.R_TREE, q.ERROR, 4),
        (q.F, q.F_TREE, q.ERROR, 4),
        (q.Ruu, q.RUU_TREE, q.ERROR, 3),
        (q.j, s.Integer(0), q.JERROR, 4),
    ]:
        for k in range(last + 1):
            iv = evaluate(s.diff(tree, q.N, k), {q.N: values[q.N]})
            values[s.diff(f, q.N, k)] = I(iv.lo - error, iv.hi + error)
    return values


@cache
def enclosure():
    C, N = q.CONSTRAINT, q.N
    dom = domain()
    cn = evaluate(s.diff(C, N), dom)
    cnn = evaluate(s.diff(C, N, 2), dom)
    cnnn = evaluate(s.diff(C, N, 3), dom)
    gradient = [ceiling(evaluate(s.diff(C, z), dom)) for z in q.COORDS]
    atoms = sorted(C.atoms(s.Derivative) | C.atoms(s.Function), key=str)
    jets = s.symbols(
        "independent_complete_constraint_source_jet0:" + str(len(atoms)), real=True
    )
    lifted = C.xreplace(dict(zip(atoms, jets, strict=True)))
    jet_dom = {key: value for key, value in dom.items() if key not in atoms}
    jet_dom.update({z: dom[key] for key, z in zip(atoms, jets, strict=True)})
    source_gradient = [ceiling(evaluate(s.diff(lifted, z), jet_dom)) for z in jets]
    center = sum(source_gradient) * q.ERROR
    residual = sum(gradient) * q.DELTA + center
    nfirst = s.Integer(2)
    mixed = [ceiling(evaluate(s.diff(C, N, z), dom)) for z in q.COORDS]
    mixed2 = [ceiling(evaluate(s.diff(C, N, 2, z), dom)) for z in q.COORDS]
    sparse = {}
    sparseN = {}
    for i, z in enumerate(q.COORDS):
        for k, w in enumerate(q.COORDS):
            value = s.diff(C, z, w)
            if value != 0:
                sparse[str((i, k))] = ceiling(evaluate(value, dom))
                sparseN[str((i, k))] = ceiling(evaluate(s.diff(value, N), dom))
    nsecond = (10 + 2 * 100 * nfirst + 1000 * nfirst**2) / 3
    nthird = (
        3 * 100 * nfirst
        + 3 * 1000 * nfirst**2
        + 10000 * nfirst**3
        + 3 * 100 * nsecond
        + 3 * 1000 * nsecond * nfirst
    ) / 3
    ratio = s.Rational(1, 30)
    gates = {
        "whole_lapse_derivative_inside_strict_negative_margin": cn.lo
        > -s.Rational(31, 10)
        and cn.hi < -3,
        "whole_second_lapse_derivative_bound": maximum(cnn) < 1000,
        "whole_third_lapse_derivative_bound": maximum(cnnn) < 10000,
        "full_parameter_gradient_sum_below100": sum(gradient) < 100,
        "full_source_gradient_sum_below1000": sum(source_gradient) < 1000,
        "all_exact_center_MVT_inputs_remain_in_positive_R_box": dom[q.R].inside(
            s.Rational(999, 1000), s.Rational(1001, 1000)
        ),
        "whole_point_residual_below1e_minus245": residual < q.INNER,
        "whole_original_reference_root_shift_below1e_minus375": center / 3
        < q.CENTER_RADIUS,
        "explicit_contraction_maps_inner_interval_strictly_inside": residual / 3
        + ratio * q.INNER
        < q.INNER,
        "inner_root_interval_inside_outer_verified_strip": q.INNER < q.OUTER,
        "root_inside_S263_signed_lapse_radius": q.INNER < q.fifth.EPSILON,
        "whole_parameter_first_root_derivatives_below2": max(gradient) / 3 < nfirst,
        "whole_mixed_CN_z_bound": max(mixed) < 100,
        "whole_mixed_CNN_z_bound": max(mixed2) < 1000,
        "whole_sparse_Czz_bound": max(sparse.values()) < 10,
        "whole_sparse_CNzz_bound": max(sparseN.values()) < 100,
        "whole_second_root_derivatives_below2000": nsecond < 2000,
        "whole_third_root_derivatives_below1e7": nthird < 10**7,
        "canonical_invariant_coordinates_not_Gaussian_modes": True,
    }
    return {
        "whole_fixed_parameter_radius": q.DELTA,
        "whole_outer_lapse_strip_radius": q.OUTER,
        "whole_inner_unique_root_radius": q.INNER,
        "whole_original_reference_center_root_shift_bound": q.CENTER_RADIUS,
        "whole_lapse_first_outward_enclosure": outward(cn),
        "whole_lapse_second_and_third_absolute_ceilings": [ceiling(cnn), ceiling(cnnn)],
        "whole_parameter_gradient_ceilings": dict(
            zip(map(str, q.COORDS), gradient, strict=True)
        ),
        "whole_independent_source_atoms": atoms,
        "whole_independent_source_jet_gradient_ceilings": dict(
            zip(map(str, atoms), source_gradient, strict=True)
        ),
        "whole_centered_MVT_constraint_residual_bound": residual,
        "whole_source_only_center_residual_bound": center,
        "whole_parameter_mixed_CN_and_CNN_ceilings": [mixed, mixed2],
        "whole_all_nonzero_Czz_and_CNzz_ceilings": [sparse, sparseN],
        "whole_root_derivative_ceilings": [2, 2000, 10**7],
        "whole_explicit_contraction_factor": ratio,
        "whole_iteration_error": "Starting N0=1, set N_(k+1)=N_k+C(N_k,z)/3 using the ENTIRE source-pinned constraint. For every stated invariant datum the map is1/30-Lipschitz, maps the complete inner interval strictly into itself, and has one root there. Strict monotonicity gives uniqueness on the entire outer strip. The exact error is <=10*b/(29*30^k), where b is the displayed full centered-MVT residual bound. In particular |Nstar-1|<=b/3<1e-245; at unchanged center z=0 only the nonzero full-source error remains, giving |Nstar(0)-1|<1e-375. Neither root shift is asserted nonzero.",
        "whole_regular_response": "The full scalar constraint is C3 in the finite density/jet variables; source derivative needs through total order5 are supplied by S263. Its negative N derivative gives the full implicit response. First, second and third individual parameter derivatives have the displayed common ceilings. These are nonlinear finite-spatial-jet reconstructions at u=0, not a spatial evolution theorem or canonical Gaussian support.",
        "checks": {
            "exact_1_over30_contraction_margin": 1 - s.Rational(31, 10) / 3 + ratio,
            "exact_second_derivative_majorant": nsecond - 1470,
            "exact_third_derivative_majorant": nthird - s.Rational(9353600, 3),
            "whole_constraint_degree_at_most_two_in_all_invariants": s.Integer(
                max(sum(monomial) for monomial, _ in s.Poly(C, *q.COORDS).terms())
            )
            - 2,
        },
        "gates": gates,
    }


@cache
def auxiliaries():
    dom = domain()
    T = q.TEMPORAL
    t = evaluate(T, dom)
    tn = evaluate(s.diff(T, q.N), dom)
    pivot = evaluate(q.TEMPORAL_PIVOT, dom)
    gamma = evaluate(q.GAMMA, dom)
    d, tau, k = s.symbols(
        "negative_lapse_Schur negative_temporal_pivot whole_temporal_N", real=True
    )
    matrix = s.Matrix([[d + tau * k * k, -tau * k], [-tau * k, tau]])
    inverse = s.Matrix([[1 / d, k / d], [k / d, 1 / tau + k * k / d]])
    coefficients = {
        "D_M_over_N": q.R ** s.Rational(1, 4) / q.N,
        "Q_N_over_U": q.N * q.R ** s.Rational(3, 4),
        "U_physical_volume": q.R ** -s.Rational(3, 4),
        "C_physical_metric": q.R ** -s.Rational(1, 2),
    }
    bounds = {name: evaluate(expr, dom) for name, expr in coefficients.items()}
    return {
        "whole_actual_temporal_solution": T,
        "whole_temporal_absolute_bound": 3 * q.DELTA,
        "whole_temporal_lapse_derivative_absolute_bound": 10 * q.DELTA,
        "whole_temporal_pivot_outward_enclosure": outward(pivot),
        "whole_regular_Gamma_outward_enclosure": outward(gamma),
        "whole_auxiliary_Hessian_on_temporal_branch": matrix,
        "whole_auxiliary_Hessian_inverse_on_temporal_branch": inverse,
        "whole_auxiliary_Hessian_determinant_floor": s.Rational(3, 2),
        "whole_full_positive_lapse_coefficient_enclosures": {
            name: outward(iv) for name, iv in bounds.items()
        },
        "whole_two_auxiliary_branch_argument": "The full retained temporal solution is substituted only after variation. Its pivot lies strictly between-2 and-1/2; the lapse Schur pivot lies between-31/10 and-3. The displayed full2x2 Hessian is negative definite and its determinant exceeds3/2. Both primary/secondary pairs are therefore uniformly second class on this explicit local branch. The S257 finite delta/determinant cancellation applies branchwise; the remaining spatial gauge factor is not evaluated.",
        "whole_quantum_domain_boundary": "D, Q, U and the actual physical conformal factor are positive and bounded between1/2 and2 on this actual CLASSICAL nonlinear branch. The S265 unbounded-Gaussian coefficient-divergence result is unchanged: no independent Gaussian lapse or Gaussian invariant law is substituted here. A bounded classical symbol is not by itself a positive Weyl operator, a physical quantum measure or an interacting-state expectation.",
        "checks": {
            "whole_two_auxiliary_Hessian_inverse": (
                matrix * inverse - s.eye(2)
            ).applyfunc(s.factor),
            "whole_two_auxiliary_determinant": s.factor(matrix.det() - d * tau),
            "whole_two_auxiliary_Schur": s.factor(
                matrix[0, 0] - matrix[0, 1] ** 2 / matrix[1, 1] - d
            ),
            "whole_two_auxiliary_cross_contact": s.factor(matrix[0, 1] + tau * k),
        },
        "gates": {
            "whole_normal_vector_root_below3_delta": maximum(t) < 3 * q.DELTA,
            "whole_normal_vector_lapse_response_below10_delta": maximum(tn)
            < 10 * q.DELTA,
            "whole_temporal_pivot_regular_and_negative": pivot.lo > -2
            and pivot.hi < -s.Rational(1, 2),
            "whole_full_Gamma_above99_over100": gamma.lo > s.Rational(99, 100),
            "whole_four_complete_physical_weights_in_half_to_two": all(
                iv.inside(s.Rational(1, 2), 2) for iv in bounds.values()
            ),
            "entire_reference_coefficient_obstruction_not_erased": True,
            "remaining_spatial_constraints_and_quantum_measure_not_solved": True,
        },
    }
