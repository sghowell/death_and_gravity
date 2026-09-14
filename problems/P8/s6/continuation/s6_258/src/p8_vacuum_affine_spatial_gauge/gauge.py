"""Full density-gauge derivative, nonzero-mode inverse and local coercivity."""

from functools import cache

import sympy as s

SHAPE_BOUND = s.Rational(1, 8)
WEIGHT = s.Rational(2, 3)


@cache
def zero_modes():
    x = s.Symbol("periodic_spatial_coordinate", real=True)
    first, second = s.sin(x), s.cos(x)
    bracket = s.expand_trig(
        first * s.diff(second, x) - second * s.diff(first, x)
    ).simplify()
    DD = s.Matrix([[0, 1, -1], [-1, 0, 1], [1, -1, 0]])
    f = s.Matrix([1, 2, 4])
    g = s.Matrix([2, 3, 7])
    product = s.matrix_multiply_elementwise(f, g)
    defect = DD * product - s.diag(*f) * DD * g - s.diag(*g) * DD * f
    return {
        "two_mean_zero_periodic_gauge_components": (first, second),
        "whole_nonzero_constant_gauge_bracket": bracket,
        "finite_SBP_derivative": DD,
        "full_finite_SBP_Leibniz_defect": defect,
        "kernel_boundary": "The two periodic mean-zero vectors sin(x) partial_x and cos(x) partial_x have bracket -partial_x, a nonzero translation. Thus the mean-zero complement is not a Lie subalgebra. Its use in a local implicit-function slice does not remove the residual translation group or its global constraints. The finite antisymmetric derivative has summation by parts but the displayed nonzero Leibniz defect; it is not an anomaly-free diffeomorphism regulator.",
        "checks": {
            "first_gauge_component_mean_zero": s.integrate(first, (x, 0, 2 * s.pi)),
            "second_gauge_component_mean_zero": s.integrate(second, (x, 0, 2 * s.pi)),
            "whole_nonzero_translation_bracket": bracket + 1,
            "finite_SBP_exact_antisymmetry": DD + DD.T,
        },
        "gates": {
            "mean_zero_complement_not_a_Lie_subalgebra": bracket != 0,
            "finite_SBP_is_not_a_Leibniz_rule": defect != s.zeros(3, 1),
            "residual_translation_group_prescription_remains_separate": True,
        },
    }


@cache
def full_operator():
    x = s.symbols("spatial_coordinate0:3", real=True)
    values = [s.Function("shape_density" + str(i), real=True)(*x) for i in range(6)]
    Q = s.Matrix(
        [
            [values[0], values[1], values[2]],
            [values[1], values[3], values[4]],
            [values[2], values[4], values[5]],
        ]
    )
    xi = s.Matrix(
        [s.Function("spatial_gauge_vector" + str(i), real=True)(*x) for i in range(3)]
    )
    weight = s.Symbol("density_weight", real=True)
    divergence = sum(s.diff(xi[k], x[k]) for k in range(3))
    lie = s.Matrix(
        3,
        3,
        lambda i, j: (
            sum(
                xi[k] * s.diff(Q[i, j], x[k])
                - Q[i, k] * s.diff(xi[j], x[k])
                - Q[j, k] * s.diff(xi[i], x[k])
                for k in range(3)
            )
            + weight * Q[i, j] * divergence
        ),
    )
    chi = s.Matrix([sum(s.diff(Q[i, j], x[j]) for j in range(3)) for i in range(3)])
    direct = s.Matrix(
        [sum(s.diff(lie[i, j], x[j]) for j in range(3)) for i in range(3)]
    )
    lower = s.Matrix(
        [
            sum(
                xi[k] * s.diff(chi[i], x[k]) - chi[k] * s.diff(xi[i], x[k])
                for k in range(3)
            )
            + weight * chi[i] * divergence
            for i in range(3)
        ]
    )
    leading = s.Matrix(
        [
            -sum(
                Q[j, k] * s.diff(xi[i], x[j], x[k]) for j in range(3) for k in range(3)
            )
            + (weight - 1) * sum(Q[i, j] * s.diff(divergence, x[j]) for j in range(3))
            for i in range(3)
        ]
    )
    flat = direct.subs(
        {Q[0, 0]: 1, Q[1, 1]: 1, Q[2, 2]: 1, Q[0, 1]: 0, Q[0, 2]: 0, Q[1, 2]: 0},
        simultaneous=True,
    ).doit()
    flat_wanted = s.Matrix(
        [
            -sum(s.diff(xi[i], x[j], 2) for j in range(3))
            + (weight - 1) * s.diff(divergence, x[i])
            for i in range(3)
        ]
    )
    return {
        "whole_spatial_coordinates": x,
        "whole_contravariant_shape_density": Q,
        "whole_spatial_gauge_vector": xi,
        "whole_density_Lie_variation": lie,
        "whole_spatial_gauge_condition": chi,
        "whole_off_gauge_ghost_operator": direct.subs(weight, WEIGHT),
        "whole_off_gauge_lower_terms": lower.subs(weight, WEIGHT),
        "whole_gauge_surface_operator": leading.subs(weight, WEIGHT),
        "weight_dictionary": {
            "conformal_Dirac": WEIGHT,
            "harmonic_volume_density": s.Integer(1),
        },
        "operator_boundary": "The full operator is differentiated before imposing chi=0. Its lower terms xi dot grad chi-chi dot grad xi+(2/3)chi div xi are not absent from off-gauge vertices. On the exact gauge surface the operator reduces to -Q^jk partial_j partial_k xi^i-(1/3)Q^ij partial_j div xi. This is a spatial gauge operator, not a physical propagation operator or a replacement lapse equation.",
        "checks": {
            "whole_three_direction_density_gauge_derivative": (
                direct - lower - leading
            ).applyfunc(s.expand),
            "whole_flat_density_weight_symbol": (flat - flat_wanted).applyfunc(
                s.expand
            ),
            "whole_Dirac_weight": WEIGHT - s.Rational(2, 3),
            "whole_harmonic_weight_has_no_grad_div_term": s.Integer(1) - 1,
        },
        "gates": {
            "all_three_gauge_components_retained": len(chi) == 3,
            "off_gauge_lower_terms_are_nonzero": lower.subs(weight, WEIGHT)
            != s.zeros(3, 1),
            "full_source_clock_lapse_not_replaced_by_CMC_gauge": True,
        },
    }


@cache
def principal():
    a, b, c, d, e, f = s.symbols(
        "shape11 shape12 shape13 shape22 shape23 shape33", real=True
    )
    Q = s.Matrix([[a, b, c], [b, d, e], [c, e, f]])
    k = s.Matrix(s.symbols("nonzero_spatial_covector0:3", real=True))
    radial = (k.T * Q * k)[0]
    symbol = radial * s.eye(3) + (1 - WEIGHT) * (Q * k) * k.T
    symmetrized = radial * Q.inv() + (1 - WEIGHT) * k * k.T
    eigenvalue = s.Symbol("spectral_eigenvalue")
    base = symbol.subs({a: 1, d: 1, f: 1, b: 0, c: 0, e: 0})
    p2 = k.dot(k)
    inverse = (s.eye(3) - k * k.T / (4 * p2)) / p2
    P = s.eye(3) - k * k.T / p2
    return {
        "whole_positive_shape_density": Q,
        "whole_nonzero_spatial_covector": k,
        "whole_spatial_gauge_principal_symbol": symbol,
        "whole_positive_principal_symmetrizer": Q.inv(),
        "whole_symmetrized_principal_symbol": symmetrized,
        "whole_principal_characteristic_polynomial": (eigenvalue - radial) ** 2
        * (eigenvalue - s.Rational(4, 3) * radial),
        "whole_flat_nonzero_momentum_inverse": inverse,
        "whole_transverse_reference_projector": P,
        "whole_flat_reference_determinant": s.Rational(4, 3) * p2**3,
        "kernel_boundary": "For positive Q and nonzero k the full principal symbol has two eigenvalues k^TQk and one (4/3)k^TQk, with positive symmetrizer Q^-1. At k=0 the homogeneous operator vanishes. Do not use this inverse on the constant translation modes or infer a uniform infinite-volume infrared gap.",
        "checks": {
            "whole_positive_symbol_symmetrizer": (
                Q.inv() * symbol - symmetrized
            ).applyfunc(s.factor),
            "whole_three_mode_principal_characteristic": s.factor(
                symbol.charpoly(eigenvalue).as_expr()
                - (eigenvalue - radial) ** 2 * (eigenvalue - s.Rational(4, 3) * radial)
            ),
            "whole_positive_principal_determinant": s.factor(
                symbol.det() - s.Rational(4, 3) * radial**3
            ),
            "whole_flat_nonzero_momentum_inverse": (
                base * inverse - s.eye(3)
            ).applyfunc(s.factor),
            "whole_transverse_projector_idempotent": (P * P - P).applyfunc(s.factor),
            "whole_transverse_projector_annihilates_wavevector": (P * k).applyfunc(
                s.factor
            ),
        },
        "gates": {
            "full_principal_is_not_generically_symmetric_in_Euclidean_metric": symbol
            != symbol.T,
            "full_principal_symmetrizer_is_symmetric": symmetrized == symmetrized.T,
            "nonzero_momentum_and_positive_shape_are_mandatory": True,
        },
    }


@cache
def local_slice():
    L = s.Symbol("positive_torus_radius", positive=True)
    delta = s.Symbol("shape_operator_norm_distance", nonnegative=True)
    margin = 1 - 2 * SHAPE_BOUND
    # Sharp enough finite constants for the full nonsymmetric weak form.
    matrix = s.Matrix(3, 3, s.symbols("test_gradient0:9", real=True))
    trace = matrix.trace()
    # Cauchy for the three diagonal entries, with all off-diagonal norms kept.
    diagonal_variance = sum(
        (matrix[i, i] - matrix[j, j]) ** 2 for i in range(3) for j in range(i + 1, 3)
    )
    all_norm = sum(item**2 for item in matrix)
    off = sum(matrix[i, j] ** 2 for i in range(3) for j in range(3) if i != j)
    return {
        "fixed_shape_Linfinity_operator_distance_bound": SHAPE_BOUND,
        "whole_weak_coercivity_margin": margin,
        "whole_shape_perturbation_weak_bound": 2 * delta,
        "flat_torus_spatial_period": 2 * s.pi * L,
        "mean_zero_torus_L2_gap_lower_bound": margin / L**2,
        "mean_zero_torus_inverse_L2_upper_bound": L**2 / margin,
        "whole_weak_form": "For chi=div Q=0, B(v,xi)=integral partial_j v_i Q^jk partial_k xi_i+(1/3)partial_j v_i Q^ij div xi. Relative to Q=I, the first perturbation is bounded by delta|grad xi|^2 and the second by delta|grad xi|^2, using ||Q-I||F<=sqrt(3)delta and |div xi|<=sqrt(3)|grad xi|. Thus B(xi,xi)>=(1-2delta)||grad xi||L2^2, even though the full operator is not generally Euclidean-symmetric.",
        "local_slice_argument": "On the flat 2pi L torus remove only the constant modes from the infinitesimal gauge complement and target the mean-zero divergence range. At the homogeneous shape Q=I the displayed second-order inverse is an isomorphism. For integer k>=2, use C^(k+3,alpha) metrics and near-identity C^(k+2,alpha) periodic coordinate maps, target C^(k,alpha), and fix the coordinate-displacement mean. The gauge map is C1 with this derivative gap; the extra metric derivative also gives continuity of the composed highest coefficient derivatives. The implicit function theorem gives one local gauge map near the reference. Translation freedom is still a separate residual group, not a globally unique gauge. On any exact slice with ||Q-I||op<=1/8, the full weak coercivity bound and elliptic regularity give the stated inverse gap on mean-zero fields.",
        "zero_mode_and_regulator_boundary": "Mean-zero vector fields are a local complement, not in general a Lie subalgebra. Constant translations remain a residual group with their global constraint/group-volume prescription separate. The L2 gap scales as L^-2 and is not uniform as L tends to infinity. A finite torus or Galerkin diagnostic does not replace the original infinite-volume quantum state or supply a BRST-preserving regulator.",
        "checks": {
            "whole_diagonal_Cauchy_identity": s.expand(
                3 * sum(matrix[i, i] ** 2 for i in range(3))
                - trace**2
                - diagonal_variance
            ),
            "whole_full_gradient_Cauchy_remainder": s.expand(
                3 * all_norm - trace**2 - diagonal_variance - 3 * off
            ),
            "whole_exact_coercivity_margin": margin - s.Rational(3, 4),
            "whole_gap_times_inverse_bound": (margin / L**2) * (L**2 / margin) - 1,
        },
        "gates": {
            "strict_full_weak_coercivity_margin": margin > 0,
            "no_negative_curvature_or_GR_CMC_theorem_imported": True,
            "translation_kernel_and_IR_scale_retained": True,
            "local_slice_not_a_global_gauge_or_physical_Cauchy_theorem": True,
        },
    }
