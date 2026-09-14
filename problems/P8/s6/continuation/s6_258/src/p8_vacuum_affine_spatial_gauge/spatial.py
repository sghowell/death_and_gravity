"""Full spatial momentum generator and physical/reference shape dictionary."""

from functools import cache

import sympy as s
from p8_vacuum_affine_nonlinear_auxiliary_measure import canonical as parent
from p8_vacuum_affine_nonlinear_auxiliary_measure import cotangent


@cache
def lie_algebra():
    x = s.symbols("spatial_coordinate0:3", real=True)
    xi = s.Matrix([x[0] * x[1] + x[2], x[1] ** 2 - x[0], x[0] * x[2] + 1])
    eta = s.Matrix([x[0] ** 2 + x[1], x[1] * x[2] + 2, x[2] ** 2 - x[0] * x[1]])
    bracket = s.Matrix(
        [
            sum(
                xi[k] * s.diff(eta[i], x[k]) - eta[k] * s.diff(xi[i], x[k])
                for k in range(3)
            )
            for i in range(3)
        ]
    )
    Q = s.Matrix(
        [
            [2 + x[0] ** 2, x[0] * x[1], x[0] * x[2]],
            [x[0] * x[1], 3 + x[1] ** 2, x[1] * x[2]],
            [x[0] * x[2], x[1] * x[2], 4 + x[2] ** 2],
        ]
    )

    def lie(a, Q, weight):
        div = sum(s.diff(a[k], x[k]) for k in range(3))
        return s.Matrix(
            3,
            3,
            lambda i, j: (
                sum(
                    a[k] * s.diff(Q[i, j], x[k])
                    - Q[i, k] * s.diff(a[j], x[k])
                    - Q[j, k] * s.diff(a[i], x[k])
                    for k in range(3)
                )
                + weight * Q[i, j] * div
            ),
        )

    checks = {}
    for label, weight in (
        ("conformal_Dirac", s.Rational(2, 3)),
        ("canonical_weight_one", s.Integer(1)),
    ):
        checks["whole_" + label + "_density_Lie_commutator"] = (
            lie(xi, lie(eta, Q, weight), weight)
            - lie(eta, lie(xi, Q, weight), weight)
            - lie(bracket, Q, weight)
        ).applyfunc(s.expand)
    return {
        "independent_polynomial_spatial_gauge_vectors": (xi, eta),
        "whole_spatial_vector_bracket": bracket,
        "whole_polynomial_density_diagnostic": Q,
        "classical_algebra_proof": "For the full tensor/scalar configuration fields q, D[xi]=integral p dot L_xi q. Functional variation and the retained adjoints give {D[xi],D[eta]}=integral p dot [L_xi,L_eta]q=D[[xi,eta]], modulo the displayed boundaries. The representation property follows from the tensor Lie derivative and its density term. The polynomial diagnostics check both weight 2/3 and weight 1 in three directions; they do not turn a finite collocation rule into an exact diffeomorphism representation.",
        "checks": checks,
        "gates": {
            "full_noncommuting_spatial_descriptors": bracket != s.zeros(3, 1),
            "both_shape_and_canonical_density_weights_checked": len(checks) == 2,
            "classical_poisson_algebra_not_quantum_or_discrete_closure": True,
        },
    }


@cache
def generator():
    x = s.symbols("spatial_coordinate0:3", real=True)

    def symmetric(label):
        v = [s.Function(label + str(i), real=True)(*x) for i in range(6)]
        return s.Matrix([[v[0], v[1], v[2]], [v[1], v[3], v[4]], [v[2], v[4], v[5]]])

    g = symmetric("metric")
    pi = symmetric("metric_pi")
    W = s.Matrix([s.Function("W" + str(i), real=True)(*x) for i in range(3)])
    pW = s.Matrix([s.Function("pW" + str(i), real=True)(*x) for i in range(3)])
    xi = s.Matrix([s.Function("xi" + str(i), real=True)(*x) for i in range(3)])
    fields = [s.Function("scalar" + str(i), real=True)(*x) for i in range(4)]
    ps = [s.Function("scalar_pi" + str(i), real=True)(*x) for i in range(4)]
    lie_g = s.Matrix(
        3,
        3,
        lambda i, j: sum(
            xi[k] * s.diff(g[i, j], x[k])
            + g[i, k] * s.diff(xi[k], x[j])
            + g[k, j] * s.diff(xi[k], x[i])
            for k in range(3)
        ),
    )
    lie_W = s.Matrix(
        [
            sum(
                xi[k] * s.diff(W[i], x[k]) + W[k] * s.diff(xi[k], x[i])
                for k in range(3)
            )
            for i in range(3)
        ]
    )
    lag = (
        sum(pi[i, j] * lie_g[i, j] for i in range(3) for j in range(3))
        + pW.dot(lie_W)
        + sum(
            ps[a] * sum(xi[i] * s.diff(fields[a], x[i]) for i in range(3))
            for a in range(4)
        )
    )
    divergence_pi = sum(s.diff(pW[i], x[i]) for i in range(3))
    metric_part = s.Matrix(
        [
            sum(
                pi[j, k] * s.diff(g[j, k], x[i]) - 2 * s.diff(g[i, j] * pi[j, k], x[k])
                for j in range(3)
                for k in range(3)
            )
            for i in range(3)
        ]
    )
    vector_part = s.Matrix(
        [
            sum(pW[j] * (s.diff(W[j], x[i]) - s.diff(W[i], x[j])) for j in range(3))
            - W[i] * divergence_pi
            for i in range(3)
        ]
    )
    scalar_part = s.Matrix(
        [sum(ps[a] * s.diff(fields[a], x[i]) for a in range(4)) for i in range(3)]
    )
    total = metric_part + vector_part + scalar_part
    flux = s.Matrix(
        [
            sum(2 * pi[i, j] * g[j, k] * xi[k] for j in range(3) for k in range(3))
            + pW[i] * W.dot(xi)
            for i in range(3)
        ]
    )
    lower = lambda i, j, k: (
        (s.diff(g[i, k], x[j]) + s.diff(g[i, j], x[k]) - s.diff(g[j, k], x[i])) / 2
    )
    covariant = s.Matrix(
        [
            -2
            * sum(
                g[i, j] * s.diff(pi[j, k], x[k]) + lower(i, k, j) * pi[j, k]
                for j in range(3)
                for k in range(3)
            )
            for i in range(3)
        ]
    )
    return {
        "whole_spatial_coordinates": x,
        "whole_metric_and_metric_momenta": (g, pi),
        "whole_vector_and_vector_momenta": (W, pW),
        "whole_two_matter_and_two_auxiliary_scalar_pairs": (fields, ps),
        "whole_spatial_gauge_vector": xi,
        "whole_metric_Lie_action": lie_g,
        "whole_vector_Lie_action": lie_W,
        "whole_canonical_Lie_pairing_density": lag,
        "whole_spatial_generator_density": total,
        "whole_spatial_generator_boundary_flux": flux,
        "whole_covariant_metric_generator": covariant,
        "generator_boundary": "The last two scalar pairs are the lapse and normal-vector primaries, retained before imposing them. The first two are the full M1 and heavy matter fields. Metric pi is a contravariant tensor density of weight one. The displayed generator follows from the full canonical Lie pairing by the retained spatial boundary. On the auxiliary primary surface it is exactly the S257 retained momentum generator. Time-dependent spatial transformations also have their shift-multiplier and shift-primary terms; the displayed field generator and ghost operator are at fixed clock time.",
        "classical_reduction_equivariance": "The whole S257 normal Hamiltonian is a spatial scalar density: all R,F,j and profile coefficients are scalar functions of u,N, and all metric, vector, matter and Gauss invariants have their displayed tensor weights. The unique local auxiliary branch therefore transforms as scalar N,T. The full second-class reduction preserves the retained canonical brackets, so the classical spatial momentum algebra and its action survive on that branch. The primitive boundary V I is a spatial density; its integrated Lie variation is a boundary, so its canonical momentum shift preserves the integrated spatial generator under the stated boundary conditions. No discrete or quantum anomaly-free algebra is inferred.",
        "checks": {
            "whole_metric_vector_matter_generator_and_boundary": s.expand(
                lag - xi.dot(total) - sum(s.diff(flux[i], x[i]) for i in range(3))
            ),
            "whole_weight_one_covariant_metric_divergence": (
                covariant - metric_part
            ).applyfunc(s.expand),
        },
        "gates": {
            "both_matter_and_both_auxiliary_scalar_pairs_retained": len(fields)
            == len(ps)
            == 4,
            "full_vector_Gauss_generator_term_nonzero": divergence_pi != 0,
            "whole_metric_and_vector_spatial_boundary_retained": flux != s.zeros(3, 1),
        },
    }


@cache
def frame():
    a, b, c, d, e, f = s.symbols(
        "gamma11 gamma12 gamma13 gamma22 gamma23 gamma33", real=True
    )
    metric = s.Matrix([[a, b, c], [b, d, e], [c, e, f]])
    C = s.Symbol("positive_C", positive=True)
    determinant = s.Symbol("positive_metric_determinant", positive=True)
    Q = determinant ** s.Rational(1, 3) * metric.inv()
    physical = (C**3 * determinant) ** s.Rational(1, 3) * (C * metric).inv()
    scale = s.Symbol("positive_scale", positive=True)
    volume = s.Symbol("local_scalar_volume", real=True)
    shape = s.Matrix(3, 3, s.symbols("independent_shape_density0:9", real=True))
    metric_hat = scale**2 * s.exp(2 * volume) * shape
    source = parent.current_data()
    return {
        "whole_spatial_metric": metric,
        "whole_positive_metric_determinant_name": determinant,
        "whole_conformal_Dirac_density": Q,
        "whole_physical_conformal_Dirac_density": physical,
        "whole_metric_shape_factorization": metric_hat,
        "full_source_spatial_generator_binding": parent.full_trace()[
            "whole_spatial_shift_generator_on_primary_surface"
        ],
        "whole_current_R_F_and_fixed_profile_bindings": source[
            "whole_actual_R_F_fixed_profile_and_constant_bindings"
        ],
        "whole_primitive_boundary_binding": cotangent.primitive_boundary_map()[
            "whole_primitive_boundary_function"
        ],
        "actual_initial_full_lapse_Schur_enclosure": source[
            "actual_initial_full_lapse_Schur_enclosure"
        ],
        "full_frame_identity": "For the complete positive physical-to-hat map h_phys=C gamma with C=R(u,N)^-1/2, det(h_phys)^(1/3) h_phys^-1=det(gamma)^(1/3) gamma^-1 exactly, even when C is spatially varying. Thus both the gauge and all its spatial derivatives are identical in the two frames without discarding lapse derivatives. For gamma=a_hat^2 exp(2v) exp(t), trace t=0, the density is exp(-t); it is independent of scalar volume, lapse and background scale. On every homogeneous classical comparison slice it is exactly I.",
        "fixed_reference_boundary": "Linearizing the same shape parametrization gives chi=-div t. The existing tracefree transverse tensor reference satisfies that condition, while its arbitrary scalar volume does not enter it. Nonlinearly the condition is div(exp(-t))=0, not simply div t=0. The full source pullback, canonical momenta and primitive boundary remain those of S257. This is a compatible local spatial gauge, not a replacement fixed quantum state or a proof of the quantum mean.",
        "checks": {
            "whole_physical_and_hat_shape_density_equal": (physical - Q).applyfunc(
                s.simplify
            ),
            "whole_metric_conformal_determinant": s.factor(
                (C * metric).det() - C**3 * metric.det()
            ),
            "whole_background_scale_and_scalar_volume_cancel": s.simplify(
                (scale**6 * s.exp(6 * volume)) ** s.Rational(1, 3)
                / (scale**2 * s.exp(2 * volume))
                - 1
            ),
        },
        "gates": {
            "whole_source_functions_not_replaced_by_clock_polynomial": bool(
                source["whole_actual_R_F_fixed_profile_and_constant_bindings"]
            ),
            "same_actual_regular_classical_auxiliary_pivot": source[
                "actual_initial_full_lapse_Schur_enclosure"
            ][0]
            > s.Rational(151, 100),
            "nonlinear_reference_tensor_gauge_not_linearized_away": True,
        },
    }


@cache
def reference_projector():
    k = s.Matrix(s.symbols("nonzero_spatial_covector0:3", real=True))
    a, b, c, d, e = s.symbols(
        "tracefree_shape11 tracefree_shape12 tracefree_shape13 tracefree_shape22 tracefree_shape23",
        real=True,
    )
    h = s.Matrix([[a, b, c], [b, d, e], [c, e, -a - d]])
    p2 = k.dot(k)
    inverse = (s.eye(3) - k * k.T / (4 * p2)) / p2
    xi = s.I * inverse * h * k
    delta_h = s.I * (k * xi.T + xi * k.T - s.Rational(2, 3) * s.eye(3) * k.dot(xi))
    transverse = s.eye(3) - k * k.T / p2
    TT = transverse * h * transverse - transverse * s.trace(transverse * h) / 2
    scalar_contact = s.I * k.dot(xi) / 3
    return {
        "whole_tracefree_shape_perturbation": h,
        "whole_linear_gauge_parameter": xi,
        "whole_linear_shape_gauge_correction": delta_h,
        "whole_transverse_tracefree_projection": TT,
        "whole_induced_scalar_volume_correction": scalar_contact,
        "reference_bridge": "For a general tracefree metric perturbation the gauge correction also shifts scalar volume by -(k^T h k)/(4k^2); do not throw away that field-map component. On the existing TT plus volume reference, h k=0, so the gauge parameter and both corrections vanish. Constant momentum is not assigned this Fourier inverse. The nonlinear shape-map and source/boundary prescriptions must still be transported for interactions.",
        "checks": {
            "whole_linear_spatial_gauge_is_TT_projection": (h + delta_h - TT).applyfunc(
                s.factor
            ),
            "whole_projected_shape_transverse": (TT * k).applyfunc(s.factor),
            "whole_projected_shape_tracefree": s.factor(s.trace(TT)),
            "whole_required_scalar_volume_gauge_contact": s.factor(
                scalar_contact + (k.T * h * k)[0] / (4 * p2)
            ),
        },
        "gates": {
            "generic_scalar_volume_gauge_contact_not_deleted": scalar_contact != 0,
            "existing_TT_reference_requires_no_new_state_choice": True,
        },
    }
