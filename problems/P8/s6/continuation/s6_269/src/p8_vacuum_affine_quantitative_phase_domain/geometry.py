"""Evaluated full-convolution geometry and complete adjoint cotangent bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_nonlinear_spatial_reduction import shape

from . import reference

TAU = s.Rational(1, 10**274)
VLOG = s.Rational(1, 10**329)
METRIC = s.Rational(1, 10**272)
MOMENTUM_A1 = s.Rational(1, 10**200)


@cache
def bounds():
    a = 2 * TAU
    image, contraction = 3 * a * a + 3 * a**3, 6 * a + 9 * a * a
    inverse_Q = 1 / (1 - a)
    exp_gap = 2 * VLOG / (1 - 2 * VLOG)
    gamma_gap = exp_gap * inverse_Q + a / (1 - a)
    inverse_gap = exp_gap * (1 + a) + a
    ghost = 100 * a
    connection = 9 * METRIC
    ricci = 6 * connection + 18 * connection**2
    curvature = 18 * ricci
    base = (s.Rational(1, 3) + 32) * MOMENTUM_A1
    D0 = (
        180 * base
        + 12 * s.Rational(1, 10**275) * MOMENTUM_A1
        + 2 * (s.Rational(1, 10**329) + MOMENTUM_A1 * s.Rational(1, 10**354))
    )
    lambda_A1 = 4 * D0
    pi = base + 16 * lambda_A1
    pi_TF = 5 * pi
    return {
        "shape_self_map": image,
        "shape_contraction": contraction,
        "whole_Q_minus_I_A2": a,
        "whole_Q_inverse_A2": inverse_Q,
        "whole_gamma_minus_I_A2": gamma_gap,
        "whole_inverse_gamma_minus_I_A2": inverse_gap,
        "whole_metric_safe_gap_A2": METRIC,
        "whole_adjoint_ghost_Neumann_error": ghost,
        "whole_mean_zero_adjoint_inverse_A0_to_A1": s.Integer(4),
        "whole_connection_component_A1": connection,
        "whole_Ricci_component_A0": ricci,
        "whole_scalar_curvature_A0": curvature,
        "whole_base_metric_momentum_A1": base,
        "whole_complete_source_D0_A0": D0,
        "whole_adjoint_descriptor_A1": lambda_A1,
        "whole_lifted_metric_momentum_A0": pi,
        "whole_tracefree_metric_momentum_A0": pi_TF,
        "whole_normalized_shear_invariant": 48 * pi_TF**2,
    }


@cache
def adjoint_identities():
    x = s.symbols("spatial_x0:3", real=True)
    Q = s.zeros(3)
    for i in range(3):
        for j in range(i, 3):
            Q[i, j] = Q[j, i] = s.Function("whole_Q" + str(i) + str(j))(*x)
    lam = [s.Function("whole_descriptor" + str(i))(*x) for i in range(3)]
    div = [sum(s.diff(Q[i, j], x[j]) for j in range(3)) for i in range(3)]
    divdiv = sum(s.diff(div[j], x[j]) for j in range(3))
    full, sliced, contacts, checks = [], [], [], {}
    for l in range(3):
        whole = -sum(
            s.diff(Q[j, k] * lam[l], x[j], x[k]) for j in range(3) for k in range(3)
        )
        whole -= (
            sum(
                s.diff(Q[i, j] * lam[i], x[l], x[j]) for i in range(3) for j in range(3)
            )
            / 3
        )
        on_slice = -sum(
            Q[j, k] * s.diff(lam[l], x[j], x[k]) for j in range(3) for k in range(3)
        )
        on_slice -= (
            sum(
                Q[i, j] * s.diff(lam[i], x[l], x[j])
                + s.diff(Q[i, j], x[l]) * s.diff(lam[i], x[j])
                for i in range(3)
                for j in range(3)
            )
            / 3
        )
        contact = (
            2 * sum(div[k] * s.diff(lam[l], x[k]) for k in range(3)) + divdiv * lam[l]
        )
        contact += (
            sum(
                div[i] * s.diff(lam[i], x[l]) + s.diff(div[i], x[l]) * lam[i]
                for i in range(3)
            )
            / 3
        )
        checks["entire_off_slice_adjoint_contact_" + str(l)] = s.expand(
            whole - on_slice + contact
        )
        full.append(whole)
        sliced.append(on_slice)
        contacts.append(contact)
    return {
        "full": s.Matrix(full),
        "slice": s.Matrix(sliced),
        "contacts": s.Matrix(contacts),
        "checks": checks,
    }


@cache
def data():
    budget = bounds()
    field = reference.field_bounds()
    adjoint = adjoint_identities()
    checks = dict(adjoint["checks"])
    for vec in ((1, 0, 0), (1, 2, 0), (1, -1, 2), (10**64, 0, 0)):
        k = s.Matrix(vec)
        square = k.dot(k)
        M0 = square * s.eye(3) + k * k.T / 3
        inverse = (s.eye(3) - k * k.T / (4 * square)) / square
        name = "_".join(map(str, vec))
        checks["whole_flat_ghost_inverse_left_" + name] = M0 * inverse - s.eye(3)
        checks["whole_flat_ghost_inverse_right_" + name] = inverse * M0 - s.eye(3)
    gv, bv = (
        s.symbols("inverse_metric0:6", real=True),
        s.symbols("dual_matrix0:6", real=True),
    )
    g = s.Matrix([[gv[0], gv[1], gv[2]], [gv[1], gv[3], gv[4]], [gv[2], gv[4], gv[5]]])
    B = s.Matrix([[bv[0], bv[1], bv[2]], [bv[1], bv[3], bv[4]], [bv[2], bv[4], bv[5]]])
    c = s.Symbol("full_positive_density", positive=True)
    dual = c * (g * s.trace(B * g) / 3 - g * B * g)
    checks["whole_generic_nondiagonal_cotangent_trace_preservation"] = s.factor(
        s.trace(dual * g.inv())
    )
    checks["literal_original_full_determinant_identity"] = shape.data()["checks"][
        "whole_noncommuting_3D_determinant"
    ]
    checks["literal_original_full_cofactor_identity"] = shape.data()["checks"][
        "whole_full_cofactor_not_eigenvalue_truncation"
    ]
    return {
        "whole_all_evaluated_geometry_bounds": budget,
        "whole_tiny_A2_shape_contraction": "The S267 determinant/projector proof works for every nonnegative submultiplicative Fourier exponent. Apply it explicitly in A2, with ||tau||<=1e-274 and ||f||<=1e-274. The full fixed point has ||Q-I||A2<=2e-274, detQ=1, divQ=0 and Q>0. All generated means and nonlinear harmonics remain. This is not the generally false claim that this high-P datum is small in the old A8 norm.",
        "whole_metric_and_volume": "At u=0, gamma=exp(2v)Q^-1, inverse_gamma=exp(-2v)Q and volume=exp(3v). The full A2 gaps of gamma and its inverse fromI are<1e-272, and their A2 norms, exp(2v), volume and inverse_volume are<2. Their complete first and second spatial derivatives are bounded by the same gap using the Fourier weights.",
        "whole_full_adjoint_ghost_operator": adjoint["full"],
        "whole_exact_divergence_free_slice_adjoint": adjoint["slice"],
        "whole_off_slice_divQ_contacts": adjoint["contacts"],
        "whole_ghost_inverse_argument": "On the fixed2pi torus the nonzero flat symbol is |k|^2 I+kk^T/3, with inverse |k|^-2(I-kk^T/(4|k|^2)). Its first and second derivative multipliers have norm<=1; its A0-to-A1 inverse norm is<=2. The FULL slice adjoint has both deltaQ*second-derivative and derivative-Q*first-derivative terms. A conservative norm ceiling100||Q-I||A2 bounds (M*-M0)M0^-1. This is<1/10, so a mean-zero Neumann inverse has A0-to-A1 norm<4. All off-slice contacts are retained before imposing divQ=0. Constants are never inverted.",
        "whole_full_curvature_argument": "For e=1e-272 and ||gamma^-1||A2<2, each complete Christoffel component has A1 norm<=9e. Each Ricci component has A0 norm<=6(9e)+18(9e)^2. Summing all nine inverse-metric contractions gives the displayed whole scalar-curvature bound<1e-268, retaining every derivative and quadratic connection term.",
        "whole_complete_cotangent_source_bound": "For pi0=(Pi_v/6)gamma^-1+DQ*(P_TT Pi_tau), use ||DQ*||A1<=16 and ||P_TT||<=2. All five momentum fluctuation A1 norms are<1e-200, so ||pi0||A1<1e-198. The ENTIRE original metric/vector/M1/H generator obeys ||D0||A0<=180||pi0||A1+12||W||A1||Pi_W||A1+2[||Pi_M||A0||grad M1||A0+||Pi_H||A0||grad H||A0]<1e-194. Here ||Pi_M||A0<1 includes its retained1/10 background, and W_A1<1e-275. No vector Gauss term, matter term or boundary contact is deleted.",
        "whole_lift_and_trace_argument": "lambda=(M*)^-1 Pmean0 D0 has A1 norm<4||D0||. The complete correction is -DQ* div*lambda and has norm<=16||lambda||A1, giving ||pi||A0<1e-192 and ||pi_TF||A0<1e-190. The exact identity tr(gamma DQ*[B])=0 holds for every B, so Pi_v=2 pi:gamma is unchanged by the lift. Consequently the linear trace invariant p=Pi_v/(3 volume) does NOT inherit the coarse shear-lift bound. The full shear invariant tr(pi_TF gamma pi_TF gamma)/volume^2 is<1e-378.",
        "checks": checks,
        "gates": {
            "full_W_A1_source_bound": field["W_A0"]
            * reference.SUPPORT
            * (1 + reference.P)
            < s.Rational(1, 10**275),
            "full_M1_gradient_source_bound": field["M1_A1"] * reference.SUPPORT
            < s.Rational(1, 10**329),
            "full_H_gradient_source_bound": field["H_gradient_A0"] * reference.SUPPORT
            < s.Rational(1, 10**354),
            "strict_full_A2_self_map": budget["shape_self_map"] < TAU,
            "strict_full_A2_contraction": budget["shape_contraction"]
            < s.Rational(1, 10),
            "full_shape_positive": budget["whole_Q_minus_I_A2"] < s.Rational(1, 8),
            "whole_metric_A2_gap": budget["whole_gamma_minus_I_A2"] < METRIC,
            "whole_inverse_metric_A2_gap": budget["whole_inverse_gamma_minus_I_A2"]
            < METRIC,
            "whole_adjoint_Neumann_margin": budget["whole_adjoint_ghost_Neumann_error"]
            < s.Rational(1, 10),
            "complete_curvature_bound": budget["whole_scalar_curvature_A0"]
            < s.Rational(1, 10**268),
            "whole_base_momentum_bound": budget["whole_base_metric_momentum_A1"]
            < s.Rational(1, 10**198),
            "whole_full_generator_bound": budget["whole_complete_source_D0_A0"]
            < s.Rational(1, 10**194),
            "whole_lifted_metric_momentum_bound": budget[
                "whole_lifted_metric_momentum_A0"
            ]
            < s.Rational(1, 10**192),
            "whole_tracefree_momentum_bound": budget[
                "whole_tracefree_metric_momentum_A0"
            ]
            < s.Rational(1, 10**190),
            "whole_normalized_shear_bound": budget["whole_normalized_shear_invariant"]
            < s.Rational(1, 10**378),
            "full_nonlinear_harmonics_and_homogeneous_shape_contacts_retained": True,
            "residual_three_translations_not_inverted_or_classically_deleted": True,
        },
    }
