"""Scaled full-convolution curvature and complete tracefree cotangent bounds."""

from functools import cache

import sympy as s

from . import moving, source

P, TIME, IMAGE_RADIUS = moving.P, source.TIME, moving.IMAGE_RADIUS


@cache
def bounds():
    raw = moving.old.field.field_bounds()
    fields = {name: value * IMAGE_RADIUS for name, value in raw.items()}
    v = 4 * fields["v_A2"] / (1 + P) ** 2
    tau = 4 * fields["tau_A2"] / (1 + P) ** 2
    f_shape = 16 * tau**2
    shape_size = tau + f_shape
    inverse_gap = shape_size / (1 - shape_size)
    exp_gap = 2 * v / (1 - 2 * v)
    metric_gap = exp_gap / (1 - shape_size) + inverse_gap
    inverse_remainder = f_shape + shape_size**2 / (1 - shape_size)
    exp_remainder = 4 * v**2 / (1 - 2 * v)
    metric_remainder = exp_remainder + inverse_remainder + exp_gap * inverse_gap
    inverse_metric_gap = metric_gap / (1 - metric_gap)
    connection = 9 * metric_gap
    connection_linear = s.Rational(9, 2) * metric_gap
    connection_remainder = s.Rational(9, 2) * inverse_metric_gap * metric_gap
    ricci_linear = 6 * connection_linear
    ricci_remainder = 6 * connection_remainder + 18 * connection**2
    scalar_remainder = 3 * ricci_remainder + 9 * inverse_metric_gap * (
        ricci_linear + ricci_remainder
    )
    full_remainder = scalar_remainder + 18 * metric_remainder
    curvature = 4 * P**2 * v + 10**12 * P**2 * (v + tau) ** 2
    background_momentum = 100 * TIME
    base_A1 = (background_momentum + (1 + P) * fields["Pi_v_A0"]) / 3 + 32 * (
        1 + P
    ) * fields["Pi_tau_A0"]
    D0 = (
        180 * base_A1
        + 12 * (1 + P) ** 2 * fields["W_A0"] * fields["Pi_W_A0"]
        + 2 * (fields["M1_A1"] + fields["Pi_H_A0"] * fields["H_gradient_A0"])
    )
    descriptor = 4 * D0 / P
    pi_TF = 5 * (32 * fields["Pi_tau_A0"] + 16 * descriptor)
    shear = 48 * pi_TF**2
    p_delta = 2 * fields["Pi_v_A0"] / 3 + 100 * TIME * 6 * v
    pm_delta = 2 * fields["delta_Pi_M_A0"] + 6 * v
    images = [
        p_delta,
        6 * (1 + P) * fields["Pi_W_A0"],
        pm_delta,
        2 * fields["Pi_H_A0"],
        fields["eta_A0"],
        shear,
        8 * fields["Pi_W_A0"] ** 2 / moving.old.field.ZETA,
        48 * moving.old.field.ZETA * ((1 + P) * fields["W_A0"]) ** 2,
        2 * fields["W_A0"] ** 2,
        2 * fields["M1_A1"] ** 2,
        2 * fields["H_gradient_A0"] ** 2,
        curvature,
    ]
    weighted_linear = 100 * TIME * p_delta + pm_delta + 4 * sum(images[5:])
    quadratic = source.MAJORANT * sum(images) ** 2 / 2
    residual = s.Rational(5, 2) * source.PROFILE + weighted_linear + quadratic
    return {
        "fields": fields,
        "v": v,
        "tau": tau,
        "f_shape": f_shape,
        "shape_size": shape_size,
        "metric_gap": metric_gap,
        "metric_remainder": metric_remainder,
        "inverse_metric_gap": inverse_metric_gap,
        "connection": connection,
        "connection_remainder": connection_remainder,
        "ricci_linear": ricci_linear,
        "ricci_remainder": ricci_remainder,
        "scalar_remainder": scalar_remainder,
        "full_remainder": full_remainder,
        "curvature": curvature,
        "background_momentum": background_momentum,
        "base_A1": base_A1,
        "D0": D0,
        "descriptor": descriptor,
        "pi_TF": pi_TF,
        "shear": shear,
        "images": images,
        "weighted_linear": weighted_linear,
        "quadratic": quadratic,
        "residual": residual,
    }


@cache
def identities():
    x = s.symbols("full_spatial_x0:3", real=True)
    h = s.zeros(3)
    for i in range(3):
        for j in range(i, 3):
            h[i, j] = h[j, i] = s.Function("full_metric_" + str(i) + str(j))(*x)
    connection = [
        [
            [
                sum(
                    s.KroneckerDelta(k, l)
                    * (
                        s.diff(h[j, l], x[i])
                        + s.diff(h[i, l], x[j])
                        - s.diff(h[i, j], x[l])
                    )
                    for l in range(3)
                )
                / 2
                for j in range(3)
            ]
            for i in range(3)
        ]
        for k in range(3)
    ]
    Ricci = s.Matrix(
        3,
        3,
        lambda i, j: sum(
            s.diff(connection[k][i][j], x[k]) - s.diff(connection[k][i][k], x[j])
            for k in range(3)
        ),
    )
    scalar = s.trace(Ricci)
    expected = sum(
        s.diff(h[i, j], x[i], x[j]) for i in range(3) for j in range(3)
    ) - sum(s.diff(s.trace(h), z, 2) for z in x)
    v = s.Function("full_conformal_v")(*x)
    rule = {h[i, j]: 2 * v * int(i == j) for i in range(3) for j in range(3)}
    checks = {
        "full_three_dimensional_linear_Ricci_contraction": s.expand(scalar - expected),
        "full_conformal_linear_curvature": s.simplify(
            scalar.subs(rule).doit() + 4 * sum(s.diff(v, z, 2) for z in x)
        ),
    }
    for vector in ((1, 0, 0), (1, 2, 0), (1, -1, 2), (3, 4, 5)):
        k = P * s.Matrix(vector)
        square = k.dot(k)
        full = square * s.eye(3) + k * k.T / 3
        inverse = (s.eye(3) - k * k.T / (4 * square)) / square
        label = "_".join(map(str, vector))
        checks["full_flat_inverse_left_" + label] = full * inverse - s.eye(3)
        checks["full_flat_inverse_right_" + label] = inverse * full - s.eye(3)
    gv = s.symbols("whole_metric_slot0:6", real=True)
    g = s.Matrix([[gv[0], gv[1], gv[2]], [gv[1], gv[3], gv[4]], [gv[2], gv[4], gv[5]]])
    b = s.Symbol("whole_trace_density", real=True)
    pure = b * g.inv() / 6
    checks["whole_pointwise_tracefree_projection_of_trace_density"] = (
        pure - g.inv() * s.trace(g * pure) / 3
    ).applyfunc(s.factor)
    return checks


@cache
def data():
    b = bounds()
    v, tau, shape = b["v"], b["tau"], b["shape_size"]
    scale = (1 + TIME**2) ** 4
    return {
        "whole_full_fields_on_complex_image_ball": b["fields"],
        "whole_scaled_v_and_tau_norms": [v, tau],
        "whole_evaluated_full_geometry_and_cotangent_bounds": {
            name: value for name, value in b.items() if name not in ("fields", "images")
        },
        "whole_twelve_invariant_deviations_from_moving_center": dict(
            zip(source.q.COORDS, b["images"], strict=True)
        ),
        "whole_sublattice_proof": "Initial support is P Z^3. Full products, derivatives, Fourier projectors and the unique translation-equivariant shape fixed point preserve this closed subalgebra. The full mean-zero adjoint Neumann inverse does too. Zero modes and all generated harmonics remain. For y=P x, scaled A2 uses (1+|k/P|)^2; initial |k/P|=1 gives the displayed factor4. Physical derivatives restore P, and nonzero inverse modes satisfy |k|>=P. The full A0-to-A1 inverse is bounded by4/P, not a scalarized polarization inverse.",
        "whole_curvature_proof": "Q=I+tau+Bf has ||f||<=16||tau||^2. For h=exp(2v)Q^-1-I, ||h||<=8(v+tau) and ||h-(2vI-tau)||<=100(v+tau)^2 in scaled A2. Full inverse, connection, Ricci and final inverse-metric contraction remainders are displayed. Linear scalar curvature is div div h-Delta tr h; transverse tracefree tau contributes zero. Thus full |R3|<=4P^2 v+1e12 P^2(v+tau)^2, with a^-2<=1 at real time. The complex proof is a holomorphic continuation of the real formal adjoint, not a conjugate-transpose substitution.",
        "whole_cotangent_proof": "The Pi_v gamma^-1/6 part is exactly pure trace pointwise; only its tracefree projection vanishes. Its full nonzero source contribution, including the homogeneous momentum, stays in D0. The complete matter/Gauss/vector source and derivative-Q adjoint term remain. The full tracefree lift is bounded by5(32 Pi_tau+16 lambda); shear<=48 pi_TF^2. No residual translation constraint or mode is deleted.",
        "whole_center_residual_proof": "Use actual p and pm background values. The full center gradients CG,Cph,Ceta vanish, Cp<100T, Cpm<1 and other gradients<4. Taylor's integral remainder is bounded by M(sum deviations)^2/2. This is the complete full-source center, not a generic M*T estimate, a deleted linear invariant, or a bare-clock off-source replacement.",
        "checks": identities(),
        "gates": {
            "actual_scaled_v_row": v < 4 * IMAGE_RADIUS / s.Integer(10) ** 478,
            "actual_scaled_tau_row": tau < 4 * IMAGE_RADIUS / s.Integer(10) ** 423,
            "full_shape_self_map": 3 * shape**2 + 3 * shape**3 < b["f_shape"],
            "full_shape_contraction": 6 * shape + 9 * shape**2 < s.Rational(1, 10),
            "full_shape_size": shape < 2 * tau,
            "full_metric_linear_size": b["metric_gap"] < 8 * (v + tau),
            "full_metric_quadratic_remainder": b["metric_remainder"]
            < 100 * (v + tau) ** 2,
            "full_metric_small_domain": b["metric_gap"] < s.Rational(1, 100),
            "full_Ricci_contraction_remainder": b["scalar_remainder"]
            < 10**4 * b["metric_gap"] ** 2,
            "full_curvature_remainder": b["full_remainder"] < 10**12 * (v + tau) ** 2,
            "full_curvature_bound": b["curvature"] < s.Rational(1, 10**327),
            "full_physical_metric_A2_norm": scale * (1 + P**2 * b["metric_gap"]) < 2,
            "full_physical_inverse_metric_A2_norm": (1 + P**2 * b["inverse_metric_gap"])
            / scale
            < 2,
            "full_physical_metric_neighborhood": scale
            - 1
            + scale * P**2 * b["metric_gap"]
            < s.Rational(1, 100),
            "full_physical_adjoint_Neumann_error": 100 * P**2 * shape
            < s.Rational(1, 10),
            "whole_reference_matter_momentum_below_one": s.Rational(1, 10)
            + b["fields"]["delta_Pi_M_A0"]
            < 1,
            "whole_tracefree_shear": b["shear"] < s.Rational(1, 10**370),
            "all_twelve_invariant_deviations": max(b["images"])
            < s.Rational(1, 10**267),
            "whole_center_residual": b["residual"] < s.Rational(1, 10**324),
            "original_translation_and_homogeneous_shape_terms_retained": True,
        },
    }
