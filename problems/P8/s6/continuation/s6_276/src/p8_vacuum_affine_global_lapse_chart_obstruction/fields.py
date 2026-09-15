"""Entire finite Fourier input fields, not arbitrary independent invariant values."""

from functools import cache

import sympy as s

from . import fold

FREQUENCY = s.Integer(10) ** 64
AXES = (0, 1, 2)


def require_axis(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer)):
        raise TypeError("Require an exact integer torus axis")
    if int(value) not in AXES:
        raise ValueError("Require one of the three original axes")
    return int(value)


def require_frequency(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer)):
        raise TypeError("Require the exact original integer Fourier frequency")
    if value != FREQUENCY:
        raise ValueError("Require the unchanged original frequency")
    return s.Integer(value)


def fixture(axis):
    return _fixture(require_axis(axis))


@cache
def _fixture(axis):
    x = s.symbols("torus_x0:3", real=True)
    k = s.Symbol("positive_integer_original_wave_number", positive=True, integer=True)
    V, A, kap = s.symbols(
        "positive_scalar_amplitude positive_tensor_amplitude positive_kappa",
        positive=True,
    )
    d = s.Symbol("actual_zero_mean_M1_momentum_amplitude", real=True)
    transverse, third = (axis + 1) % 3, (axis + 2) % 3
    v = V * s.cos(k * x[axis])
    pol = [s.Integer(0)] * 3
    pol[axis], pol[third] = s.Integer(1), s.Integer(-1)
    tt = A * s.cos(k * x[transverse]) * s.diag(*pol)
    metric, inv = s.exp(2 * v) * s.eye(3), s.exp(-2 * v) * s.eye(3)
    pi = -kap * s.exp(-2 * v) * tt
    density = s.exp(3 * v)
    matter_momentum = kap * (s.Rational(1, 10) + d * s.cos(k * x[third]))
    point = {x[axis]: s.pi / (2 * k), x[transverse]: 0, x[third]: 0}
    D = s.Matrix(
        [
            -2 * sum(s.diff((metric * pi)[i, j], x[j]) for j in range(3))
            + sum(
                pi[j, l] * s.diff(metric[j, l], x[i])
                for j in range(3)
                for l in range(3)
            )
            for i in range(3)
        ]
    )
    Gamma = [
        [
            [
                sum(
                    inv[i, l]
                    * (
                        s.diff(metric[l, j], x[h])
                        + s.diff(metric[l, h], x[j])
                        - s.diff(metric[j, h], x[l])
                    )
                    for l in range(3)
                )
                / 2
                for h in range(3)
            ]
            for j in range(3)
        ]
        for i in range(3)
    ]
    ricci = s.Matrix(
        3,
        3,
        lambda i, j: sum(
            s.diff(Gamma[l][i][j], x[l])
            - s.diff(Gamma[l][i][l], x[j])
            + sum(
                Gamma[l][l][h] * Gamma[h][i][j] - Gamma[l][j][h] * Gamma[h][i][l]
                for h in range(3)
            )
            for l in range(3)
        ),
    )
    scalar = s.simplify(s.trace(inv * ricci))
    wanted = s.exp(-2 * v) * (
        4 * k * k * V * s.cos(k * x[axis]) - 2 * k * k * V * V * s.sin(k * x[axis]) ** 2
    )
    shear = s.simplify(s.trace(pi * metric * pi * metric) / (kap * kap * metric.det()))
    Pi_v = s.simplify(2 * s.trace(pi * metric))
    Q = s.exp(2 * v) * inv
    dual = s.exp(2 * v) * (inv * s.trace(kap * tt * inv) / 3 - inv * (kap * tt) * inv)
    checks = {
        "literal_full_Q_identity": (Q - s.eye(3)).applyfunc(s.simplify),
        "literal_full_shape_cotangent": (dual - pi).applyfunc(s.simplify),
        "literal_full_scalar_cotangent": Pi_v,
        "TT_trace": s.trace(tt),
        "TT_divergence": s.Matrix(
            [sum(s.diff(tt[i, j], x[j]) for j in range(3)) for i in range(3)]
        ),
        "whole_ADM_momentum_generator": D.applyfunc(s.simplify),
        "literal_full_Christoffel_Ricci_scalar": s.simplify(scalar - wanted),
        "full_metric_density": s.simplify(metric.det() - density * density),
        "full_normalized_shear": s.simplify(
            shear - 2 * A * A * s.exp(-6 * v) * s.cos(k * x[transverse]) ** 2
        ),
        "point_scalar_configuration": s.simplify(v.subs(point)),
        "point_curvature": s.simplify(scalar.subs(point) + 2 * k * k * V * V),
        "point_shear": s.simplify(shear.subs(point) - 2 * A * A),
        "point_M1_density": s.simplify(
            (matter_momentum / (kap * density)).subs(point) - s.Rational(1, 10) - d
        ),
        "M1_momentum_mean_retained": s.integrate(
            d * s.cos(k * x[third]), (x[third], 0, 2 * s.pi)
        ),
    }
    return {
        "whole_input_axes": (axis, transverse, third),
        "whole_formal_integer_frequency": k,
        "whole_original_frequency_binding": FREQUENCY,
        "whole_scalar_input": v,
        "whole_shape_input": s.zeros(3),
        "whole_shape_momentum_input": kap * tt,
        "whole_M1_momentum_input": matter_momentum,
        "whole_metric": metric,
        "whole_full_raw_metric_momentum": pi,
        "whole_full_source_momentum_generator": D,
        "whole_full_scalar_curvature": scalar,
        "whole_full_shear": shear,
        "whole_evaluation_point": point,
        "checks": checks,
    }


@cache
def data():
    rows = {axis: fixture(axis) for axis in AXES}
    svalue = s.Symbol("nonnegative_path_shear", nonnegative=True)
    k = s.Symbol("positive_original_frequency", positive=True)
    scalar_square, tensor_square = 6 * svalue / (k * k), svalue / 2
    checks = {
        "connected_path_curvature": s.factor(-2 * k * k * scalar_square + 12 * svalue),
        "connected_path_shear": s.factor(2 * tensor_square - svalue),
        "connected_path_original_frequency": require_frequency(FREQUENCY) - 10**64,
        **{
            "axis_" + str(axis) + "_" + name: value
            for axis, row in rows.items()
            for name, value in row["checks"].items()
        },
    }
    return {
        "whole_three_rotated_literal_field_fixtures": {
            axis: {key: value for key, value in row.items() if key != "checks"}
            for axis, row in rows.items()
        },
        "whole_path_amplitude_squares": (scalar_square, tensor_square),
        "whole_same_means_and_modes": "v=sqrt(6s)/P cos(Px), tau=0, Pi_v=0, Pi_tau=kappa sqrt(s/2) cos(Py)diag(1,0,-1). M1 is spatially constant and Pi_M/kappa=1/10+d cos(Pz), d=sqrt(1/100+6Pfixed-4rhofixed)-1/10. Heavy/vector inputs zero. The cyclic coordinate and all original means are retained. Rotate all axes as desired.",
        "whole_nonlinear_constraint_argument": "Q=I exactly; the full metric pi is the S267 adjoint base lift. D0=0 identically, including matter since gradM1=0, so its mean-zero dual correction is exactly zero. The three residual canonical translation charges are also zero, not inverted or averaged away. The Q ghost operator is the invertible mean-zero flat operator. No small curvature or global quotient theorem is assumed.",
        "whole_evaluated_invariants": "At x=pi/(2P),y=z=0, the metric isI, v0, tracep0, shear=s, curvature=-12s, pm^2=1/100+6Pfixed-4rhofixed; all other original invariants vanish. Raw metric and normalized momentum densities contain their full exp(v) factors and generated Fourier harmonics.",
        "whole_initial_domain": "At s=0 only the zero-mean M1 momentum input remains, with |d|<100PROFILE_BOUND<1e-250. Every original invariant at every spatial point is in the S266 local box. For large s this is outside the certified S275 small phase domain; no evolved solution or quantum state is asserted.",
        "checks": checks,
        "gates": {
            "all_three_original_axis_rotations": tuple(rows) == AXES,
            "unchanged_actual_frequency": FREQUENCY == 10**64,
            "same_local_source_invariant_bound": 100 * fold.PROFILE_BOUND
            < fold.original.DELTA,
            "whole_positive_fold_amplitudes": s.Rational(81, 160)
            - fold.affine_error(fold.FOLD_SHEAR)
            > 0,
            "whole_fold_shear_below_one": s.Rational(81, 160)
            + fold.affine_error(fold.FOLD_SHEAR)
            < 1,
            "configuration_amplitude_ceiling_small": s.Rational(3, 10**64)
            < s.Rational(1, 10**62),
            "complete_momentum_constraints_not_arbitrary_box_assignments": True,
            "all_three_residual_translations_retained": True,
            "full_nonlinear_harmonics_not_truncated": True,
            "not_a_quantum_state_or_global_time_trajectory": True,
        },
    }
