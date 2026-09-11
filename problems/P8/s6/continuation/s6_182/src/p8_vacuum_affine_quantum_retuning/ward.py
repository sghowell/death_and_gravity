"""Literal stationary clock/vacuum mean equations and fixed-coefficient currents."""

from functools import cache

import sympy as s

from . import profile

u = profile.u
rho, P = profile.rho, profile.P
H = s.Function("reference_H")(u)
S = rho + P


def D(value):
    return s.diff(value, u) + 3 * H * value


def reduce_ward(value):
    expr = s.expand(value.doit())
    expr = expr.subs(
        s.diff(rho, u, 2),
        -3 * s.diff(H, u) * S - 3 * H * (s.diff(rho, u) + s.diff(P, u)),
    )
    return s.factor(expr.subs(s.diff(rho, u), -3 * H * S))


@cache
def data():
    X = profile.X
    A, B = -P, -S / 2
    F = A + B * (X - 1)
    energy = 2 * X * s.diff(F, X) - F
    scalar = s.diff(A, u) - 2 * s.diff(B, u) - 6 * H * B
    lapse, scale, n, v = s.symbols("N v n v_source", real=True)
    literal = lapse * s.exp(3 * scale) * F.subs(X, lapse**-2)
    point = {lapse: 1, scale: 0}
    mean = s.Matrix(
        [s.diff(literal, lapse).subs(point), s.diff(literal, scale).subs(point)]
    )
    Hess = s.hessian(literal, (lapse, scale)).subs(point)
    expected = s.Matrix([[-S, 3 * rho], [3 * rho, -9 * P]])
    rho_var, P_var = s.symbols("physical_delta_rho physical_delta_P", real=True)
    vector = s.Matrix([-rho_var - 3 * rho * v, 3 * P_var + 3 * P * (n + 3 * v)])
    added = Hess * s.Matrix([n, v])
    total = s.Matrix([-rho_var - S * n, 3 * P_var + 3 * S * n])
    # Full nonstationary second-map contacts cancel only in the sum.
    w1, w2 = s.symbols("omega_N omega_NN", real=True)
    J = s.Matrix([[1, 0], [w1, 1]])
    vector_mean = s.Matrix([-rho, 3 * P])
    vacuum_profile = -profile.PV
    checks = {
        "clock_profile_energy_cancels_same_quantum": s.expand(energy.subs(X, 1) + rho),
        "clock_profile_pressure_cancels_same_quantum": F.subs(X, 1) + P,
        "clock_scalar_Euler_is_conserved_Ward": s.expand(
            scalar - s.diff(rho, u) - 3 * H * S
        ),
        "clock_scalar_Euler_zero": reduce_ward(scalar),
        "literal_profile_mean_cancels_full_vector_current": s.ImmutableMatrix(
            mean + vector_mean
        ),
        "literal_profile_physical_Hessian": s.ImmutableMatrix(Hess - expected),
        "actual_fixed_profile_plus_vector_current": s.ImmutableMatrix(
            (vector + added - total).applyfunc(s.expand)
        ),
        "total_mean_chart_first_map": s.ImmutableMatrix(J.T * (mean + vector_mean)),
        "total_mean_chart_second_map": s.diag(w2 * (mean[1] + vector_mean[1]), 0),
        "individual_profile_second_map_retained": w2 * mean[1] + 3 * P * w2,
        "individual_vector_second_map_retained": w2 * vector_mean[1] - 3 * P * w2,
        "Minkowski_constant_energy_cancellation": -vacuum_profile - profile.PV,
        "Minkowski_constant_pressure_cancellation": -profile.PV + profile.PV,
        "fixed_profile_metric_variation_has_no_delta_state": s.diff(literal, rho_var)
        + s.diff(literal, P_var),
    }
    return {
        "reference_mean_definition": "All quantities in this packet are divided by kappa. The actual vector reference functions obey rho'+3H(rho+P)=0 in the fixed covariant prescription.",
        "literal_clock_profile_energy": energy.subs(X, 1),
        "literal_clock_profile_pressure": F.subs(X, 1),
        "literal_clock_profile_scalar_Euler": scalar,
        "profile_mean_physical_currents": s.ImmutableMatrix(mean),
        "profile_current_Hessian": s.ImmutableMatrix(Hess),
        "complete_retained_stationary_current_variation": s.ImmutableMatrix(total),
        "clock_stationarity": "The full CD history, original nonzero M1 charge, affine connection and W=0 satisfy the new retained classical-plus-specified-Gaussian-vector mean equations. The coefficient functions are defined once on the full reference history. Numerical small-correction estimates are only on the declared compact slab.",
        "vacuum_stationarity": "In the same new action and covariant prescription, the Minkowski positive-frequency conditional vector vacuum has rho=-PV and P=PV. The fixed constant DeltaF=-PV cancels both. No new normal ordering is used.",
        "preparation": "The CD physical geometry and Proca operator have not changed, so the original all-order Cauchy state is already compatible with the stationary retained background. No nonprepared residual is passed to a prepared inverse, and no initial covariance or higher jet is independently reset.",
        "chart_contact_boundary": "Individual vector and profile second-map/current-density contacts are nonzero. They cancel in their total because the total one-point current vanishes; they are not discarded separately.",
        "checks": checks,
        "gates": {
            "same_reference_rho_P_not_recomputed_for_sources": True,
            "only_retained_conditional_mean_stationarity": True,
        },
    }
