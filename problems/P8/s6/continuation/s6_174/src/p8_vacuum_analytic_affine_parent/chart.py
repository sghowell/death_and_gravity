"""Regular covariant metric chart and exact generic-R ADM cancellations."""

from functools import cache

import sympy as s
from p8_affine_vacuum_domain import family
from p8_vacuum_flat_dirac_hadamard.symbols import rational


def require_domain(X, R, zeta):
    x, r, z = map(rational, (X, R, zeta))
    if not -s.Rational(1, 4096) < x < s.Rational(6, 5):
        raise ValueError("Outside the declared analytic coefficient strip")
    if not s.Rational(1, 2) < r < s.Rational(6, 5):
        raise ValueError("Outside the proven positive tensor-factor enclosure")
    if not 0 < z <= s.Rational(1, 2000):
        raise ValueError("Require the stated positive vector frequency regime")
    return x, r, z


@cache
def data():
    u, X = s.symbols("u X", real=True)
    R = s.Function("positive_R")(u, X)
    RX, Ru = s.diff(R, X), s.diff(R, u)
    # g_physical = C*g_hat + D*du*du, with X_hat=X_physical.
    C = R ** -s.Rational(1, 2)
    disformal = (1 - C) / X
    omega = -s.log(R) / 4
    B = -R / 2
    c = -s.sqrt(X) * RX
    d = -3 * X * RX**2 / (4 * R)
    e = -X * RX + 7 * X**2 * RX**2 / (4 * R)
    fphi = -Ru / 2
    v, K, KK, V = s.symbols("positive_s Khat Khat_ij_squared V", real=True)
    shift = v * s.diff(omega, u) + 2 * v * s.diff(omega, X) * V
    physicalK = K + 3 * shift
    physicalKK = KK + 2 * shift * K + 3 * shift**2
    literal = (
        B * (physicalK**2 - physicalKK)
        + c * physicalK * V
        + d * V**2
        + 2 * v * fphi * physicalK
    )
    literal = literal.subs(s.sqrt(X), v).subs(X, v * v)
    bx = (4 * B * v * s.diff(omega, u) + 2 * v * fphi).subs(X, v * v)
    fx = (
        6 * B * v * v * s.diff(omega, u) ** 2 + 6 * v * v * fphi * s.diff(omega, u)
    ).subs(X, v * v)
    lv = (12 * v * v * fphi * s.diff(omega, X)).subs(X, v * v)
    expected = B.subs(X, v * v) * (K * K - KK) + bx * K + fx + lv * V
    delta = -2 * X * s.diff(omega, X)
    B4 = R / 2
    gradient = e + 2 * B4 * delta**2 + 4 * (B4 - 2 * X * s.diff(B4, X)) * delta
    y = s.symbols("positive_sqrt_R", positive=True)
    a, b = s.symbols("a b", real=True)
    xsmall = s.symbols("xsmall", real=True)
    germ = 1 + a * xsmall * xsmall + b * xsmall**3
    cg = s.series(germ ** -s.Rational(1, 2), xsmall, 0, 4).removeO()
    dg = s.cancel((1 - cg) / xsmall)
    return {
        "metric_map": "g_physical=C*g_hat+D*du*du, C=R^(-1/2), D=(1-C)/X; the physical matter metric remains g_physical",
        "C": C,
        "D": disformal,
        "inverse_map": "g_hat=C^(-1)*g_physical+(1-C^(-1))*du*du/X; both coefficient quotients extend analytically through X=0",
        "determinant_ratio_physical_to_hat": R ** -s.Rational(3, 2),
        "X_map": "X_physical=X_hat exactly; the algebraic metric map has explicit inverse and nonzero field-space Jacobian",
        "metric_jacobian_determinant_at_fixed_du": C**9,
        "time_linear_boundary_primitive": "I(u,s)=int_1^s R(u,r^2)^(-3/4)*[3*r^2*R_u*R_X/(2R)](u,r^2) dr",
        "time_boundary_rule": "sqrt(hhat)*I_s*(dot s-N^i partial_i s)=total_divergence-N*sqrt(hhat)*[I*Khat+s*I_u]",
        "remaining_linear_V": lv,
        "transformed_trace_coefficient_per_hat_volume": -(R ** s.Rational(1, 4)) / 3,
        "spatial_lapse_gradient_after_full_boundary": s.simplify(gradient),
        "vacuum_C_germ": cg,
        "vacuum_D_germ": dg,
        "checks": {
            "norm_denominator_exactly_one": s.simplify(C + X * disformal - 1),
            "inverse_metric_coefficient": s.simplify((1 - 1 / C) / X + disformal / C),
            "inverse_norm_derivative_one": s.simplify(
                s.diff(X / (C + X * disformal), X) - 1
            ),
            "metric_rank_one_jacobian_factor": s.simplify(
                C - X * s.diff(C, X) - X * X * s.diff(disformal, X) - 1
            ),
            "full_generic_R_ADM_time_action": s.simplify(literal - expected),
            "mixed_K_V_removed": s.simplify(s.diff(literal, K, V)),
            "V_squared_removed": s.simplify(s.diff(literal, V, 2)),
            "generic_R_spatial_gradient_removed": s.simplify(gradient),
            "regular_D_rationalization": s.factor(
                (1 - 1 / y) / X - (y * y - 1) / (X * y * (y + 1))
            ),
            "vacuum_D_zero": dg.subs(xsmall, 0),
            "vacuum_D_first_germ": s.diff(dg, xsmall).subs(xsmall, 0) - a / 2,
            "volume_trace_coefficient": s.simplify(
                R ** -s.Rational(3, 4) * (-R / 2) * s.Rational(2, 3)
                + R ** s.Rational(1, 4) / 3
            ),
        },
    }


@cache
def target_jets():
    d = family.data()
    u, X = family.u, family.X
    h = (1 + u * u) ** 3
    old = family.original.data()["original_retuned_tree_scalar"]
    checks = {}
    for j in range(4):
        checks["R_clock_jet_" + str(j)] = s.simplify(
            s.diff(d["R"] - (1 + (X - 1) / h), X, j).subs(X, 1)
        )
    for j in range(3):
        checks["F_clock_jet_" + str(j)] = s.simplify(
            s.diff(d["F"] - old, X, j).subs(X, 1)
        )
    checks["R_vacuum_value"] = d["R"].subs(X, 0) - 1
    checks["R_vacuum_first_X"] = s.diff(d["R"], X).subs(X, 0)
    checks["R_vacuum_second_X"] = s.simplify(
        s.diff(d["R"], X, 2).subs(X, 0) + 2 * family.N / h
    )
    vacuum = {u: 0, X: 0}
    checks["full_F_vacuum_value"] = d["F"].subs(vacuum)
    checks["full_F_vacuum_tadpole"] = s.diff(d["F"], u).subs(vacuum)
    checks["full_F_vacuum_mass_one"] = s.diff(d["F"], u, 2).subs(vacuum) + 1
    checks["full_F_vacuum_positive_kinetic"] = s.diff(d["F"], X).subs(
        vacuum
    ) - s.Rational(1, 2)
    return {
        "checks": checks,
        "same_target": "Unchanged S6.109 full analytic CD_matter action with original rolling M1 and retuned margin1/200",
        "global_domain": "Every real u, -1/4096<X<6/5, with 1/2<R<6/5",
        "clock_jet_order_R": 3,
        "clock_jet_order_F": 2,
        "vacuum_R_minus_one_has_X_squared_factor": True,
    }
