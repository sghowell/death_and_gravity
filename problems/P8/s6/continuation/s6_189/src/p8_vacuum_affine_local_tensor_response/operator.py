"""Full spatial local tensor Euler operator and its actual weighted adjoint."""

from functools import cache

import sympy as s
from p8_vacuum_affine_proca_gaussian.bridge import KAPPA, MASS

t = s.Symbol("t", real=True)
k = s.Symbol("spatial_momentum", nonnegative=True)
a = (1 + t * t) ** 2
H = s.diff(a, t) / a
R0 = 6 * (s.diff(H, t) + 2 * H * H)
A = s.Rational(5, 12) * MASS**2 - R0 / 36


@cache
def data():
    h = s.Function("tensor")(t)
    q = k * k / a**2
    D = lambda x: s.diff(x, t, 2) + H * s.diff(x, t) + q * x
    Dadj = lambda x: (
        s.diff(x, t, 2) + 5 * H * s.diff(x, t) + (2 * s.diff(H, t) + 6 * H * H + q) * x
    )
    L = lambda x: s.diff(x, t, 2) + 3 * H * s.diff(x, t) + q * x
    LA = A * L(h) + s.diff(A, t) * s.diff(h, t)
    composition = Dadj(D(h))
    hp = s.diff(H, t)
    hpp = s.diff(H, t, 2)
    expanded = (
        s.diff(h, t, 4)
        + 6 * H * s.diff(h, t, 3)
        + (4 * hp + 11 * H * H + 2 * q) * s.diff(h, t, 2)
        + (hpp + 7 * H * hp + 6 * H**3 + 2 * H * q) * s.diff(h, t)
        + q * q * h
    )
    action = a**3 * (A * (s.diff(h, t) ** 2 - q * h * h) - D(h) ** 2 / 60)
    EL = (
        s.diff(action, h)
        - s.diff(s.diff(action, s.diff(h, t)), t)
        + s.diff(s.diff(action, s.diff(h, t, 2)), t, 2)
    )
    f = s.Function("adjoint_test")(t)
    # Direct formal adjoint in the a^3 measure, including time derivatives of H.
    direct_adj = (
        s.diff(a**3 * f, t, 2) - s.diff(a**3 * H * f, t) + a**3 * q * f
    ) / a**3
    canonical = s.Symbol("canonical_tensor", real=True)
    gamma = 2 * canonical / s.sqrt(KAPPA)
    checks = {
        "actual_weighted_adjoint": s.simplify(direct_adj - Dadj(f)),
        "complete_fourth_order_spatial_composition": s.simplify(composition - expanded),
        "weighted_A_wave_operator": s.simplify(
            LA - (s.diff(a**3 * A * s.diff(h, t), t) / a**3 + A * q * h)
        ),
        "independent_literal_action_Euler_variation": s.simplify(
            EL + 2 * a**3 * (LA + composition / 60)
        ),
        "canonical_quadratic_action_factor": s.simplify(
            gamma * gamma / (64 * s.pi**2)
            - canonical * canonical / (16 * s.pi**2 * KAPPA)
        ),
        "fixed_local_A_curvature_term": s.simplify(
            A - (s.Rational(5, 12) * MASS**2 - R0 / 36)
        ),
        "actual_background_old_curvature": s.simplify(
            R0 - 24 * (1 + 7 * t * t) / (1 + t * t) ** 2
        ),
    }
    return {
        "quadratic_gamma_action": "(1/(64pi^2)) integral a^3 {A tr[gamma_t^2-a^-2(grad gamma)^2]-(1/60)tr[(D gamma)^2]}, A=5m^2/12-R0/36",
        "actual_A": s.factor(A),
        "actual_R0": s.factor(R0),
        "wave_operators": "L=tt+3H t-a^-2 Delta, D=tt+H t-a^-2 Delta, Ddag=tt+5H t+2H'+6H^2-a^-2 Delta in the a^3 measure",
        "full_canonical_local_response": "Qloc=(1/(8pi^2 kappa))[A L+A'partial_t+(1/60)Ddag D]; the total Euler equation reads -Lh-Qloc h+source=0",
        "full_DdagD": "t^4+6H t^3+(4H'+11H^2)t^2+(H''+7HH'+6H^3)t-2a^-2 Delta t^2-2H a^-2 Delta t+a^-4 Delta^2; the pure Delta coefficient cancels",
        "response_boundary": "This is the complete Hessian of the specified finite local matching density, not the full determinant response, its inverse, an order-reduced solution or a discarded-pole prescription.",
        "checks": checks,
    }
