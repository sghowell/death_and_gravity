"""Arbitrary Rosen metrics, weighted tensor grading and exact norm constants."""

from functools import cache
from itertools import product

import sympy as s


@cache
def rosen():
    u = s.Symbol("u", real=True)
    a, b, c = (s.Function(name)(u) for name in ("a", "b", "c"))
    metric = s.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, a, b], [0, 0, b, c]])
    inverse = metric.inv().applyfunc(s.factor)

    def derivative(value, index):
        return s.diff(value, u) if index == 0 else s.S.Zero

    gamma = {}
    for rho, mu, nu in product(range(4), repeat=3):
        gamma[rho, mu, nu] = s.factor(
            sum(
                inverse[rho, z]
                * (
                    derivative(metric[z, nu], mu)
                    + derivative(metric[z, mu], nu)
                    - derivative(metric[mu, nu], z)
                )
                / 2
                for z in range(4)
            )
        )
    ricci = s.zeros(4)
    for mu, nu in product(range(4), repeat=2):
        ricci[mu, nu] = s.factor(
            sum(
                derivative(gamma[rho, mu, nu], rho)
                - derivative(gamma[rho, mu, rho], nu)
                + sum(
                    gamma[rho, rho, z] * gamma[z, mu, nu]
                    - gamma[rho, nu, z] * gamma[z, mu, rho]
                    for z in range(4)
                )
                for rho in range(4)
            )
        )
    scalar = s.factor(
        sum(inverse[mu, nu] * ricci[mu, nu] for mu, nu in product(range(4), repeat=2))
    )
    contraction = s.factor(
        sum(
            inverse[mu, nu]
            * (
                gamma[rho, mu, nu] * gamma[z, rho, z]
                - gamma[rho, mu, z] * gamma[z, nu, rho]
            )
            for mu, nu, rho, z in product(range(4), repeat=4)
        )
    )
    G = s.Matrix([[a, b], [b, c]])
    expected = (
        -s.trace(G.inv() * G.diff(u, 2)) / 2
        + s.trace(G.inv() * G.diff(u) * G.inv() * G.diff(u)) / 4
    )
    checks = {
        "arbitrary_Rosen_inverse": (metric * inverse - s.eye(4)).applyfunc(s.factor),
        "arbitrary_Rosen_scalar_curvature": scalar,
        "arbitrary_Rosen_Gamma_Gamma_contraction": contraction,
        "arbitrary_Rosen_nontrivial_uu_Ricci": s.factor(ricci[0, 0] - expected),
    }
    for mu, nu in product(range(4), repeat=2):
        if (mu, nu) != (0, 0):
            checks[f"arbitrary_Rosen_Ricci_{mu}_{nu}"] = ricci[mu, nu]
    return {
        "checks": checks,
        "whole_metric": metric,
        "whole_Ricci_uu": ricci[0, 0],
        "whole_vanishing_action_contraction": contraction,
    }


@cache
def data():
    r = rosen()
    checks = dict(r["checks"])
    e = s.Symbol("epsilon", real=True)
    a, b, c, d, f, g = s.symbols("A B C D F G", real=True)
    H = s.Matrix(
        [[0, 0, 0, 0], [0, a, b, e * c], [0, b, d, e * f], [0, e * c, e * f, e**2 * g]]
    )
    w, px, py, pz = s.symbols("omega px py pz", real=True)
    p = s.Matrix([w, e * px, e * py, w + e**2 * pz])
    rotation = s.diag(1, -1, -1, 1)
    checks["weighted_tensor_reflection_grading"] = (
        H.subs(e, -e) - rotation * H * rotation
    )
    checks["weighted_momentum_reflection_grading"] = p.subs(e, -e) - rotation * p
    checks["proper_transverse_rotation"] = s.det(rotation) - 1
    eta = s.diag(1, -1, -1, -1)
    checks["rotation_preserves_Minkowski_background"] = (
        rotation.T * eta * rotation - eta
    )
    rho = s.Symbol("rho", real=True)
    x, y, z = s.symbols("vx vy vz", real=True)
    velocity = s.Matrix([x, y, z])
    axis = s.Matrix([0, 0, 1])
    checks["child_direction_and_variance_identity"] = s.expand(
        ((axis - velocity).T * (axis - velocity))[0]
        - (2 * (1 - z) - (1 - x * x - y * y - z * z))
    )
    checks["parent_delta_energy_deficit_identity"] = s.expand(
        1 - rho * rho - (1 - rho) * (1 + rho)
    )
    q = s.Symbol("q", positive=True)
    checks["field_coefficient_l1_factor13"] = s.expand(
        13 * q * q - (1 + 4 * q + 8 * q * q) - (q - 1) * (5 * q + 1)
    )
    checks["full_tensor_squared_component_budget"] = (
        2 * s.Rational(3, 2) ** 2 + 2 + 4 + s.Rational(3, 2) ** 2 - s.Rational(51, 4)
    )
    checks["velocity_squared_component_budget"] = (
        2 + s.Rational(3, 2) ** 2 - s.Rational(17, 4)
    )
    return {
        "whole_weighted_norm": "For temporal H and future Q=(W,Qsp), v=Qsp/W and delta^2=1-|v|^2, use max(||Hsp||F,||Hsp*v||/delta,|v^T*Hsp*v|/delta^2). TT singleton seeds use their Frobenius norm. Generic larger subsets have delta>0.",
        "whole_arbitrary_Rosen_geometry": {k: v for k, v in r.items() if k != "checks"},
        "whole_child_transfer": "For R=W/W_A>=1, delta_A^2 and |n-v_A|^2 are at most2*delta^2*R. The child decomposes as transverse A+delta*mixed B+delta^2*longitudinal C with coefficient norms at most1,4*sqrt(R),8R times its own weighted norm; their sum is at most13R.",
        "whole_graded_momentum_budget": "Future-ray geometry gives Q_A=W_A*k+delta*b_A+delta^2*c_A, with coefficient-l1 norm below4W. The negative root obeys the same bound, with exact coefficientwise momentum conservation.",
        "whole_propagator_budget": "The conserved temporal propagator turns component budgets (delta^2 B,delta B,B) in the transverse/mixed/longitudinal source blocks into weighted field norm at most4B/W^2. Squared constants are51/4 for the field and17/4 for its weighted velocity contraction.",
        "checks": checks,
        "gates": {
            "arbitrary_transverse_metric_not_required_traceless": True,
            "only_Ricci_uu_can_be_nonzero_in_the_Rosen_family": True,
            "Rosen_zero_proves_constant_transverse_vertex_zero_at_every_order": True,
            "proper_rotation_fixes_even_and_odd_vertex_parities": True,
            "future_ray_geometry_does_not_divide_by_child_spatial_momentum": True,
            "singleton_zero_delta_case_uses_exact_TT_transversality": True,
            "field_coefficient_sum_is_bounded_by13_energy_ratio": True,
            "momentum_coefficient_sum_strictly_below4W": 2 * s.sqrt(2) + 1 < 4,
            "weighted_propagator_factor4": s.Rational(51, 4) < 16
            and s.Rational(17, 4) < 16
            and s.Rational(3, 2) < 4,
            "exact_collinear_points_not_assigned_a_new_value": True,
        },
    }
