"""Independent full four-dimensional plus-polarization curvature calculation."""

from functools import cache

import sympy as s

t, z = s.symbols("t z", real=True)
a = s.Function("a", positive=True)(t)
g = s.Function("gamma_plus")(t, z)


@cache
def curvature():
    metric = s.diag(1, -a * a * s.exp(g), -a * a * s.exp(-g), -a * a)
    inverse = metric.inv()

    def partial(f, index):
        return s.diff(f, t) if index == 0 else s.diff(f, z) if index == 3 else s.S.Zero

    connection = [
        [
            [
                s.simplify(
                    sum(
                        inverse[i, l]
                        * (
                            partial(metric[l, j], k)
                            + partial(metric[l, k], j)
                            - partial(metric[j, k], l)
                        )
                        / 2
                        for l in range(4)
                    )
                )
                for k in range(4)
            ]
            for j in range(4)
        ]
        for i in range(4)
    ]
    riemann = {}
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    value = (
                        partial(connection[i][j][l], k)
                        - partial(connection[i][j][k], l)
                        + sum(
                            connection[i][m][k] * connection[m][j][l]
                            - connection[i][m][l] * connection[m][j][k]
                            for m in range(4)
                        )
                    )
                    value = s.simplify(metric[i, i] * value)
                    if value != 0:
                        riemann[i, j, k, l] = value
    ricci = s.Matrix(
        4,
        4,
        lambda j, l: s.simplify(
            sum(inverse[i, i] * riemann.get((i, j, i, l), 0) for i in range(4))
        ),
    )
    scalar = s.simplify(sum(inverse[i, i] * ricci[i, i] for i in range(4)))
    ricci2 = s.expand(
        sum(
            s.cancel(inverse[i, i] * inverse[j, j] * ricci[i, j] ** 2)
            for i in range(4)
            for j in range(4)
        )
    )
    riemann2 = s.expand(
        sum(
            s.cancel(
                inverse[i, i] * inverse[j, j] * inverse[k, k] * inverse[l, l] * v * v
            )
            for (i, j, k, l), v in riemann.items()
        )
    )
    Weyl2 = s.simplify(riemann2 - 2 * ricci2 + scalar * scalar / 3)
    epsilon = s.Symbol("epsilon", real=True)
    quadratic = s.expand(Weyl2.subs(g, epsilon * g).doit()).coeff(epsilon, 2)
    H = s.diff(a, t) / a
    R0 = 6 * (s.diff(H, t) + 2 * H * H)
    expected_scalar = -R0 - (s.diff(g, t) ** 2 - s.diff(g, z) ** 2 / a**2) / 2
    expected_Weyl = (
        s.diff(g, t, 2) + H * s.diff(g, t) + s.diff(g, z, 2) / a**2
    ) ** 2 - 4 * s.diff(g, t, z) ** 2 / a**2
    box = s.diff(g, t, 2) + H * s.diff(g, t) - s.diff(g, z, 2) / a**2
    boundary = 4 * (
        s.diff(a * s.diff(g, t) * s.diff(g, z, 2), t)
        - s.diff(a * s.diff(g, t) * s.diff(g, t, z), z)
    )
    return {
        "scalar_P8": scalar,
        "quadratic_Weyl_squared": quadratic,
        "R_old_background": R0,
        "checks": {
            "full_nonlinear_unimodular_volume": s.simplify(metric.det() + a**6),
            "exact_dynamic_spatial_scalar_curvature": s.simplify(
                scalar - expected_scalar
            ),
            "direct_dynamic_Weyl_quadratic": s.simplify(quadratic - expected_Weyl),
            "full_cosmic_Weyl_action_boundary": s.simplify(
                a**3 * (quadratic - box * box) - boundary
            ),
        },
    }


@cache
def data():
    p = curvature()
    return {
        "direct_metric": "diag(1,-a^2 exp(g),-a^2 exp(-g),-a^2), g=g(t,z); full dynamic four-dimensional connection and curvature are contracted before expansion",
        "P8_scalar": p["scalar_P8"],
        "quadratic_Weyl_squared": p["quadratic_Weyl_squared"],
        "tensor_scalar_rule": "For compact TT gamma, the second scalar-curvature density against a time-only coefficient is (1/4)a^3 tr[gamma_t^2-a^-2(grad gamma)^2]; R_old=-R_P8 and its first TT variation is zero.",
        "tensor_Weyl_rule": "The compact quadratic integral of C^2 is (1/2) integral a^3 tr[(gamma_tt+H gamma_t-a^-2 Delta gamma)^2]. The displayed plus polarization has Frobenius norm squared2; both orthonormal TT polarizations are retained by spatial rotations and bilinear orthogonality.",
        "checks": p["checks"],
    }
