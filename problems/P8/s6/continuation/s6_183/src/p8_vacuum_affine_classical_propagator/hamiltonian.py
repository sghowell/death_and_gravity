"""Source-normalized Hamiltonian inverse of the unchanged classical reference."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_response import classical

u = classical.u
H, delta, ell = classical.H, classical.delta, classical.ell
theta = u * (4 * u**6 + 12 * u**4 + 12 * u**2 + 3) / (1 + u * u) ** 4
w0 = ell * (3 * delta - 1)
L = 3 * ell * w0
P = (
    3200 * u**34
    + 48000 * u**32
    + 332800 * u**30
    + 1410400 * u**28
    + 4068800 * u**26
    + 8413616 * u**24
    + 12776592 * u**22
    + 14351856 * u**20
    + 11894320 * u**18
    + 7327520 * u**16
    + 3706272 * u**14
    + 2123580 * u**12
    + 1571848 * u**10
    + 1038260 * u**8
    + 474252 * u**6
    + 134232 * u**4
    + 20604 * u**2
    + 1215
)
J = P / (800 * (1 + u * u) ** 18)
WEIGHT = s.diag(1, 10)


@cache
def data():
    th, j, c, e, d, h = s.symbols("theta J L ell delta H", real=True)
    n, v, vd, p, pd, gn, gv = s.symbols("n vhat vhatdot p pdot gN gV", real=True)
    eff = (
        -3 * vd**2
        + 6 * th * n * vd
        + (j - 3 * th**2) * n**2
        - c * n * v
        - s.Rational(9, 2) * e**2 * v**2
    )
    source = gn * n + gv * (v + d * n)
    momentum = s.diff(eff, vd)
    velocity = th * n - p / 6
    ham = (
        -(p**2) / 12
        + n * (th * p + c * v)
        - j * n**2
        + s.Rational(9, 2) * e**2 * v**2
        + source
    )
    lapse = (th * p + c * v + gn + d * gv) / (2 * j)
    reduced = (
        -(p**2) / 12
        + (th * p + c * v + gn + d * gv) ** 2 / (4 * j)
        + s.Rational(9, 2) * e**2 * v**2
        + gv * v
    )
    vflow = -p / 6 + th * (th * p + c * v + gn + d * gv) / (2 * j)
    pflow = (
        -3 * h * p - c * (th * p + c * v + gn + d * gv) / (2 * j) - 9 * e**2 * v - gv
    )
    A = s.Matrix(
        [
            [th * c / (2 * j), th**2 / (2 * j) - s.Rational(1, 6)],
            [-(c**2) / (2 * j) - 9 * e**2, -3 * h - th * c / (2 * j)],
        ]
    )
    B = s.Matrix(
        [[th / (2 * j), d * th / (2 * j)], [-c / (2 * j), -1 - d * c / (2 * j)]]
    )
    nr = s.Matrix([[c / (2 * j), th / (20 * j)]])
    ng = s.Matrix([[1 / (2 * j), d / (2 * j)]])
    checks = {
        "literal_momentum": momentum + 6 * vd - 6 * th * n,
        "literal_Legendre_with_physical_force_source": (p * vd - eff + source).subs(
            vd, velocity
        )
        - ham,
        "algebraic_constraint": s.diff(ham, n) - 2 * j * (lapse - n),
        "constraint_eliminated_Hamiltonian": ham.subs(n, lapse) - reduced,
        "weighted_Hamilton_velocity": s.diff(reduced, p) - vflow,
        "weighted_Hamilton_momentum": -s.diff(reduced, v) - 3 * h * p - pflow,
        "original_lapse_Euler_with_source": (
            s.diff(eff - source, n).subs(vd, velocity)
        ).subs(n, lapse),
        "original_scale_Euler_with_source": (
            -c * n - 9 * e**2 * v - pd - 3 * h * p - gv
        ).subs({n: lapse, pd: pflow}),
        "complete_physical_force_matrix": s.ImmutableMatrix(
            s.Matrix([vflow, pflow]) - A * s.Matrix([v, p]) - B * s.Matrix([gn, gv])
        ),
        "weighted_lapse_reconstruction": (
            nr * WEIGHT * s.Matrix([v, p]) + ng * s.Matrix([gn, gv])
        )[0]
        - lapse,
        "unforced_original_matter_rate": -w0 * n
        - 3 * ell * v
        + classical.quadratic()["w"] * n
        + 3 * ell * v,
        "actual_full_target_theta": theta - classical.quadratic()["theta"],
        "actual_full_target_J": J - classical.quadratic()["J"],
        "actual_full_target_matter_mix": w0 - classical.quadratic()["w"],
        "theta_alternative": theta - H + u / (1 + u * u) ** 4,
        "same_original_matter_background_charge": s.diff(ell, u) + 3 * H * ell,
    }
    sub = {th: theta, j: J, c: L, e: ell, d: delta, h: H}
    actualA = s.ImmutableMatrix(A.subs(sub))
    actualB = s.ImmutableMatrix(B.subs(sub))
    return {
        "original_classical_reference": "The unchanged S6.180 classical T0, not the new retuned classical action in isolation",
        "physical_forces": "gN,gV are original lapse and log-scale Euler forces divided by kappa*a0^3",
        "original_charge": "Original nonzero M1 charge retained; independent perturbation of that charge is zero",
        "reduced_forced_Hamiltonian": reduced,
        "generic_lapse_reconstruction": lapse,
        "actual_phase_matrix": actualA,
        "actual_physical_force_matrix": actualB,
        "weighted_phase_matrix": s.ImmutableMatrix(WEIGHT * actualA * WEIGHT.inv()),
        "weighted_force_matrix": s.ImmutableMatrix(WEIGHT * actualB),
        "actual_lapse_state_row": s.ImmutableMatrix(nr.subs(sub)),
        "actual_lapse_force_row": s.ImmutableMatrix(ng.subs(sub)),
        "J_positive_polynomial": P,
        "J": J,
        "generic_regular_denominator": 2 * j,
        "checks": {
            key: s.ImmutableMatrix(value.applyfunc(s.factor))
            if isinstance(value, s.MatrixBase)
            else s.factor(value)
            for key, value in checks.items()
        },
        "gates": {
            "all_J_polynomial_coefficients_positive": all(
                x > 0 for x in s.Poly(P, u).all_coeffs() if x != 0
            ),
            "J_polynomial_even": all(
                power[0] % 2 == 0 for power, _ in s.Poly(P, u).terms()
            ),
            "J_constant_positive": P.subs(u, 0) == 1215,
            "no_division_by_J_minus_three_theta_squared": s.denom(s.factor(lapse))
            == 2 * j,
        },
    }
