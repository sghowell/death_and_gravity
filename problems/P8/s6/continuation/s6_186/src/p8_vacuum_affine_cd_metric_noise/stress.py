"""Literal physical Proca stress and all ten constrained mode readouts."""

from functools import cache

import sympy as s
from p8_vacuum_affine_proca_gaussian import bridge


@cache
def data():
    E = s.Matrix(s.symbols("E0:3", real=True))
    B = s.Matrix(s.symbols("B0:3", real=True))
    massA = s.Matrix(s.symbols("mA0:4", real=True))
    fields = s.Matrix([*E, *B, *massA])
    eta = s.diag(1, -1, -1, -1)
    F = s.zeros(4)
    for i in range(3):
        F[0, i + 1] = E[i]
        F[i + 1, 0] = -E[i]
        for j in range(3):
            F[i + 1, j + 1] = sum(s.LeviCivita(i, j, k) * B[k] for k in range(3))
    square = sum(
        F[a, b] * eta[a, a] * eta[b, b] * F[a, b] for a in range(4) for b in range(4)
    )
    mass_square = (massA.T * eta * massA)[0]
    T = s.Matrix(
        4,
        4,
        lambda a, b: s.expand(
            -sum(F[a, c] * eta[c, c] * F[b, c] for c in range(4))
            + eta[a, b] * square / 4
            + massA[a] * massA[b]
            - eta[a, b] * mass_square / 2
        ),
    )
    target = s.zeros(4)
    target[0, 0] = (E.dot(E) + B.dot(B) + massA.dot(massA)) / 2
    flux = E.cross(B) + massA[0] * s.Matrix(massA[1:])
    for i in range(3):
        target[0, i + 1] = target[i + 1, 0] = flux[i]
        for j in range(3):
            target[i + 1, j + 1] = (
                -E[i] * E[j] - B[i] * B[j] + massA[i + 1] * massA[j + 1]
            )
            if i == j:
                target[i + 1, j + 1] += (
                    E.dot(E)
                    + B.dot(B)
                    + massA[0] ** 2
                    - sum(massA[k] ** 2 for k in range(1, 4))
                ) / 2
    matrices = {
        (a, b): s.hessian(T[a, b], fields) / 2 for a in range(4) for b in range(4)
    }
    checks = {"literal_covariant_full_stress": s.simplify(T - target)}
    for (a, b), M in matrices.items():
        checks[f"stress_quadratic_form_{a}_{b}"] = s.expand(
            (fields.T * M * fields)[0] - T[a, b]
        )
    for a in range(4):
        for b in range(4):
            M = matrices[a, b]
            assert M == M.T
            assert max(
                sum(abs(M[i, j]) for j in range(10)) for i in range(10)
            ) <= s.Rational(1, 2)
    scale, m, k, omega = s.symbols("a m k omega", positive=True)
    f, p = s.symbols("f p")
    transverse = s.Matrix([p, 0, 0, 0, s.I * k / scale * f, 0, 0, m * f, 0, 0])
    longitudinal = s.Matrix(
        [0, 0, m / omega * p, 0, 0, 0, -s.I * k / (scale * omega) * p, 0, 0, omega * f]
    )
    # Pick k along the third axis and one transverse polarization along the first.
    # A common a^(-3/2) factor is outside these phase-space readouts.
    fr, fi, pr, pi = s.symbols("f_real f_imag p_real p_imag", real=True)
    complex_readout = {f: fr + s.I * fi, p: pr + s.I * pi}
    Tc = transverse.subs(complex_readout)
    Lc = longitudinal.subs(complex_readout)
    energyT = (pr**2 + pi**2 + (m * m + k * k / scale**2) * (fr**2 + fi**2)) / 2
    energyL = (pr**2 + pi**2 + omega**2 * (fr**2 + fi**2)) / 2
    Tnorm = s.expand((Tc.conjugate().T * Tc)[0] / 2)
    Lnorm = s.expand((Lc.conjugate().T * Lc)[0] / 2)
    checks["transverse_physical_energy"] = s.expand(Tnorm - energyT)
    checks["longitudinal_physical_energy"] = s.factor(
        (Lnorm - energyL).subs(omega, s.sqrt(m * m + k * k / scale**2))
    )
    checks["same_mass"] = bridge.MASS - 1000
    gT = s.sqrt(scale)
    gL = s.sqrt(scale) * m / omega
    V = scale ** s.Rational(3, 2)
    checks.update(
        {
            "transverse_electric_constraint_readout": s.simplify(
                V * gT * p / scale**2 - transverse[0]
            ),
            "transverse_magnetic_constraint_readout": s.simplify(
                V * s.I * k * f / (scale**2 * gT) - transverse[4]
            ),
            "transverse_mass_spatial_readout": s.simplify(
                V * m * f / (scale * gT) - transverse[7]
            ),
            "longitudinal_electric_constraint_readout": s.simplify(
                V * gL * p / scale**2 - longitudinal[2]
            ),
            "longitudinal_temporal_constraint_readout": s.simplify(
                -V * s.I * k * gL * p / (scale**3 * m) - longitudinal[6]
            ),
            "longitudinal_mass_spatial_readout": s.simplify(
                V * m * f / (scale * gL) - longitudinal[9]
            ),
        }
    )
    checks["exact_volume_cancellation"] = (
        scale**3 * scale ** s.Rational(-3, 2) * scale ** s.Rational(-3, 2) - 1
    )
    return {
        "field_order": list(map(str, fields)),
        "physical_stress": T,
        "quadratic_matrices": matrices,
        "each_stress_matrix_operator_norm_upper": s.Rational(1, 2),
        "contracted_matrix_bound": "||M_f||op <= (1/2) sum_(a,b)|f_ab| <=2 ||f||F; all16 tensor entries",
        "phase_stripped_mode_readouts": {
            "one_transverse": transverse,
            "longitudinal": longitudinal,
        },
        "common_volume_factor": "a^(-3/2) in every field mode; a^3 in the local stress smearing",
        "checks": checks,
        "gates": {
            "all_sixteen_stress_matrices_symmetric": all(
                M == M.T for M in matrices.values()
            ),
            "all_sixteen_operator_norms_at_most_half": all(
                max(sum(abs(M[i, j]) for j in range(10)) for i in range(10))
                <= s.Rational(1, 2)
                for M in matrices.values()
            ),
            "complete_positive_energy_matrix": matrices[0, 0] == s.eye(10) / 2,
        },
    }
