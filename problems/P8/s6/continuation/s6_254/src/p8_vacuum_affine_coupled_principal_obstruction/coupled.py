"""Whole metric, two-matter and longitudinal-vector quadratic reduction."""

from functools import cache

import sympy as s

D, Z, J, K, P, a, zeta, Y, Yv = s.symbols("D Z J K P a zeta Y Yv", positive=True)
r, rN, th, dH, H = s.symbols("r rN theta dH H", real=True)
L0, L2, C, Vvv = s.symbols("L0 L2 C Vvv", real=True)
cs = s.Matrix(s.symbols("c0:2", real=True))
ws = s.Matrix(s.symbols("w0:2", real=True))
ds = s.Matrix(s.symbols("d0:2", real=True))
vs = s.Matrix(s.symbols("v0:2", real=True))
m00, m01, m11 = s.symbols("mass00 mass01 mass11", real=True)
mass = s.Matrix([[m00, m01], [m01, m11]])
v, A, pv, pA, n, b, A0, vd, Ad = s.symbols("v A pv pA n b A0 vd Ad", real=True)
sigma = s.Matrix(s.symbols("sigma0:2", real=True))
sd = s.Matrix(s.symbols("sd0:2", real=True))
ps = s.Matrix(s.symbols("ps0:2", real=True))
PHASE = s.Matrix([v, *sigma, A, pv, *ps, pA])
RATES = s.Matrix([vd, *sd, Ad])
MOMENTA = s.Matrix([pv, *ps, pA])
q = P**2 / a**2
Lnv = L0 + L2 * q
GAMMA = 1 - 3 * Z * r**2 / (2 * D)
OMEGA = s.zeros(8)
OMEGA[:4, 4:] = s.eye(4)
OMEGA[4:, :4] = -s.eye(4)
Qb, QL, Pb, PL = s.symbols("Qb QL Pb PL", real=True)
CENTRAL = s.Matrix([Qb, *sigma, QL, Pb, *ps, PL])


def lagrangian():
    charge = cs.dot(sigma)
    source = r * (3 * vd - b) + 3 * rN * dH * n
    return (
        -3 * D * vd**2
        + Z * sd.dot(sd) / 2
        + 6 * th * n * vd
        + n * ws.dot(sd)
        + (J - 3 * th**2 / D + ws.dot(ws) / (2 * Z)) * n**2
        - 3 * vd * charge
        + Lnv * n * v
        + n * ds.dot(sigma)
        + Vvv * v**2
        + v * vs.dot(sigma)
        + (sigma.T * mass * sigma)[0] / 2
        + b * (2 * D * vd - 2 * th * n + charge)
        + 2 * C * q * v**2
        - Y * q * sigma.dot(sigma) / 2
        + K * (Ad - P * A0) ** 2 / 2
        + Z * (A0 + source) ** 2 / 2
        - Yv * A**2 / 2
    )


def lapse_numerator():
    return (
        th * (pv - 3 * r * P * pA + 3 * cs.dot(sigma)) / D
        - ws.dot(ps) / Z
        - Lnv * v
        - ds.dot(sigma)
        - 3 * rN * dH * P * pA
    )


def hamiltonian():
    charge = cs.dot(sigma)
    return (
        ps.dot(ps) / (2 * Z)
        - (pv - 3 * r * P * pA) * charge / (2 * D)
        - 3 * charge**2 / (4 * D)
        + Y * q * sigma.dot(sigma) / 2
        - (2 * C * q + Vvv) * v**2
        - v * vs.dot(sigma)
        - (sigma.T * mass * sigma)[0] / 2
        + lapse_numerator() ** 2 / (4 * J)
        + r * P * pA * pv / (2 * D)
        - 3 * r**2 * P**2 * pA**2 / (4 * D)
        + pA**2 / (2 * K)
        + P**2 * pA**2 / (2 * Z)
        + Yv * A**2 / 2
    )


@cache
def reduction():
    lag = lagrangian()
    kinetic = s.hessian(lag, RATES)
    linear = s.Matrix(
        [s.diff(lag, x).subs(dict.fromkeys(RATES, 0), simultaneous=True) for x in RATES]
    )
    velocity_solution = kinetic.inv() * (MOMENTA - linear)
    unred = s.expand(
        (MOMENTA.dot(RATES) - lag).subs(
            dict(zip(RATES, velocity_solution)), simultaneous=True
        )
    )
    aux = s.hessian(unred, [n, b, A0]).applyfunc(s.factor)
    A0star = s.factor(-(s.diff(unred, A0) - aux[2, 2] * A0) / aux[2, 2])
    solution = {n: lapse_numerator() / (2 * J), b: pv / (2 * D)}
    solution[A0] = A0star.subs(solution, simultaneous=True)
    checks = {
        "whole_four_velocity_Hessian": kinetic - s.diag(-6 * D * GAMMA, Z, Z, K),
        "whole_three_auxiliary_Hessian_determinant": aux.det()
        + 4 * D * Z * J / (3 * GAMMA),
        "whole_joint_configuration_measure_square": kinetic.det() * aux.det()
        - 8 * D**2 * Z**3 * J * K,
        "whole_shift_pivot": aux[1, 1] + 2 * D / 3,
        "whole_temporal_vector_pivot": aux[2, 2] + Z / GAMMA,
        "whole_last_lapse_pivot": aux[0, 0] - aux[0, 2] ** 2 / aux[2, 2] + 2 * J,
        "entire_eight_phase_reduced_Hamiltonian": unred.subs(
            solution, simultaneous=True
        )
        - hamiltonian(),
    }
    checks.update(
        {
            "complete_auxiliary_equation_" + str(x): s.diff(unred, x).subs(
                solution, simultaneous=True
            )
            for x in (n, b, A0)
        }
    )
    return {
        "whole_lagrangian": lag,
        "whole_unreduced_Hamiltonian": unred,
        "whole_reduced_Hamiltonian": hamiltonian(),
        "whole_velocity_solution": velocity_solution,
        "whole_auxiliary_solution": solution,
        "whole_auxiliary_Hessian": aux,
        "configuration_density_squared": 8 * D**2 * Z**3 * J * K,
        "checks": checks,
        "gates": {
            "both_matter_fields_and_all_source_contacts_retained": all(
                lag.has(x) for x in (*ds, *vs, m00, m01, m11)
            ),
            "no_heavy_reference_zero_substituted_off_clock": True,
            "joint_reduction_not_separate_vector_elimination": True,
        },
    }


def central_hamiltonian():
    charge = cs.dot(sigma)
    t = r * th / D + rN * dH
    linear = (
        -2 * a * th * P**2 * Qb / D
        + 3 * a**3 * th * charge / D
        - (L2 / 2 + L0 * a**2 / (2 * P**2)) * Pb
        - ws.dot(ps) / Z
        - a**3 * ds.dot(sigma)
        - 3 * s.sqrt(zeta) * t * PL
    )
    return (
        linear**2 / (4 * J * a**3)
        + ps.dot(ps) / (2 * Z * a**3)
        + charge * a * P**2 * Qb / D
        + 3 * charge * r * s.sqrt(zeta) * PL / (2 * D)
        - 3 * a**3 * charge**2 / (4 * D)
        + a * Y * P**2 * sigma.dot(sigma) / 2
        - C * Pb**2 / (2 * a * P**2)
        - a * Vvv * Pb**2 / (4 * P**4)
        - a**2 * vs.dot(sigma) * Pb / (2 * P**2)
        - a**3 * (sigma.T * mass * sigma)[0] / 2
        - r * s.sqrt(zeta) * P**2 * Qb * PL / (D * a**2)
        + zeta * GAMMA * PL**2 / (2 * Z * a**3)
        + zeta * PL**2 / (2 * K * P**2 * a**3)
        + a**3 * Yv * P**2 * QL**2 / (2 * zeta)
        - H * Qb * Pb
    )


@cache
def chart():
    normalization = s.diag(1, 1, 1, s.sqrt(zeta), a**3, a**3, a**3, a**3 / s.sqrt(zeta))
    physical = (
        a**3
        * normalization.inv().T
        * s.hessian(hamiltonian(), PHASE)
        * normalization.inv()
    )
    T = s.eye(8)
    T[0, 0] = T[4, 4] = 0
    T[0, 4], T[4, 0] = -1 / (2 * a * P**2), 2 * a * P**2
    T[3, 3], T[7, 7] = 1 / P, P
    connection = -OMEGA * H * a * T.diff(a) * T.inv()
    transformed = T.inv().T * physical * T.inv() + connection
    target = s.hessian(central_hamiltonian(), CENTRAL)
    return {
        "physical_normalization_without_common_sqrt_kappa": normalization,
        "whole_physical_Hessian": physical,
        "whole_central_symplectic_map": T,
        "whole_time_connection": connection,
        "whole_central_Hessian": target,
        "checks": {
            "whole_eight_phase_symplectic_map": T * OMEGA * T.T - OMEGA,
            "whole_physical_density_normalization": normalization
            * OMEGA
            * normalization.T
            - a**3 * OMEGA,
            "entire_time_dependent_central_Hamiltonian": transformed - target,
            "entire_central_connection": (CENTRAL.T * connection * CENTRAL)[0] / 2
            + H * Qb * Pb,
            "physical_density_leaves_full_determinant_unchanged": a**24
            / normalization.det() ** 2
            - 1,
        },
        "gates": {
            "comoving_P_strictly_positive": P.is_positive,
            "central_chart_does_not_divide_by_Theta": not T.has(th),
            "time_connection_not_deleted": connection != s.zeros(8),
        },
    }
