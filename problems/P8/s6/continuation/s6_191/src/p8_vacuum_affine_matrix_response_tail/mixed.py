"""Mixed time/amplitude derivative bounds without a root condition number."""

from functools import cache
from math import comb

import sympy as s
from p8_vacuum_affine_matrix_adiabatic import jets as prior

TIME = prior.JETS
PARAMETERS = 2


def subindices(index, exclude_ends=False):
    j, a = index
    for n in range(j + 1):
        for b in range(a + 1):
            if exclude_ends and (n, b) in ((0, 0), index):
                continue
            yield (n, b), (j - n, a - b), comb(j, n) * comb(a, b)


def ordered_indices(max_time=TIME, max_parameters=PARAMETERS):
    return sorted(
        ((j, a) for j in range(max_time + 1) for a in range(max_parameters + 1)),
        key=lambda x: (sum(x), x),
    )


@cache
def constants():
    zero = (0, 0)
    exponential = {}
    q, square = {}, {}
    for j, a in ordered_indices():
        value = sum(
            comb(j, l) * s.Integer(a) ** (j - l) * prior.bell(l) for l in range(j + 1)
        )
        exponential[j, a] = s.Integer(1) if (j, a) == zero else 2 * value
        if (j, a) == zero:
            q[j, a] = square[j, a] = s.Integer(1)
            continue
        first = sum(
            comb(j, l) * prior.scale_bound(1, j - l) * exponential[l, a]
            for l in range(j + 1)
        )
        q[j, a] = max(prior.scale_bound(3, j), first) if a == 0 else first
        square[j, a] = sum(
            comb(j, l) * prior.scale_bound(2, j - l) * exponential[l, a]
            for l in range(j + 1)
        )
    q[1, 0] = s.Integer(5)
    square[1, 0] = s.Integer(4)
    root, inverse, omega, reciprocal = ({zero: s.Integer(1)} for _ in range(4))
    for index in ordered_indices():
        if index == zero:
            continue
        root[index] = 2 * (
            q[index]
            + sum(
                f * root[left] * root[right]
                for left, right, f in subindices(index, True)
            )
        )
        inverse[index] = sum(
            f * inverse[left] * root[right]
            for left, right, f in subindices(index)
            if left != index
        )
        omega[index] = (
            square[index]
            + sum(
                f * omega[left] * omega[right]
                for left, right, f in subindices(index, True)
            )
        ) / 2
        reciprocal[index] = sum(
            f * reciprocal[left] * omega[right]
            for left, right, f in subindices(index)
            if left != index
        )
    L, rho = {}, {}
    for j, a in ordered_indices(TIME - 1):
        L[j, a] = sum(
            f * inverse[left] * root[right[0] + 1, right[1]]
            for left, right, f in subindices((j, a))
        )
        rho[j, a] = (
            sum(
                f * reciprocal[left] * omega[right[0] + 1, right[1]]
                for left, right, f in subindices((j, a))
            )
            / 2
        )
    return {
        "exp": exponential,
        "K": q,
        "omega_squared": square,
        "root": root,
        "inverse_root": inverse,
        "omega": omega,
        "inverse_omega": reciprocal,
        "L": L,
        "rho": rho,
        "R": L,
        "S": {index: L[index] + rho[index] for index in L},
    }


@cache
def data():
    c = constants()
    p = prior.constants()
    checks = {}
    mapping = {
        "exp": "matrix_exponential_relative",
        "K": "K_relative",
        "omega_squared": "omega_squared_relative",
        "root": "square_root_balanced",
        "inverse_root": "inverse_square_root_balanced",
        "omega": "omega_relative",
        "inverse_omega": "omega_reciprocal_relative",
        "L": "L_derivatives",
        "rho": "rho_derivatives",
        "R": "R_derivatives",
        "S": "S_derivatives",
    }
    for now, old in mapping.items():
        for j, value in enumerate(p[old]):
            checks[f"same_zero_parameter_{now}_{j}"] = c[now][j, 0] - value
    checks["first_relative_frequency_parameter_bound"] = c["omega"][0, 1] - 1
    checks["second_relative_frequency_parameter_bound"] = c["omega"][0, 2] - 2
    return {
        "domain": "gamma=epsilon Gamma, ||Gamma^(j)||<=1 through j12 and |epsilon|<=.01, with a common zero neighborhood of the initial Cauchy surface. Mixed distinct directions have the same product majorants.",
        "mixed_constants": c,
        "checks": checks,
        "gates": {
            "same_low_time_generator_bounds": c["R"][0, 0] == 10 and c["S"][0, 0] == 11,
            "all_mixed_constants_positive": all(
                v > 0 for row in c.values() for v in row.values()
            ),
            "parameter_derivatives_do_not_raise_time_jet_requirement": TIME >= 11,
        },
    }
