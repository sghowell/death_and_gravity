"""Parameter derivatives of the full finite ordered reference and residual."""

from functools import cache

import sympy as s
from p8_vacuum_affine_matrix_adiabatic import jets as old_jets
from p8_vacuum_affine_matrix_adiabatic import riccati as old_reference

from . import mixed

ORDER = old_jets.ORDER
MASS = old_jets.MASS


def triples(index):
    for left, rest, f in mixed.subindices(index):
        for middle, right, g in mixed.subindices(rest):
            yield left, middle, right, f * g


def parameter_reference(omega, rotation, squeeze, order=ORDER):
    if type(order) is not int or not 1 <= order <= ORDER:
        raise ValueError("Require a native finite reference order1..10")
    needed = {(j, a) for j in range(order + 1) for a in range(3)}
    if any(not needed.issubset(row) for row in (omega, rotation, squeeze)):
        raise ValueError("Require complete mixed time/parameter jets")
    zero_index = (0, 0)
    inverse = {zero_index: 1 / omega[zero_index]}
    for index in mixed.ordered_indices(order):
        if index == zero_index:
            continue
        inverse[index] = -inverse[zero_index] * sum(
            f * omega[left] * inverse[right]
            for left, right, f in mixed.subindices(index)
            if left != zero_index
        )
    zero = s.zeros(squeeze[zero_index].rows)
    b = {}
    for d in range(order + 1):
        for a in range(3):
            b[1, d, a] = (
                s.I
                * sum(
                    (
                        f * inverse[left] * squeeze[right]
                        for left, right, f in mixed.subindices((d, a))
                    ),
                    zero,
                )
                / 2
            ).applyfunc(s.expand)
    for n in range(1, order):

        def right(d, a, n=n):
            value = b[n, d + 1, a] - sum(
                (
                    f
                    * (
                        rotation[left] * b[n, rest[0], rest[1]]
                        - b[n, rest[0], rest[1]] * rotation[left]
                    )
                    for left, rest, f in mixed.subindices((d, a))
                ),
                zero,
            )
            for j in range(1, n):
                l = n - j
                value += sum(
                    (
                        f
                        * b[j, left[0], left[1]]
                        * squeeze[middle]
                        * b[l, last[0], last[1]]
                        for left, middle, last, f in triples((d, a))
                    ),
                    zero,
                )
            return value

        for d in range(order - n + 1):
            for a in range(3):
                b[n + 1, d, a] = (
                    -s.I
                    * sum(
                        (
                            f * inverse[left] * right(rest[0], rest[1])
                            for left, rest, f in mixed.subindices((d, a))
                        ),
                        zero,
                    )
                    / 2
                ).applyfunc(s.expand)
    return b


@cache
def constants():
    c = mixed.constants()
    v, R, S = c["inverse_omega"], c["R"], c["S"]
    b = {}
    for d in range(ORDER + 1):
        for a in range(3):
            b[1, d, a] = (
                sum(
                    f * v[left] * S[right]
                    for left, right, f in mixed.subindices((d, a))
                )
                / 2
            )
    for n in range(1, ORDER):

        def right(d, a, n=n):
            value = b[n, d + 1, a] + 2 * sum(
                f * R[left] * b[n, right[0], right[1]]
                for left, right, f in mixed.subindices((d, a))
            )
            for j in range(1, n):
                l = n - j
                value += sum(
                    f * b[j, left[0], left[1]] * S[middle] * b[l, last[0], last[1]]
                    for left, middle, last, f in triples((d, a))
                )
            return value

        for d in range(ORDER - n + 1):
            for a in range(3):
                b[n + 1, d, a] = (
                    sum(
                        f * v[left] * right(rest[0], rest[1])
                        for left, rest, f in mixed.subindices((d, a))
                    )
                    / 2
                )
    size = {}
    residual = {}
    for a in range(3):
        size[a] = sum(b[n, 0, a] / MASS ** (n - 1) for n in range(1, ORDER + 1))
        value = b[ORDER, 1, a] + 2 * sum(
            f * R[left] * b[ORDER, right[0], right[1]]
            for left, right, f in mixed.subindices((0, a))
        )
        for j in range(1, ORDER + 1):
            for l in range(1, ORDER + 1):
                if j + l < ORDER:
                    continue
                value += sum(
                    f * b[j, left[0], left[1]] * S[middle] * b[l, last[0], last[1]]
                    for left, middle, last, f in triples((0, a))
                ) / MASS ** (j + l - ORDER)
        residual[a] = value
    return {
        "ordered_mixed_coefficient_bounds": b,
        "reference_parameter_norm_over_inverse_frequency": size,
        "residual_parameter_norm_over_inverse_frequency_tenth": residual,
    }


@cache
def data():
    c = constants()
    old = old_reference.constants()
    checks = {}
    for n, row in old["coefficient_derivative_bounds"].items():
        for d, value in enumerate(row):
            checks[f"same_undifferentiated_reference_{n}_{d}"] = (
                c["ordered_mixed_coefficient_bounds"][n, d, 0] - value
            )
    checks["same_exact_undifferentiated_residual"] = (
        c["residual_parameter_norm_over_inverse_frequency_tenth"][0]
        - old["residual_over_inverse_frequency_tenth"]
    )
    checks["same_exact_undifferentiated_reference_size"] = (
        c["reference_parameter_norm_over_inverse_frequency"][0]
        - MASS * old["reference_norm_at_mass_floor"]
    )
    return {
        "reference": "Differentiate the full ordered S6.190 recurrence through two amplitude parameters, retaining all mixed time jets and inverse-frequency derivatives. No commuting polarization substitution is made.",
        "constants": c,
        "bounds": "||partial_epsilon^a rhat||<=A_a nu_minus^-1 and ||partial_epsilon^a F||<=C_a nu_minus^-10 for a=0,1,2. Parameter differentiation preserves reference frequency order.",
        "checks": checks,
        "gates": {
            "all_mixed_reference_constants_positive": all(
                x > 0 for row in c.values() for x in row.values()
            ),
            "same_reference_smallness": c[
                "reference_parameter_norm_over_inverse_frequency"
            ][0]
            / MASS
            < s.Rational(1, 100),
            "no_extra_time_jet_demand": ORDER + 1 <= mixed.TIME,
        },
    }
