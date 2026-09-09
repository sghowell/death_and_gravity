"""Global Cartesian polarization projectors; no singular angular frame."""

from functools import cache

import sympy as sp


@cache
def data():
    k = sp.Matrix(sp.symbols("kx ky kz", real=True))
    q = k.dot(k)
    PL = k * k.T / q
    PT = sp.eye(3) - PL
    Pi = sp.Matrix(
        9,
        9,
        lambda row, col: (
            (
                PT[row // 3, col // 3] * PT[row % 3, col % 3]
                + PT[row // 3, col % 3] * PT[row % 3, col // 3]
                - PT[row // 3, row % 3] * PT[col // 3, col % 3]
            )
            / 2
        ),
    )
    clean = lambda M: M.applyfunc(sp.factor)
    trace = sp.Matrix(1, 9, lambda _, j: int(j // 3 == j % 3))
    div = sp.Matrix(3, 9, lambda i, j: k[j % 3] if j // 3 == i else 0)
    E1 = sp.diag(0, 1, -1)
    E2 = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]])
    e1, e2 = (sp.Matrix(list(E)) for E in (E1, E2))
    aligned = Pi.subs({k[0]: 1, k[1]: 0, k[2]: 0})
    return {
        "spatial_momentum": k,
        "longitudinal": PL,
        "transverse": PT,
        "TT": Pi,
        "checks": {
            "Cartesian_longitudinal_projector": clean(PL * PL - PL),
            "Cartesian_transverse_projector": clean(PT * PT - PT),
            "complementary_Cartesian_Proca_polarizations": clean(PL * PT),
            "longitudinal_polarization_rank": sp.factor(sp.trace(PL) - 1),
            "transverse_polarization_rank": sp.factor(sp.trace(PT) - 2),
            "Cartesian_TT_projector": clean(Pi * Pi - Pi),
            "Cartesian_TT_projector_selfadjoint": clean(Pi - Pi.T),
            "Cartesian_TT_trace_constraint": clean(trace * Pi),
            "Cartesian_TT_spatial_constraint": clean(div * Pi),
            "two_physical_tensor_polarizations": sp.factor(sp.trace(Pi) - 2),
            "actual_E_norm_two_polarization_completion": clean(
                aligned - (e1 * e1.T + e2 * e2.T) / 2
            ),
        },
    }
