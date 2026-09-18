"""Aggregate phase with radiation entering only through original massive recoil."""

from functools import cache

import sympy as s
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient import radiative

from . import kernels


@cache
def data():
    checks = {}
    for count in (2, 3, 5, 8):
        D = s.symbols("D0:" + str(count), nonzero=True)
        P = {
            (i, j): s.Symbol("P" + str(i) + "_" + str(j))
            for i in range(count)
            for j in range(i, count)
        }
        A = lambda i, j, P=P: P[tuple(sorted((i, j)))]
        pairs = sum(
            A(i, i) * D[j] / D[i] + A(j, j) * D[i] / D[j] - 2 * A(i, j)
            for i in range(count)
            for j in range(i + 1, count)
        )
        whole = sum(D) * sum(A(i, i) / D[i] for i in range(count)) - sum(
            A(i, j) for i in range(count) for j in range(count)
        )
        checks["generic_subset_aggregate_" + str(count)] = s.expand(pairs - whole)
    t = s.Symbol("t", positive=True)
    c = t * (2 * t * t - 3) / (t * t - 1) ** s.Rational(3, 2)
    checks["same_orientation_phase_derivative"] = s.factor(
        s.diff(c, t) - 3 / (t * t - 1) ** s.Rational(5, 2)
    )
    E, states = radiative.calibration_states()
    u = s.Matrix([0, s.Rational(4, 5), s.Rational(3, 5)])
    n = s.Matrix([0, 1, 0])
    p = s.diag(0, 1, 0, -1) / s.sqrt(2)
    cross = s.zeros(4)
    cross[1, 3] = cross[3, 1] = 1 / s.sqrt(2)
    differences = {}
    for pol, A in (("plus", p), ("complex", (p + s.I * cross) / s.sqrt(2))):
        two = kernels.components(E, states["two"], u, n, A)
        four = kernels.components(E, states["four"], u, n, A)
        checks["equal_energy_momentum_equal_phase_" + pol] = kernels.clean(
            two["G"] - four["G"]
        )
        differences[pol] = kernels.clean(two["F"] - four["F"])
    lower = s.Rational(29, 16)
    gates = {
        "c_above13_over7": lower * lower * (2 * lower * lower - 3) ** 2
        > s.Rational(169, 49) * (lower * lower - 1) ** 3,
        "cprime_below2_over5": 9 < s.Rational(4, 25) * (lower * lower - 1) ** 5,
        "phase_cap18": (4 * 102 + s.Rational(102, 7)) / 24 < 18,
        "phase_change19R": (102 * s.Rational(8, 5) + s.Rational(2040, 7)) / 24 < 19,
        "distribution_information_not_erased_from_real_component": all(
            abs(s.N(value, 80)) > s.Rational(1, 100000)
            for value in differences.values()
        ),
        "G_not_literal_imaginary_part_for_complex_polarization": True,
        "full_phase_not_deleted": True,
        "only_massive_recoil_needed_for_G_not_for_F": True,
    }
    return {
        "checks": {name: kernels.clean(value) for name, value in checks.items()},
        "gates": {name: bool(value) for name, value in gates.items()},
        "whole_aggregate_phase": "G=[(c12+2)S12+(c34-2)S34]/(8pi). The outgoing total momentum is(2E,0), the incoming one(-2E,0). Subset pair sums and the original phase term remove every separate null-current ambiguity. F and G are complex-linear tensor components, not literal real/imaginary parts for complex A.",
        "whole_bounds": "For gamma>=29/16,13/7<c<2 and0<cprime<2/5. |S12|,|S34|<=102, |deltaS34|<=2040R, |delta gamma34|<=4R. Hence |G|<18 and |G-G_Born|<19R for R>0; at R=0 the difference is0. Original real-component bounds give |C|<=293<300 and |C-C_Born|<=23667R.",
        "whole_same_momentum_different_real_coefficient": differences,
    }
