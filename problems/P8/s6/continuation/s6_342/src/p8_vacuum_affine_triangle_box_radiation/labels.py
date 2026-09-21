"""All24 labelled kernels and the six distinct mass-ordered flat boxes."""

import itertools
from functools import cache
from itertools import permutations

import sympy as s
from p8_vacuum_affine_heavy_scalar_four_point_loop import amplitude as old

from . import source
from .radiation import dot

PERMUTATIONS = tuple(permutations(range(4)))
C_OFF = s.Function("Cbar_offshell")
D_OFF = s.Function("Dbar_ordered_offshell")


def endpoint(value):
    return source.CONTACT + source.CUBIC**2 / (source.HEAVY_MASS2 - value)


def hard(momenta, triangle_value=C_OFF, box_value=D_OFF):
    if len(momenta) != 4 or sum(momenta, s.zeros(4, 1)) != s.zeros(4, 1):
        raise ValueError("Require a conserved four-point hard loop configuration")
    result = 0
    for a, b, c, d in PERMUTATIONS:
        v = dot(momenta[a] + momenta[b], momenta[a] + momenta[b])
        w = dot(momenta[b] + momenta[c], momenta[b] + momenta[c])
        virtualities = tuple(dot(momenta[j], momenta[j]) for j in (a, b, c, d))
        result += (
            source.CUBIC**2
            * endpoint(v)
            * triangle_value(v, virtualities[3], virtualities[2])
            / 4
        )
        result += source.CUBIC**4 * box_value(v, w, *virtualities) / 4
    return result


@cache
def data():
    checks = {}
    gates = {}
    S, T, U = old.svar, old.tvar, old.uvar
    pairs = {
        frozenset((0, 1)): S,
        frozenset((2, 3)): S,
        frozenset((0, 2)): T,
        frozenset((1, 3)): T,
        frozenset((0, 3)): U,
        frozenset((1, 2)): U,
    }
    tri_flat = 0
    box_flat = 0
    counts = {}
    for a, b, c, d in itertools.permutations(range(4)):
        v = pairs[frozenset((a, b))]
        w = pairs[frozenset((b, c))]
        tri_flat += old.g**2 * old.fixed_vertex(v) * old.Cfun(v) / 4
        box_flat += old.g**4 * old.Dfun(v, w) / 4
        counts[(v, w)] = counts.get((v, w), 0) + 1
    checks["all24_triangles_exact_S235_grouping"] = s.expand(
        tri_flat
        - sum(2 * old.g**2 * old.fixed_vertex(v) * old.Cfun(v) for v in (S, T, U))
    )
    checks["all24_boxes_exact_S235_ordered_grouping"] = s.expand(
        box_flat
        - sum(
            old.g**4 * old.Dfun(v, w) for v, w in itertools.permutations((S, T, U), 2)
        )
    )
    checks["full_S235_minus_bubbles_is_exact_labelled_TD"] = s.expand(
        16 * s.pi**2 * old.complete_loop()
        - sum(old.fixed_vertex(v) ** 2 * old.Bfun(v) / 2 for v in (S, T, U))
        - tri_flat
        - box_flat
    )
    gates["six_ordered_boxes_each_four_labels"] = len(counts) == 6 and set(
        counts.values()
    ) == {4}
    gates["mass_assignment_swap_is_not_generic_identity"] = (
        old.Dfun(S, T) - old.Dfun(T, S) != 0
    )

    return {
        "checks": {key: s.expand(value) for key, value in checks.items()},
        "gates": {key: bool(value) for key, value in gates.items()},
        "whole_label_count": len(PERMUTATIONS),
        "whole_ordered_box_multiplicities": counts,
        "whole_flat_triangle_sum": tri_flat,
        "whole_flat_ordered_box_sum": box_flat,
        "whole_full_offshell_definition": "F_T=g^2/4 sum_perm A(v_ab) Cbar(v_ab;p_d^2,p_c^2); F_D=g^4/4 sum_perm Dbar(v_ab,v_bc;p_a^2,p_b^2,p_c^2,p_d^2). Every master is the complete original-mass S341 Feynman contour integral, not a Taylor truncation. Restore1/(16pi^2). Cbar's two light legs are exchange-symmetric; ordered Dbar's two opposite mass cuts are not.",
        "whole_fixed_symmetric_subtraction": "Subtract F_sym=F_T+F_D at all scalar virtualities1 and s=t=u=4/3 once. This is the selected sector's linear contribution to the existing fixed OS4 condition; neither the original Born contact nor a curvature coefficient is retuned.",
    }
