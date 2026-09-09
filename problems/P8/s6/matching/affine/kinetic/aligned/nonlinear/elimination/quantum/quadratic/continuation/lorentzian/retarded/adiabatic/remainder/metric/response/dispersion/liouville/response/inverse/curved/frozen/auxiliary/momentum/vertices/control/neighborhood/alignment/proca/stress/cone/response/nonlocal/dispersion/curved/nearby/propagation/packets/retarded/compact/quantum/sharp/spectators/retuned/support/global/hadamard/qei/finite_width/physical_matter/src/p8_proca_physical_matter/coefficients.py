"""Fresh retuned constant-Proca lapse jets; no frozen source is modified."""

from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model as old
from p8_constant_proca import model as proca
from p8_proca_retuned_margin import model as retuned


@cache
def data():
    rows = dict(proca.jets())
    zero = (0,) * len(old.VARIABLES)
    c = old.coefficients()
    N = old.N
    delta = (
        -retuned.INCREMENT * N * c["U"] * (N**-2 - 1) ** 2 / c["background"]["h"] ** 2
    )
    increments = tuple(sp.factor(sp.diff(delta, N, j).subs(N, 1)) for j in range(5))
    before = rows[zero]
    rows[zero] = tuple(
        sp.factor(a + b) for a, b in zip(before, increments, strict=True)
    )
    Jnew = c["background"]["J"] + 4 * retuned.NEW_MARGIN / c["background"]["h"] ** 2
    return {
        "rows": rows,
        "old_background_row": before,
        "added_background_row": increments,
        "new_background_row": rows[zero],
        "Jnew": Jnew,
        "literal_retuned_Hamiltonian_increment": delta,
        "checks": {
            "retuned_actual_background_lapse_force_zero": rows[zero][1],
            "retuned_actual_nonzero_lapse_Hessian": sp.factor(rows[zero][2] + 2 * Jnew),
            "unchanged_background_value_of_retuning": increments[0],
            "unchanged_background_force_of_retuning": increments[1],
            "actual_covariant_margin_lapse_Hessian_increment": sp.factor(
                increments[2] + 8 * retuned.INCREMENT / c["background"]["h"] ** 2
            ),
            "retuned_ordinary_Proca_temporal_coefficient": sp.factor(
                rows[proca.powers(2, 2)][1] - sp.diff(N / (2 * c["U"]), N).subs(N, 1)
            ),
            "retuned_ordinary_Proca_spatial_mass_coefficient": sp.factor(
                rows[proca.powers(6, 1)][1]
                - sp.diff(N * c["e_omega"] / 2, N).subs(N, 1)
            ),
        },
    }
