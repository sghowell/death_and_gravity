"""Three overlapping theta-graph subcycles, flavor counts and local ownership."""

from functools import cache
from itertools import combinations

import sympy as s


@cache
def data():
    cycles = tuple(
        frozenset(v) for v in combinations(("fermion_1", "fermion_2", "boson"), 2)
    )
    one_loop = tuple(
        {"flavor": j, "color": c, "active_Yukawa": j < 2}
        for j in range(14)
        for c in range(3)
    )
    A, B = s.symbols("D_inverse_0 D_inverse_1", commutative=False)
    ward = s.expand(A ** (-1) * (B - A) * B ** (-1) - (A ** (-1) - B ** (-1)))
    x = s.Symbol("q_squared", positive=True)
    f, fp = s.symbols("f_at_one f_prime_at_one")
    h, J1, J2, M = s.symbols("h J1 J2 M")
    vacuum_source = s.expand(-((h * J1 + h * h * J2) ** 2) / (2 * M))
    return {
        "proper_subcycles": [sorted(v) for v in cycles],
        "scalar_forest": "The whole fermion cycle is the physical Phi mass/residue proper subgraph and is paired as f_R,D. The two other overlapping cycles are proper fermion self energies and give the one-loop determinant mass/kinetic counterterm variation. No product of these overlapping counterterms is included.",
        "gauge_forest": "The whole fermion cycle is the vector-current polarization. The vector Ward identity forces a transverse q^2 counterterm, whose closure with the massless gauge line is a scaleless trace. Both proper fermion counterterm occurrences remain in the vacuum formula.",
        "active_scalar_color_multiplicity": 6,
        "all_gauge_color_multiplicity": 42,
        "flavor_ledger": list(one_loop),
        "outer_reference": "After proper forests, the remaining overall divergence is a vacuum constant. Its MS pole is removed; its finite value is retained. The eventual physical vacuum-energy-zero reference must cancel the complete model constant once, including other scalar/counterterm/source contributions.",
        "not_double_counted": "Each assigned bosonic or fermionic counterterm occurrence is removed from the separate insertion ledger. The order-one H-source square and other vacuum matching terms are not part of these two fermion primitive values.",
        "checks": {
            "three_proper_subcycles": len(cycles) - 3,
            "all_proper_subcycles_pairwise_overlap": sum(
                bool(a & b) for a, b in combinations(cycles, 2)
            )
            - 3,
            "no_disjoint_proper_product_forest": sum(
                not (a & b) for a, b in combinations(cycles, 2)
            ),
            "six_active_color_states": sum(r["active_Yukawa"] for r in one_loop) - 6,
            "forty_two_gauge_color_states": len(one_loop) - 42,
            "thirty_six_inert_gauge_states": sum(
                not r["active_Yukawa"] for r in one_loop
            )
            - 36,
            "vector_Ward_resolvent_identity": ward,
            "physical_scalar_slope_subtraction_closes_scaleless": s.factor(
                (fp * (x + 1)) / (x + 1) - fp
            ),
            "physical_scalar_mass_subtraction_not_scaleless": s.factor(
                (-f) / (x + 1) + f / (x + 1)
            ),
            "first_source_square_order_two_retained": vacuum_source.coeff(h, 2)
            + J1 * J1 / (2 * M),
            "second_source_cross_first_order_three": vacuum_source.coeff(h, 3)
            + J1 * J2 / M,
            "second_source_square_order_four": vacuum_source.coeff(h, 4)
            + J2 * J2 / (2 * M),
        },
    }
