"""Necessary finite-order support, not a count of connected mixed diagrams."""

from functools import cache

import sympy as s
from p8_vacuum_quantum_map import interactions as first


def max_map_degree(external, loops):
    if type(external) is not int or external not in (0, 2, 4):
        raise ValueError("Need native physical Phi source count zero, two or four")
    if type(loops) is not int or not 0 <= loops <= 2:
        raise ValueError("Need native loop order zero through two")
    return max(0, loops - 1 + external // 2)


def support_rows():
    # (name, Phi slots, other slots). Free quadratic p=0 is not an
    # interaction at j=0; its generated p>0 terms are interactions.
    classes = (
        ("Phi_quadratic", 2, 0),
        ("Phi_quartic", 4, 0),
        ("H_Phi_squared", 2, 1),
        ("Yukawa", 1, 2),
        ("H_quadratic", 0, 2),
        ("fermion_quadratic", 0, 2),
        ("gauge_quadratic", 0, 2),
        ("gauge_fermion", 0, 3),
        ("gauge_cubic", 0, 3),
        ("gauge_quartic", 0, 4),
        ("ghost_gauge", 0, 3),
    )
    rows = []
    for name, phi, other in classes:
        for j in range(3):
            for p in range(phi + 1):
                generated = phi + other + 2 * p
                if j == 0 and phi + other == 2 and p == 0:
                    continue
                if generated + 2 * j <= 8:
                    rows.append(
                        {
                            "operator": name,
                            "parent_valence": phi + other,
                            "map_degree": p,
                            "generated_valence": generated,
                            "counterterm_loop_weight": j,
                        }
                    )
    return rows


def validate_support(rows):
    if rows != support_rows():
        raise ValueError("The finite-order generated-operator support changed")
    return True


@cache
def data():
    E, N, V, J, P = s.symbols("E total_parent_valence V counterterm_weight map_degree")
    # Includes physical source endpoints as vertices; P includes their
    # cubic replacements as well as action replacements.
    edges = (N + E + 2 * P) / 2
    loops = edges - V - E + 1 + J
    formula = 1 - E / 2 + (N - 2 * V) / 2 + P + J
    return {
        "effective_loop_degree_formula": formula,
        "physical_source_map_degree_caps": {
            f"E{e}_L{l}": max_map_degree(e, l) for e in (0, 2, 4) for l in range(3)
        },
        "four_light_two_loop_scalar_topologies": first.topologies(4, 2),
        "necessary_mixed_operator_support": support_rows(),
        "scope": "The condition generated_valence+2j<=8 is necessary support for four physical light sources through two loops; not every listed tuple is a connected mixed diagram. The full action, not this support list, defines the theory. Field-independent vacuum and linear H-source counterterms have map degree zero and are kept separately.",
        "checks": {
            "Euler_and_physical_source_halfedges": s.expand(loops - formula),
            "four_sources_two_loops_need_map_order_three": max_map_degree(4, 2) - 3,
            "two_sources_two_loops_need_map_order_two": max_map_degree(2, 2) - 2,
            "vacuum_two_loops_need_map_order_one": max_map_degree(0, 2) - 1,
            "pure_scalar_two_loop_patterns": sum(
                v not in ((3, 0, 0, 0, 0), (1, 1, 0, 0, 0), (0, 0, 1, 0, 0))
                for v in first.topologies(4, 2)
            ),
            "three_scalar_two_loop_patterns": len(first.topologies(4, 2)) - 3,
            "generated_octic_is_required": sum(
                r["operator"] == "Phi_quartic"
                and r["map_degree"] == 2
                and r["counterterm_loop_weight"] == 0
                for r in support_rows()
            )
            - 1,
            "first_quadratic_counterterm_sextic_retained": sum(
                r["operator"] == "Phi_quadratic"
                and r["map_degree"] == 2
                and r["counterterm_loop_weight"] == 1
                for r in support_rows()
            )
            - 1,
            "first_Yukawa_counterterm_cubic_Phi_retained": sum(
                r["operator"] == "Yukawa"
                and r["map_degree"] == 1
                and r["counterterm_loop_weight"] == 1
                for r in support_rows()
            )
            - 1,
            "all_support_weighted_valences_at_most_eight": sum(
                r["generated_valence"] + 2 * r["counterterm_loop_weight"] > 8
                for r in support_rows()
            ),
        },
    }
