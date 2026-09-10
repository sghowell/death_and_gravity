"""Complete two-point one-loop counterterm insertions and outer on-shell grouping."""

from functools import cache

import sympy as sp
from p8_vacuum_forward_loop import subtraction as inherited

from . import denominators, graphs


def on_shell(expr, s):
    return expr - expr.subs(s, 1) - (s - 1) * sp.diff(expr, s).subs(s, 1)


@cache
def data():
    L, g, M, I0 = sp.symbols(
        "quartic_L cubic_squared heavy_mass_squared regulated_zero_bubble", real=True
    )
    s = sp.Symbol("Minkowski_invariant")
    h = sp.symbols("first_heavy_inverse second_heavy_inverse third_heavy_inverse")
    T, B, B2 = sp.symbols(
        "regulated_light_tadpole mixed_one_loop_bubble mixed_one_loop_squared_heavy"
    )
    alpha, R, Q = sp.symbols("fixed_alpha decaying_R full_Q")
    finite = sp.Symbol("once_fixed_finite_quartic_contact")
    d = inherited.data()
    dL, dg, dM = [
        d[k]
        for k in (
            "delta_polynomial_quartic_UV",
            "delta_cubic_squared_UV",
            "delta_heavy_mass_squared_UV",
        )
    ]
    Ibar = I0 / (16 * sp.pi**2)
    one_loop = (-L + g / M) * T / 2 + g * B
    local_ct = (-dL + dg / M - g * dM / M**2) * T / 2
    p_ct = dg * B - g * dM * B2
    bare_proper_uv = -L * g * Ibar * B + g**2 * Ibar * B2 / 2
    full_square = (-L + g * sum(h)) ** 2 / 6
    local = L**2 / 6
    mixed = -L * g * sum(h) / 3
    diagonal = g**2 * sum(v * v for v in h) / 6
    different = g**2 * sum(h[i] * h[j] for i in range(3) for j in range(i + 1, 3)) / 3
    hbar, v, H1, H2 = sp.symbols(
        "loop_order light_shell_distance one_loop_analytic_H two_loop_analytic_H"
    )
    inverse = v + hbar * v**2 * H1 + hbar**2 * v**2 * H2
    c0, c1 = sp.symbols("arbitrary_local_mass arbitrary_local_kinetic")
    f = sp.Function("regulated_self_energy")(s)
    checks = {
        "same_entire_regulated_reference_cubic_counterterm": sp.expand(
            dg - L * g * Ibar
        ),
        "same_entire_regulated_reference_heavy_mass_counterterm": sp.expand(
            dM - g * Ibar / 2
        ),
        "complete_one_loop_parameter_variation": sp.expand(
            dL * sp.diff(one_loop, L)
            + dg * sp.diff(one_loop, g)
            + dM * (sp.diff(one_loop, M) - g * B2)
            - local_ct
            - p_ct
        ),
        "all_momentum_dependent_proper_UV_terms_cancel": sp.expand(
            p_ct + bare_proper_uv
        ),
        "mixed_proper_subtraction": sp.expand(-L * g * B * Ibar + dg * B),
        "same_heavy_proper_subtraction": sp.expand(g**2 * B2 * Ibar / 2 - g * dM * B2),
        "full_sunset_vertex_square_partition": sp.expand(
            full_square - local - mixed - diagonal - different
        ),
        "summed_mixed_sunset_weight": sp.expand(
            mixed.subs(dict.fromkeys(h, 1)) + L * g
        ),
        "summed_same_heavy_sunset_weight": sp.expand(
            diagonal.subs(dict.fromkeys(h, 1)) - g**2 / 2
        ),
        "summed_different_heavy_sunset_weight": sp.expand(
            different.subs(dict.fromkeys(h, 1)) - g**2
        ),
        "outer_OS_annihilates_all_local_parameter_CT": on_shell(local_ct, s),
        "finite_quartic_insertion_is_p_independent": sp.diff(
            finite * sp.diff(one_loop, L), s
        ),
        "finite_quartic_insertion_outer_OS_zero": on_shell(
            finite * sp.diff(one_loop, L), s
        ),
        "arbitrary_local_mass_and_kinetic_OS_zero": sp.expand(on_shell(c0 + c1 * s, s)),
        "nested_constant_alpha_multiplies_full_OS_one_loop": sp.expand(
            on_shell(alpha * f, s) - alpha * on_shell(f, s)
        ),
        "nested_Q_split": sp.expand(
            (
                ((-L + g / M) + 2 * g * h[0]) * Q / 2
                - (((-L + g / M) + 2 * g * h[0]) * (alpha + R) / 2)
            ).subs(Q, alpha + R)
        ),
        "OS_self_energy_value": sp.simplify(on_shell(f, s).subs(s, 1)),
        "OS_self_energy_derivative": sp.simplify(sp.diff(on_shell(f, s), s).subs(s, 1)),
        "two_loop_inverse_factorization": sp.expand(
            inverse - v * (1 + v * (hbar * H1 + hbar**2 * H2))
        ),
        "two_loop_unit_light_residue": sp.limit(v / inverse, v, 0) - 1,
    }
    for kind, choice in graphs.cases():
        gph = graphs.refine(kind, choice)
        label = graphs.group(kind, choice)
        key = kind + "_" + "".join(map(str, choice))
        if label == "outer_constant":
            checks[key + "_all_external_momentum_absent"] = denominators.point(
                kind, choice
            )["P"]
        if kind == "sunset" and choice == (0, 0):
            source = next(n for n in gph["vertices"] if 0 in gph["external"][n])
            sink = next(n for n in gph["vertices"] if 1 in gph["external"][n])
            for index, row in enumerate(graphs.uv(gph)):
                if row[1] == 1:
                    vertices = {n for j in row[0] for n in gph["edges"][j]}
                    checks[
                        "local_sunset_contracted_proper_core_"
                        + str(index)
                        + "_external_coincidence"
                    ] = int(source in vertices and sink in vertices) - 1
    return {
        "same_entire_regulated_zero_bubble": I0,
        "one_loop_two_point_form": one_loop,
        "complete_local_parameter_counterterm_piece": local_ct,
        "complete_momentum_dependent_parameter_counterterm_piece": p_ct,
        "sunset_proper_UV_piece_before_subtraction": bare_proper_uv,
        "sunset_vertex_product": full_square,
        "sunset_partition": {
            "local": local,
            "mixed": mixed,
            "same_heavy": diagonal,
            "different_heavy": different,
        },
        "complete_renormalized_group_ledger": [
            "local sunset: proper local contractions are p independent; outer OS is evaluated using its twice differentiated integrable Schwinger representation",
            "six mixed sunset refinements: inherited cubic counterterm replaces the proper light bubble I by I-I0bar; subtract one local zero-external reference before the OS bound",
            "three same-heavy sunset refinements: inherited heavy-mass counterterm replaces the proper light bubble I by I-I0bar",
            "six different-heavy sunset refinements: all UV cores absent; integrate the full Schwinger domain before OS",
            "nested family: inner constant light and heavy tadpoles combine with fixed inner OS counterterms; Q=alpha+R; outer constant pieces vanish only after outer OS; remaining mixed inner pieces use the actual inherited kernel",
            "one-loop quartic reference and once-fixed finite potential contact insertions are local tadpoles and vanish under outer OS; heavy one-point counterterms also contribute only local light quadratic terms",
            "new two-loop light mass and kinetic counterterms impose the same mass-one, unit-residue OS conditions; total canonical interactions retain the inherited reference scheme",
        ],
        "two_loop_inverse_propagator": inverse,
        "scope": "The full fixed-order two-loop canonical light two-point function, with all 32 raw refinements and the required one-loop parameter and light quadratic insertions grouped before estimates. This is not yet an independent complete two-loop four-point/LSZ ledger, a two-loop source-aware derivative map or an all-orders propagator theorem.",
        "checks": checks,
    }
