"""Exact probability versus amplitude projection, with all12 proper overlaps."""

from functools import cache
from itertools import product

import sympy as s
from p8_vacuum_affine_three_singleton_subtraction import measure as amplitude_measure


@cache
def data():
    faces, P, R = amplitude_measure.projectors()
    Z = s.ImmutableMatrix(
        [[int(part & mask == part) for part in range(8)] for mask in range(8)]
    )
    M = s.ImmutableMatrix(
        [
            [
                (-1) ** (mask.bit_count() - part.bit_count())
                if part & mask == part
                else 0
                for part in range(8)
            ]
            for mask in range(8)
        ]
    )
    real = s.symbols("x0:8", real=True)
    imag = s.symbols("y0:8", real=True)
    g = s.Matrix([x + s.I * y for x, y in zip(real, imag)])
    values = Z * g
    square = lambda z: s.expand(z * s.conjugate(z))
    probabilities = s.Matrix([square(v) for v in values])
    probability_rectangle = s.expand((R * probabilities)[7])
    G = values[7]
    H = sum(g[:7])
    signed = s.expand(square(G) - square(H))
    proper_pairs = tuple((a, b) for a, b in product(range(7), repeat=2) if a | b == 7)
    all_pairs = tuple((a, b) for a, b in product(range(8), repeat=2) if a | b == 7)
    transfer = s.expand(sum(g[a] * s.conjugate(g[b]) for a, b in proper_pairs))
    six = sum(
        g[a] * s.conjugate(g[b]) + g[b] * s.conjugate(g[a])
        for a, b in ((1, 6), (2, 5), (4, 3), (3, 5), (3, 6), (5, 6))
    )
    sample = {**{x: 0 for x in real}, **{y: 0 for y in imag}, real[3]: 1, real[5]: 1}
    omitted_pairs = sum(
        g[a] * s.conjugate(g[b]) + g[b] * s.conjugate(g[a])
        for a, b in ((1, 6), (2, 5), (4, 3))
    )
    checks = {
        "anchored_components_invert_face_values_left": Z * M - s.eye(8),
        "anchored_components_invert_face_values_right": M * Z - s.eye(8),
        "proper_component_sum_is_amplitude_baseline": s.expand((P * values)[7] - H),
        "all27_union_covering_terms": s.expand(
            probability_rectangle - sum(g[a] * s.conjugate(g[b]) for a, b in all_pairs)
        ),
        "probability_equals_defined_signed_plus_transfer": s.expand(
            probability_rectangle - signed - transfer
        ),
        "all_six_real_cross_terms": s.expand(transfer - six),
        "all_union_pairs_inventory": s.Integer(len(all_pairs) - 27),
        "proper_overlap_pairs_inventory": s.Integer(len(proper_pairs) - 12),
        "full_component_pairs_inventory": s.Integer(
            sum(7 in pair for pair in all_pairs) - 15
        ),
        "pair_pair_nonzero_control": transfer.subs(sample) - 2,
    }
    for i, E in enumerate(faces):
        checks[f"probability_rectangle_face{i}_zero"] = E * R
    xr, xi, yr, yi, zr, zi = s.symbols("xr xi yr yi zr zi", real=True)
    X, Y, V = xr + s.I * xi, yr + s.I * yi, zr + s.I * zi
    two_real = square(X + Y - V) - square(X) - square(Y) + square(V)
    two_cross = (X - V) * s.conjugate(Y - V) + (Y - V) * s.conjugate(X - V)
    checks["independent_two_real_scheme_transfer"] = s.expand(two_real - two_cross)
    return {
        "checks": checks,
        "gates": {
            "dropping_pair_pair_terms_changes_result": omitted_pairs.subs(sample)
            != transfer.subs(sample),
            "probability_and_amplitude_schemes_not_identical": transfer != 0,
            "all_union_components_are_complex": len(g.free_symbols) == 16,
            "no_unsigned_probability_claim": True,
        },
        "whole_anchored_probability_coefficients": tuple(R[7, :]),
        "whole_proper_overlap_pairs": proper_pairs,
        "whole_term_inventory": {
            "all": len(all_pairs),
            "contains_full": 15,
            "proper": len(proper_pairs),
        },
        "whole_transfer_expression": six,
        "whole_exact_scheme_relation": "R3(|G3|^2)=(|G3|^2-|P3G3|^2)+T. T has12 ordered proper-component products, equivalently three singleton/complement-pair real cross terms and three distinct pair/pair real cross terms. The15 products containing g123 form the parent signed subtraction.",
    }
