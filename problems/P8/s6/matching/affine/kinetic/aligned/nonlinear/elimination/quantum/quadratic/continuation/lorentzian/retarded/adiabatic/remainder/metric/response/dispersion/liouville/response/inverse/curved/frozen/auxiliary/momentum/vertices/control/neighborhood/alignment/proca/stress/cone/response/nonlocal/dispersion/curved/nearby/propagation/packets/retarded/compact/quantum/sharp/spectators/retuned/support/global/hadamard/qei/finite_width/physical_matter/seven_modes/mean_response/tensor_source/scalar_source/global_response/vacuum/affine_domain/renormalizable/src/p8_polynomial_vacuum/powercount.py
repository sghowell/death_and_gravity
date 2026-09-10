"""All-order superficial power counting of the Gaussian-heavy model.

The exhaustive table is an exact integer consequence, not a loop computation.
The accompanying proof covers unbounded graph loop order and subgraphs.
"""

from fractions import Fraction

import sympy as sp


def exact_nonnegative(value):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction, sp.Rational)):
        raise TypeError("An exact nonnegative integer is required")
    q = sp.Rational(value)
    if q.q != 1 or q < 0:
        raise ValueError("An exact nonnegative integer is required")
    return int(q)


def degree(light_external, heavy_external, cubic_vertices):
    e, h, v = map(exact_nonnegative, (light_external, heavy_external, cubic_vertices))
    if e % 2 or v < h or (v - h) % 2:
        raise ValueError("Light parity and exact heavy-leg counting are required")
    return 4 - e - h - v


def candidates():
    # If D>=0, e+h+v<=4; v>=h. This finite enumeration is exhaustive
    # at every loop order, not a bound on the number of quartic vertices.
    return {
        (e, h, v): degree(e, h, v)
        for e in range(0, 5, 2)
        for h in range(5)
        for v in range(5)
        if v >= h and (v - h) % 2 == 0 and e + h + v <= 4
    }


def data():
    e, h, v, w, ip, ih = sp.symbols(
        "E_light E_heavy V_cubic V_quartic I_light I_heavy", integer=True
    )
    loop = ip + ih - v - w + 1
    raw = 4 * loop - 2 * (ip + ih)
    incidences = {ip: (2 * v + 4 * w - e) / 2, ih: (v - h) / 2}
    rows = candidates()
    heavy = {
        (a, b): max(D for (c, d, _), D in rows.items() if (c, d) == (a, b))
        for a, b, _ in rows
        if b > 0
    }
    return {
        "superficial_degree": 4 - e - h - v,
        "candidate_degree_table": [
            {
                "external_light": a,
                "external_heavy": b,
                "cubic_vertices": c,
                "max_momentum_degree": d,
            }
            for (a, b, c), d in sorted(rows.items())
        ],
        "heavy_dependent_counterterm_max_degree": {
            str(k): v for k, v in sorted(heavy.items())
        },
        "required_basis": [
            "vacuum constant",
            "H tadpole",
            "H squared mass",
            "H Phi squared cubic",
            "Phi squared mass",
            "light kinetic term",
            "Phi fourth power",
        ],
        "excluded_divergent_monomials": [
            "heavy kinetic term",
            "H cubed",
            "H fourth power",
            "H squared Phi squared",
        ],
        "qualifiers": "Only potentially divergent structures are listed; scaleless or kinematically absent diagrams can remove, not add, entries. Mass and tadpole insertions lower superficial degree. A finite constant heavy field normalization remains Gaussian. No all-energy continuum existence or loop-size conclusion follows.",
        "checks": {
            "topological_superficial_degree": sp.expand(
                raw.subs(incidences) - (4 - e - h - v)
            ),
            "heavy_leg_count": sp.expand((2 * ih + h - v).subs(incidences)),
            "light_leg_count": sp.expand((2 * ip + e - 2 * v - 4 * w).subs(incidences)),
            "heavy_onepoint_max_degree": sp.Integer(heavy[(0, 1)] - 2),
            "heavy_twopoint_max_degree": sp.Integer(heavy[(0, 2)]),
            "heavy_light_cubic_max_degree": sp.Integer(heavy[(2, 1)]),
        },
        "bounds": {
            "all_heavy_dependent_candidates_are_at_most_quadratic_in_H": set(heavy)
            == {(0, 1), (0, 2), (2, 1)},
            "heavy_kinetic_divergence_excluded_by_degree": heavy[(0, 2)] < 2,
            "only_light_wavefunction_divergence_allowed": max(
                D for (a, b, _), D in rows.items() if a == 2 and b == 0
            )
            == 2,
            "light_quartic_degree_zero": max(
                D for (a, b, _), D in rows.items() if a == 4 and b == 0
            )
            == 0,
        },
    }


def bad_cases():
    bad = (
        True,
        False,
        1.0,
        sp.Float(1),
        "1",
        sp.I,
        sp.oo,
        -sp.oo,
        sp.zoo,
        sp.nan,
        -1,
        sp.Rational(1, 2),
    )
    rows = [
        (
            "type_or_sign_" + str(j) + "_" + str(i),
            degree,
            tuple(v if k == j else (2, 0, 0)[k] for k in range(3)),
        )
        for j in range(3)
        for i, v in enumerate(bad)
    ]
    rows += [
        ("incidence_" + str(i), degree, v)
        for i, v in enumerate(((1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 2, 1)))
    ]
    return rows
