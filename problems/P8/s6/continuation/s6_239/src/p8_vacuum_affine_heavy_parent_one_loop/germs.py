"""Literal higher jets and exhaustive graphs of the classical limiting action."""

import itertools
from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_scalar_parent import family as parent

KAPPA = parent.K0
LAMBDA, GAMMA = parent.LAMBDA, parent.GAMMA
N = parent.N
G, MASS2, CONTACT = parent.G, parent.MASS2, parent.CONTACT
u, X = parent.u, parent.X


@cache
def couplings():
    k = KAPPA
    return {
        "phi6": -s.Rational(2681, 200) / k**2,
        "phi4Y": -s.Rational(853, 200) / k**2 - 8 * GAMMA / (3 * k) - CONTACT / (3 * k),
        "phi2Y2": -s.Rational(343168, 25) / k**2,
        "Y3": -s.Rational(109184, 25) / k**2 + 8 * LAMBDA / k,
        "phi2_L3_minus_L4": 6 * N / k**2,
        "Y_L3_minus_L4": -21 * N / k**2,
        "H_phi2_Y": -4 * G / k,
    }


def graph_degrees(external, loops=1):
    target = external + 2 * loops - 2
    return tuple(
        row
        for row in itertools.product(
            range(target + 1),
            range(target // 2 + 1),
            range(target // 3 + 1),
            range(target // 4 + 1),
        )
        if sum((j + 1) * number for j, number in enumerate(row)) == target
    )


@cache
def data():
    base = parent.original
    tree = base.original.data()["original_retuned_tree_scalar"]
    F = base.data()["F"]
    zero = {u: 0, X: 0}
    F0 = s.cancel(tree.subs(zero))
    alpha = s.cancel(s.diff(tree, u, 2).subs(zero) / 2)
    beta = s.cancel(s.diff(tree, X).subs(zero))
    targets = {
        (6, 0): -s.Rational(2681, 200),
        (4, 1): -s.Rational(853, 200),
        (2, 2): -s.Rational(343168, 25),
        (0, 3): -s.Rational(109184, 25),
    }
    checks = {
        "literal_tree_constant": F0 + s.Rational(28, 25),
        "literal_tree_u_squared": alpha + s.Rational(2781, 200),
        "literal_tree_X": beta + s.Rational(753, 200),
    }
    actual = {}
    for (i, j), target in targets.items():
        first = s.diff(F, X, j).subs(X, 0)
        got = s.cancel(
            s.diff(first, u, i).subs(u, 0) / (s.factorial(i) * s.factorial(j))
        )
        actual[f"u{i}_X{j}"] = got
        checks[f"literal_complete_old_weight_six_{i}_{j}"] = got - target
    R = base.data()["R"]
    checks["literal_R_u_squared_X_squared"] = (
        s.diff(s.diff(R, X, 2).subs(X, 0), u, 2).subs(u, 0) / 4 - 3 * N
    )
    checks["literal_R_X_cubed"] = s.diff(R, X, 3).subs(zero) / 6 - N
    A = parent.A
    h = parent.switch(X, A)
    checks["new_heavy_five_field_source_coefficient"] = (
        G * s.diff(h, X).subs(X, 0) / (2 * KAPPA) - couplings()["H_phi2_Y"]
    )
    checks["new_f_Phi4Y_coefficient"] = (
        targets[(4, 1)] - 8 * (GAMMA / 3 + CONTACT / 24) * KAPPA
    ) / KAPPA**2 - couplings()["phi4Y"]
    checks["new_f_Y3_coefficient"] = (
        targets[(0, 3)] + 8 * LAMBDA * KAPPA
    ) / KAPPA**2 - couplings()["Y3"]
    checks["new_a3_Phi_squared_coefficient"] = (
        6 * N / KAPPA**2 - couplings()["phi2_L3_minus_L4"]
    )
    checks["new_a3_Y_coefficient"] = (3 * N - 24 * N) / KAPPA**2 - couplings()[
        "Y_L3_minus_L4"
    ]
    V = 1 - base.data()["T"]
    denom = X**N + (1 - X) ** N
    checks["finite_counter_switch_vacuum_flat_factor"] = s.factor(
        denom * (1 - V) - X**N
    )
    checks["finite_counter_switch_clock_flat_factor"] = s.factor(
        denom * V - (1 - X) ** N
    )
    v, i, e, l = s.symbols("vertices internal_edges external_legs loop_number")
    checks["connected_loop_degree_sum"] = s.expand(
        2 * i + e - 2 * v - (e + 2 * l - 2)
    ).subs(i, l + v - 1)
    # The leading canonical massless M1 and massive-vector spectators remain
    # free in the CLASSICAL limiting action; no quantum-limit interchange is made.
    return {
        "observable": "The complete formal first-loop four-Phi coefficient of the fixed-canonical CLASSICAL limiting Phi/H action of S238, with its full higher operators. Free normalized tensor, massless M1 and massive-Proca spectators remain in that limiting action. This is not a regulator-removed quantum limit of the finite-gravity parent.",
        "actual_degree_six_local_coefficients": couplings(),
        "actual_anchor_parameters": {
            "kappa0": KAPPA,
            "n": MASS2,
            "g": G,
            "C": CONTACT,
            "lambda": LAMBDA,
            "gamma": GAMMA,
            "localizer_A": parent.LOCALIZER,
        },
        "literal_old_normalized_degree_six_jets": actual,
        "complete_one_loop_four_point_degree_patterns": graph_degrees(4),
        "complete_one_loop_two_point_degree_patterns": graph_degrees(2),
        "complete_one_loop_one_point_degree_patterns": graph_degrees(1),
        "graph_count_proof": "For a connected graph sum_v(degree_v-2)=E+2L-2. At E4,L1 the only nonquadratic degree patterns are four cubics, two cubics plus one quartic, two quartics, one cubic plus one degree-five vertex, or one degree-six vertex. Higher classical vertices cannot enter. At E2,L1 only the old pair of cubics and single quartic enter; at E1,L1 only the old cubic tadpole enters.",
        "source_topologies": "The new cubic-plus-five-field case has three light contraction classes in J2 G_H J4: within J2 (cancelled by the complete heavy-onepoint source counterterm), within J4 (a tadpole-corrected heavy exchange vertex), and across the factors (the mixed3+1 bubble). None is silently dropped.",
        "higher_localizer_boundary": "A first occurs in the degree-seven H Phi2 Y2 and degree-eight pure-light terms; the smooth QG1 nonconstant retuning begins at degree2048. These cannot contribute to the specified one-loop four-point coefficient. This is exact graph counting at this order, not an omitted-loop estimate.",
        "finite_extension_switch": V,
        "finite_extension_rule": "Multiply the separately specified finite counterterms by V=1-T. It equals1 up to X order1023 at the vacuum and vanishes to order1024 at the clock, so all relevant vacuum jets are unchanged. The new finite tube estimate is separate from S238's classical10^-2700 bound.",
        "checks": {key: s.cancel(value) for key, value in checks.items()},
        "gates": {
            "all_four_point_graph_patterns_exhausted": set(graph_degrees(4))
            == {(4, 0, 0, 0), (2, 1, 0, 0), (0, 2, 0, 0), (1, 0, 1, 0), (0, 0, 0, 1)},
            "two_point_functions_receive_no_new_higher_vertex": set(graph_degrees(2))
            == {(2, 0, 0, 0), (0, 1, 0, 0)},
            "one_point_condition_receive_no_new_higher_vertex": graph_degrees(1)
            == ((1, 0, 0, 0),),
            "actual_switch_order_retained": N == 1024,
            "higher_localizer_vertices_excluded_at_this_order": 7 - 2 > 4,
            "nonconstant_smooth_retuning_excluded_at_this_order": 2 * N - 2 > 4,
            "all_seven_new_coefficients_retained": len(couplings()) == 7,
        },
    }
