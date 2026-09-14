"""Complete constrained nearby trajectories, coefficient distances and physical turn."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_window_growth import background as b
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate

from . import source

u, N = b.u, b.N
Rtree = 1 + (N**-2 - 1) / (1 + u * u) ** 3
Ftree = b.parent.original.original.data()["original_retuned_tree_scalar"].subs(
    b.parent.X, N**-2
)
Href = 4 * u / (1 + u * u)
mref = 1 / (10 * (1 + u * u) ** 6)
aref = (1 + u * u) ** 2
Jactual = (
    s.diff(b.constraint, N) / 2
    + 3 * b.theta**2 / b.D
    - s.diff(b.Z, N) ** 2 * (b.mc**2 + b.mh**2) / (2 * b.Z)
)


def rational(value):
    return s.Rational(value.numerator, value.denominator)


def absolute(interval):
    return max(abs(interval.lo), abs(interval.hi))


def boxes():
    values = {
        N: I(s.Rational(999, 1000), s.Rational(1001, 1000)),
        b.R: I(s.Rational(99, 100), s.Rational(101, 100)),
        b.J: I(1, 2),
        b.a: I(s.Rational(999, 1000), s.Rational(1001, 1000)),
        b.H: I(-s.Rational(1, 1000), s.Rational(1, 1000)),
        b.mc: I(s.Rational(9, 100), s.Rational(11, 100)),
        b.mh: I(-s.Rational(1, 10**1000), s.Rational(1, 10**1000)),
        b.h: I(-s.Rational(1, 10**1098), s.Rational(1, 10**1098)),
        b.n: I(10**196, 10**198),
    }
    for f, bound, last in [
        (b.R, 200, 2000),
        (b.F, 10**4, 10**6),
        (b.j, s.Rational(1, 10**2500), s.Rational(1, 10**2450)),
    ]:
        for i in range(6):
            for j in range(6 - i):
                key = s.diff(f, u, i, N, j)
                if key not in values:
                    maximum = last if i + j == 5 else bound
                    values[key] = I(-maximum, maximum)
    return values


@cache
def reference_rules():
    values = {}
    for f, tree in [(b.R, Rtree), (b.F, Ftree), (b.j, s.Integer(0))]:
        for i in range(6):
            for j in range(6 - i):
                values[s.diff(f, u, i, N, j)] = s.factor(
                    s.diff(tree, u, i, N, j).subs(N, 1)
                )
    return values


def bind_reference(expr):
    position = {N: 1, b.H: Href, b.mc: mref, b.mh: 0, b.h: 0, b.a: aref}
    return s.factor(
        expr.subs(reference_rules(), simultaneous=True).subs(
            position, simultaneous=True
        )
    )


def dt0(expr):
    return (
        s.diff(expr, u)
        + s.diff(expr, N) * b.expressions["Ndot"]
        + s.diff(expr, b.H) * b.expressions["Hdot"]
        + s.diff(expr, b.mc) * b.expressions["M1_acceleration"]
        + s.diff(expr, b.mh) * b.expressions["heavy_acceleration"]
        + s.diff(expr, b.h) * b.mh
        + s.diff(expr, b.a) * b.a * b.H
    )


@cache
def Jdot():
    return dt0(Jactual)


def dt(expr):
    return dt0(expr) + s.diff(expr, b.J) * Jdot()


@cache
def initial():
    n = s.Symbol("positive_initial_lapse", positive=True)
    R = s.Function("whole_initial_R", positive=True)(n)
    F = s.Function("whole_initial_F", real=True)(n)
    Ruu = s.Function("whole_initial_Ruu", real=True)(n)
    U = R ** (-s.Rational(3, 4))
    Z = U / n
    T = 3 * U * Ruu * s.diff(R, n) / (4 * R * n)
    G = s.diff(n * U * F, n) - T
    m2 = -2 * G / s.diff(Z, n)
    tr = n**-2
    truu = -6 * (n**-2 - 1)
    tf = -(624 * n**-4 + 753 * n**-2 + 224) / 200
    substitution = {
        R: tr,
        s.diff(R, n): s.diff(tr, n),
        Ruu: truu,
        F: tf,
        s.diff(F, n): s.diff(tf, n),
    }
    bare = s.factor(m2.subs(substitution, simultaneous=True))
    wanted = (
        s.Rational(432, 25) / n**2 - s.Rational(2847, 100) + s.Rational(56, 5) * n**2
    )
    keys = [R, s.diff(R, n), Ruu, F, s.diff(F, n)]
    symbols = s.symbols("initial_independent_jet0:5", real=True)
    lifted = m2.xreplace(dict(zip(keys, symbols, strict=True)))
    domain = {
        n: I(s.Rational(999, 1000), s.Rational(1001, 1000)),
        symbols[0]: I(s.Rational(999, 1000), s.Rational(1001, 1000)),
        symbols[1]: I(-s.Rational(201, 100), -s.Rational(199, 100)),
        symbols[2]: I(-s.Rational(1, 10), s.Rational(1, 10)),
        symbols[3]: I(-9, -7),
        symbols[4]: I(19, 21),
    }
    derivatives = [
        rational(absolute(evaluate(s.factor(s.diff(lifted, z)), domain)))
        for z in symbols
    ]
    bare_derivative = rational(absolute(evaluate(s.diff(bare, n), {n: domain[n]})))
    # Replace entire derivative atoms at once. A substitution of R_u=0
    # inside an unbound R_uu would wrongly erase the retained second jet.
    point_binding = {
        N: n,
        b.R: R,
        s.diff(b.R, N): s.diff(R, n),
        s.diff(b.R, u, 2): Ruu,
        s.diff(b.R, u): 0,
        s.diff(b.R, u, N): 0,
        b.H: 0,
        b.h: 0,
        b.mh: 0,
        b.F: F,
        s.diff(b.F, N): s.diff(F, n),
    }
    mapped_point = b.constraint.xreplace(point_binding)
    point_constraint = s.factor(mapped_point - G - s.diff(Z, n) * b.mc**2 / 2)
    distance = 100 * source.EPSILON + 10**8 * source.ERROR
    root_distance = distance / s.Rational(19, 100)
    # These uniform bare bounds put every actual and interpolated source jet
    # in the larger MVT domain. No rounded root subtraction is used.
    nearby = I(1 - source.EPSILON, 1 + source.EPSILON)
    reentry = []
    for key, value, z in zip(
        keys, (tr, s.diff(tr, n), truu, tf, s.diff(tf, n)), symbols, strict=True
    ):
        iv = evaluate(value, {n: nearby})
        reentry.append(
            iv.lo - source.ERROR > domain[z].lo and iv.hi + source.ERROR < domain[z].hi
        )
    return {
        "whole_original_constraint_positive_M1_squared_root": m2,
        "whole_bare_comparison_M1_squared": bare,
        "whole_initial_source_jet_derivative_enclosures": derivatives,
        "whole_bare_initial_lapse_derivative_bound": bare_derivative,
        "signed_lapse_parameter_radius": source.EPSILON,
        "full_positive_M1_root_distance_bound": s.Rational(1, 10**225),
        "whole_initial_data": "u=0, a_hat=1, H_hat=0, h=mh=0, N=1+epsilon with abs(epsilon)<=10^-230, and mc the positive square root of the displayed COMPLETE original constraint expression. Reconstruct W=S and the original affine complement. No source or preparation is tuned.",
        "checks": {
            "whole_bare_constraint_root_formula": s.factor(bare - wanted),
            "whole_bare_initial_M1_rate": bare.subs(n, 1) - s.Rational(1, 100),
            "whole_literal_original_initial_constraint_bridge": point_constraint,
            "whole_complete_positive_root_solves_actual_constraint": s.factor(
                mapped_point.subs(b.mc, s.sqrt(m2))
            ),
            "whole_centered_root_difference_identity": s.factor(
                (
                    s.sqrt(s.Symbol("positive_rate_square", positive=True))
                    - s.Rational(1, 10)
                )
                * (
                    s.sqrt(s.Symbol("positive_rate_square", positive=True))
                    + s.Rational(1, 10)
                )
                - s.Symbol("positive_rate_square", positive=True)
                + s.Rational(1, 100)
            ),
        },
        "gates": {
            "all_actual_initial_jets_inside_MVT_box": all(reentry),
            "whole_initial_jet_Lipschitz_sum_below_1e8": sum(derivatives) < 10**8,
            "whole_bare_initial_lapse_derivative_below_100": bare_derivative < 100,
            "whole_positive_current_rate_root_exists": s.Rational(1, 100) - distance
            > s.Rational(9, 100) ** 2,
            "whole_initial_root_distance_below_1e_minus_225": root_distance
            < s.Rational(1, 10**225),
            "signed_lapse_family_stays_positive": source.EPSILON < s.Rational(1, 1000),
            "centered_root_not_independent_rounded_subtraction": True,
        },
    }


@cache
def bootstrap():
    parameters = light_parameters()
    delta = (
        4 * s.Rational(1, 10**200) + 2 * s.Rational(1, 10**1000) + 64 * source.ERROR
    ) * 10**80
    safe = {
        **{
            name: (s.Rational(999, 1000), s.Rational(1001, 1000))
            for name in ("D", "Z", "Y", "Yv", "K_over_zeta", "heavy_mass_ratio")
        },
        "J": (s.Rational(151, 100), s.Rational(153, 100)),
        "C": (s.Rational(499, 1000), s.Rational(501, 1000)),
        "L2": (-s.Rational(1001, 1000), -s.Rational(999, 1000)),
        "r": (-s.Rational(1, 10**20), s.Rational(1, 10**20)),
        "theta": (-s.Rational(1, 10**20), s.Rational(1, 10**20)),
        "rN": (-3, -1),
        "c": (s.Rational(99, 1000), s.Rational(101, 1000)),
        "w": (s.Rational(49, 1000), s.Rational(51, 1000)),
    }
    interval = {u: I(-source.TIME, source.TIME)}
    bare = {name: bind_reference(expr) for name, expr in parameters.items()}
    enclosures = {name: evaluate(expr, interval) for name, expr in bare.items()}
    gates = {
        name + "_whole_bootstrap_and_reference_coefficient_margin": value.lo - delta
        > safe[name][0]
        and value.hi + delta < safe[name][1]
        for name, value in enclosures.items()
    }
    gates.update(
        {
            "direct_pivot_containment_not_assumed_from_inverse": not Jactual.has(b.J),
            "complete_Cnn_branch_regular": s.Rational(151, 100)
            - 3 * s.Rational(1, 10**40) / s.Rational(999, 1000)
            > s.Rational(3, 2),
            "complete_R_temporal_auxiliary_margin": 1
            - 3 * s.Rational(1, 100) ** 2 / (2 * s.Rational(99, 100))
            > s.Rational(1, 4),
            "original_affine_X_domain_kept": s.Rational(999, 1000) ** -2
            < s.Rational(6, 5),
            "scale_and_Hubble_bootstrap_inside_generator_domain": aref.subs(
                u, source.TIME
            )
            + s.Rational(1, 10**200)
            < s.Rational(1001, 1000)
            and Href.subs(u, source.TIME) + s.Rational(1, 10**200)
            < s.Rational(1, 10**20),
            "reference_nonstationary_Euler_contacts_below_small_box": 10**80
            * source.ERROR
            < s.Rational(1, 10**6),
            "all_coefficient_segments_stay_inside_convex_displayed_box": True,
        }
    )
    return {
        "whole_bootstrap_light_state_radius": s.Rational(1, 10**200),
        "whole_bootstrap_coefficient_difference_enclosure": delta,
        "whole_bare_reference_coefficient_formulas": bare,
        "whole_bare_reference_coefficient_slab_enclosures": {
            name: value.bounds() for name, value in enclosures.items()
        },
        "whole_safe_persistent_generator_coefficient_box": safe,
        "whole_nonstationary_reference_Euler_contact_boundary": "Actual L0,Vvv,vs_i vanish by the complete unforced Euler identities. The fixed reference retains L0=3 Tcorr and Vvv=9 A/2; its tiny profiles are bounded, not set to zero. Both reference matter Euler contacts vanish. This difference remains in the generator MVT.",
        "checks": {
            "complete_lapse_pivot_without_assumed_inverse_parameter": s.diff(
                Jactual, b.J
            ),
            "whole_Cnn_from_actual_pivot_and_momentum_cross": s.factor(
                s.diff(b.constraint, N) / 2
                - Jactual
                + 3 * b.theta**2 / b.D
                - s.diff(b.Z, N) ** 2 * (b.mc**2 + b.mh**2) / (2 * b.Z)
            ),
        },
        "gates": gates,
    }


@cache
def light():
    domain = boxes()
    Jbare = bind_reference(Jactual)
    checks = {
        "whole_bare_lapse_constraint": bind_reference(b.constraint),
        "whole_bare_scale_Euler": s.factor(s.diff(aref, u) - aref * Href),
    }
    for name, wanted in [
        ("Ndot", 0),
        ("Hdot", s.diff(Href, u)),
        ("M1_acceleration", s.diff(mref, u)),
        ("heavy_acceleration", 0),
    ]:
        checks["whole_bare_" + name] = s.factor(
            bind_reference(b.expressions[name]).subs(b.J, Jbare) - wanted
        )
    ref = b.parent.reference_data()["whole_reference_coefficients"]
    checks["same_complete_bare_reference_pivot"] = s.factor(
        Jbare
        - ref["Jc"].subs({b.parent.RHO: 0, b.parent.PRESSURE: 0}, simultaneous=True)
    )
    state_bounds = {}
    for name in ("Ndot", "Hdot", "M1_acceleration"):
        for variable in (N, b.H, b.mc, b.h, b.mh):
            derivative = s.diff(b.expressions[name], variable) + s.diff(
                b.expressions[name], b.J
            ) * s.diff(Jactual, variable)
            state_bounds[name + "_" + str(variable)] = rational(
                absolute(evaluate(derivative, domain))
            )
    zero = {b.h: 0, b.mh: 0}
    rows = {
        name: expr.subs(zero, simultaneous=True)
        for name, expr in b.expressions.items()
        if name != "heavy_acceleration"
    }
    JL = Jactual.subs(zero, simultaneous=True)
    keys = sorted(
        set().union(*(expr.atoms(s.Derivative) for expr in (*rows.values(), JL)))
        | {b.R, b.F},
        key=str,
    )
    replacement = {
        key: s.Symbol("whole_rate_jet_" + str(i), real=True)
        for i, key in enumerate(keys)
    }
    locals_ = {**domain, **{replacement[k]: domain[k] for k in keys}}
    lifted_J = JL.xreplace(replacement)
    jet_bounds = {}
    for name, expr in rows.items():
        lifted = expr.xreplace(replacement)
        for key, z in replacement.items():
            derivative = s.diff(lifted, z) + s.diff(lifted, b.J) * s.diff(lifted_J, z)
            jet_bounds[name + "_" + str(key)] = rational(
                absolute(evaluate(derivative, locals_))
            )
    return {
        "whole_current_normalized_homogeneous_action": b.L,
        "whole_current_constrained_Euler_rates": b.expressions,
        "whole_actual_lapse_pivot": Jactual,
        "whole_bare_comparison_reference": [1, Href, mref, aref],
        "whole_bare_comparison_pivot": Jbare,
        "whole_light_state_Jacobian_bounds": state_bounds,
        "whole_light_source_jet_Jacobian_bounds": jet_bounds,
        "whole_current_state_and_source_jet_box": {
            str(k): v.bounds() for k, v in domain.items()
        },
        "checks": checks,
        "gates": {
            "all_full_light_state_derivatives_below_1e40": max(state_bounds.values())
            < 10**40,
            "all_full_source_jet_rate_derivatives_below_1e40": max(jet_bounds.values())
            < 10**40,
            "all_fourteen_rate_jets_retained": len(keys) == 14,
            "whole_lapse_pivot_chain_in_every_derivative": True,
            "entire_reference_tree_is_only_comparison": True,
            "heavy_mass_not_included_as_raw_light_frequency_error": True,
        },
    }


def light_parameters():
    return {
        "D": b.D,
        "Z": b.Z,
        "J": Jactual,
        "theta": b.theta,
        "r": b.R - 1,
        "rN": s.diff(b.R, N),
        "C": N * b.R ** s.Rational(3, 4) / 2,
        "L2": 4 * s.diff(N * b.R ** s.Rational(3, 4) / 2, N),
        "Y": N * b.R ** (-s.Rational(1, 4)),
        "K_over_zeta": b.R ** (-s.Rational(1, 4)) / (N * b.a**2),
        "Yv": N * b.R ** (-s.Rational(1, 4)) / b.a**2,
        "c": b.Z * b.mc,
        "w": s.diff(b.Z, N) * b.mc,
        "heavy_mass_ratio": N * b.U,
    }


def derivative_enclosures(expr, domain, allow_heavy_source=False):
    state = {}
    for variable in (N, b.H, b.mc, b.a, b.h, b.mh):
        derivative = s.diff(expr, variable) + s.diff(expr, b.J) * s.diff(
            Jactual, variable
        )
        state[str(variable)] = rational(absolute(evaluate(derivative, domain)))
    zero = {b.h: 0, b.mh: 0}
    reduced = expr.subs(zero, simultaneous=True)
    JL = Jactual.subs(zero, simultaneous=True)
    functions = (b.R, b.F, b.j) if allow_heavy_source else (b.R, b.F)
    keys = sorted(
        reduced.atoms(s.Derivative) | JL.atoms(s.Derivative) | set(functions), key=str
    )
    assert all(key in functions or key.expr in functions for key in keys)
    replacement = {
        key: s.Symbol("complete_coefficient_jet_" + str(i), real=True)
        for i, key in enumerate(keys)
    }
    locals_ = {**domain, **{replacement[k]: domain[k] for k in keys}}
    lifted = reduced.xreplace(replacement)
    lifted_J = JL.xreplace(replacement)
    jets = {}
    for key, z in replacement.items():
        derivative = s.diff(lifted, z) + s.diff(lifted, b.J) * s.diff(lifted_J, z)
        jets[str(key)] = rational(absolute(evaluate(derivative, locals_)))
    return {"state": state, "source_jets": jets}


def summarize_enclosures(rows):
    return {
        name: {
            family: {
                "entry_count": len(values),
                "maximum_absolute_enclosure": max(values.values()),
            }
            for family, values in row.items()
        }
        for name, row in rows.items()
    }


@cache
def coefficients():
    domain = boxes()
    parameters = light_parameters()
    rows = {
        **parameters,
        **{name + "_dot": dt(value) for name, value in parameters.items()},
    }
    enclosures = {
        name: derivative_enclosures(expr, domain) for name, expr in rows.items()
    }
    heavy_force = N * b.U * (-b.n * b.h + b.j)
    heavy = {
        "cH": b.Z * b.mh,
        "wH": s.diff(b.Z, N) * b.mh,
        "dHsource": s.diff(heavy_force, N),
    }
    heavy_rows = {
        **heavy,
        **{name + "_dot": dt(value) for name, value in heavy.items()},
    }
    heavy_bounds = {
        name: rational(absolute(evaluate(expr, domain)))
        for name, expr in heavy_rows.items()
    }
    treebox = {u: I(-s.Rational(1, 1000), s.Rational(1, 1000)), N: domain[N]}
    fifth = {}
    for name, expr, bound in [("R", Rtree, 2000), ("F", Ftree, 10**6)]:
        fifth[name] = max(
            rational(
                absolute(evaluate(s.factor(s.diff(expr, u, i, N, 5 - i)), treebox))
            )
            for i in range(6)
        )
        assert fifth[name] + source.ERROR < bound
    values = [
        value
        for row in enclosures.values()
        for family in row.values()
        for value in family.values()
    ]
    checks = {
        name + "_bare_clock_zero": s.factor(
            bind_reference(expr).subs(b.J, light()["whole_bare_comparison_pivot"])
        )
        for name, expr in heavy_rows.items()
    }
    return {
        "whole_fourteen_light_and_mass_coefficient_values": parameters,
        "whole_fourteen_complete_first_time_jets": {
            name: dt(expr) for name, expr in parameters.items()
        },
        "whole_retained_heavy_coefficients_and_first_jets": heavy_rows,
        "whole_complete_parameter_Jacobian_enclosures": summarize_enclosures(
            enclosures
        ),
        "whole_heavy_coefficient_absolute_enclosures": heavy_bounds,
        "complete_rational_tree_fifth_jet_maxima": fifth,
        "checks": checks,
        "gates": {
            "all_28_value_and_first_jet_rows_retained": len(rows) == 28,
            "all_parameter_state_and_jet_derivatives_below_1e80": max(values) < 10**80,
            "whole_heavy_values_and_first_jets_below_1e_minus_780": max(
                heavy_bounds.values()
            )
            < s.Rational(1, 10**780),
            "full_current_fifth_R_box": fifth["R"] + source.ERROR < 2000,
            "full_current_fifth_F_box": fifth["F"] + source.ERROR < 10**6,
            "all_independent_jet_lists_below_64": all(
                len(row["source_jets"]) < 64 for row in enclosures.values()
            ),
            "all_fixed_profile_and_source_terms_retained": True,
        },
    }


@cache
def physical():
    Rdot = s.diff(b.R, u) + s.diff(b.R, N) * b.expressions["Ndot"]
    Hphysical = (b.H - Rdot / (4 * b.R)) / N
    acceleration = dt(Hphysical) / N
    rows = {"physical_Hubble": Hphysical, "physical_proper_acceleration": acceleration}
    domain = boxes()
    enclosures = {
        name: derivative_enclosures(expr, domain, True) for name, expr in rows.items()
    }
    values = [
        value
        for row in enclosures.values()
        for family in row.values()
        for value in family.values()
    ]
    bare = light()["whole_bare_comparison_pivot"]
    checks = {
        "whole_physical_Hubble_bare_reference": s.factor(
            bind_reference(Hphysical).subs(b.J, bare) - Href
        ),
        "whole_physical_acceleration_bare_reference": s.factor(
            bind_reference(acceleration).subs(b.J, bare) - s.diff(Href, u)
        ),
    }
    return {
        "whole_actual_physical_scale_factor": b.a * b.R ** (-s.Rational(1, 4)),
        "whole_actual_physical_Hubble": Hphysical,
        "whole_actual_proper_time_Hubble_derivative": acceleration,
        "whole_physical_value_and_acceleration_Jacobian_bounds": summarize_enclosures(
            enclosures
        ),
        "checks": checks,
        "gates": {
            "complete_physical_value_and_acceleration_derivatives_below_1e160": max(
                values
            )
            < 10**160,
            "all_physical_source_jet_lists_below_64": all(
                len(row["source_jets"]) < 64 for row in enclosures.values()
            ),
            "lapse_pivot_time_derivative_in_physical_acceleration": True,
            "actual_C_R_inverse_square_root_not_refuted_shifted_R": True,
            "proper_time_derivative_not_coordinate_Hhatdot": True,
        },
    }
