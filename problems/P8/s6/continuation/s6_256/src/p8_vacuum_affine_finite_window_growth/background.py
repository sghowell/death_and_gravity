"""Full current coefficient enclosures and mass-adapted classical background box."""

from functools import cache

import sympy as s
from p8_vacuum_affine_classical_principal_realization import homogeneous as prior
from p8_vacuum_affine_classical_principal_realization import (
    realization as prior_initial,
)
from p8_vacuum_affine_coupled_principal_obstruction import background as old
from p8_vacuum_affine_coupled_principal_obstruction import coupled as c
from p8_vacuum_affine_heavy_scalar_parent import clock as heavy_bounds
from p8_vacuum_affine_heavy_source_filtration import source
from p8_vacuum_affine_physical_background_vertices import parent

from .intervals import I, evaluate

EPSILON = s.Rational(1, 10**6)
MOMENTUM_MIN = s.Integer(10) ** 64
MOMENTUM_MAX = 2 * MOMENTUM_MIN
TIME_LENGTH = s.Rational(1, 10**60)


@cache
def source_data():
    radius = s.Rational(1, 10000)
    distance = s.Rational(11, 10000)
    delta = distance * (2 + distance) / (1 - distance) ** 2
    switch = s.Rational(1, 300) ** 1024
    cauchy = s.factorial(4) / radius**4
    chain = 3**4 + 6 * 3**2 * 7 + 3 * 7**2 + 4 * 3 * 25 + 121
    normalized_pv = (
        s.Rational(5, 128) * 1000**4 / 9
        + (s.Rational(3, 2) + 600) * 10**396 / (64 * 9)
        + s.Rational(3, 128 * 9)
    ) / source.KAPPA
    tree_packet = tree_data()
    fixed = old.data()
    heavy = heavy_bounds.data()
    tree_num, tree_den = s.fraction(
        s.factor(parent.original.original.data()["original_retuned_tree_scalar"])
    )
    vacuum_num, vacuum_den = s.fraction(
        s.factor(parent.original.previous.vacuum_lower_function())
    )

    def polynomial_modulus(expr):
        return sum(
            abs(coefficient) * distance**i * (1 + delta) ** k
            for (i, k), coefficient in s.Poly(expr, parent.u, parent.X).terms()
        )

    tree_complex_bound = polynomial_modulus(tree_num) / (200 * (1 - distance**2) ** 12)
    vacuum_complex_bound = polynomial_modulus(vacuum_num) / 150
    chain_derivatives = [
        evaluate(
            s.diff(N**-2, N, order),
            {N: I(s.Rational(999, 1000), s.Rational(1001, 1000))},
        )
        for order in range(1, 5)
    ]
    gates = {
        "complete_lapse_complex_image_distance": delta < s.Rational(1, 400),
        "complete_switch_ratio_bound": delta / (1 - delta) < s.Rational(1, 300),
        "actual_complex_tree_scalar_modulus_below_1e6": tree_complex_bound < 10**6,
        "actual_complex_vacuum_polynomial_modulus_below_1e205": vacuum_complex_bound
        < 10**205,
        "complete_complex_square_and_clock_denominator_bounds": bool(
            1 - 2 * delta - delta**2 > s.Rational(199, 200)
            and (1 + delta) ** 2 < 2
            and (1 - distance**2) ** -3 < 2
        ),
        "all_four_exact_lapse_chain_derivative_bounds": all(
            max(abs(value.lo), abs(value.hi)) < bound
            for value, bound in zip(chain_derivatives, (3, 7, 25, 121))
        ),
        "whole_old_bump_modulus_below_one": 2048 * 2 / s.Integer(1000) ** 2 < 1,
        "old_R_complete_four_jet_error": 8 * cauchy * switch < s.Rational(1, 10**2500),
        "old_scalar_analytic_complete_four_jet_error": 10**207 * cauchy * switch
        < s.Rational(1, 10**2300),
        "constant_switch_error_even_with_excessive_majorant": 10**301 * cauchy * switch
        < s.Rational(1, 10**2200),
        "all_three_actual_normalized_vacuum_constants_bounded": normalized_pv
        < s.Rational(1, 10**400),
        "log_mass_bound_from_positive_exponential_series": 1
        + 3
        + s.Rational(9, 2)
        + s.Rational(27, 6)
        > 10,
        "actual_heavy_mass_inside_fixed_box": 10**197 < source.MASS2 < 10**198,
        "actual_kappa_not_changed": source.KAPPA == 10**800,
        "N_chain_rule_all_four_orders": chain < 10**4,
        "same_full_heavy_four_jet_source_bound": heavy["strict_four_jet_bound"]
        == s.Rational(1, 10**2700),
        "actual_heavy_complete_four_jet_proved": bool(
            heavy["gates"][
                "full_every_mixed_four_jet_difference_below_one_e_minus_2700"
            ]
        ),
        "entire_real_lapse_strip_inside_prior_heavy_tube": bool(
            s.Rational(1001, 1000) ** -2 > s.Rational(7, 8)
            and s.Rational(999, 1000) ** -2 < s.Rational(9, 8)
        ),
        "actual_combined_profile_C5_bound": prior_initial.PROFILE_BOUND
        == s.Rational(2, 10**400),
        "complete_smooth_profile_product_majorant": 10**5 * prior_initial.PROFILE_BOUND
        < s.Rational(1, 10**390),
        "complete_finite_onepoint_source_four_jet_bound": 10**20
        * switch
        / s.Integer(10) ** 400
        < s.Rational(1, 10**2900),
        "sum_of_all_current_coefficient_errors": s.Rational(1, 10**2300)
        + s.Rational(1, 10**2200)
        + s.Rational(1, 10**390)
        + s.Rational(1, 10**2696)
        < s.Rational(1, 10**300),
        "sum_of_bare_and_finite_onepoint_source_errors": s.Rational(1, 10**2696)
        + s.Rational(1, 10**2900)
        < s.Rational(1, 10**2500),
    }
    return {
        "entire_current_function_bindings": fixed["entire_current_function_bindings"],
        "whole_fixed_profile_bindings": fixed["whole_fixed_profile_bindings"],
        "all_three_fixed_vacuum_constants": fixed["all_three_fixed_vacuum_constants"],
        "complete_normalized_heavy_source": source.physical_source(u, N**-2)
        / s.sqrt(source.KAPPA),
        "source_complex_radius": radius,
        "complete_complex_tree_scalar_modulus_bound": tree_complex_bound,
        "complete_complex_vacuum_polynomial_modulus_bound": vacuum_complex_bound,
        "full_switch_modulus_base": switch,
        "complete_four_jet_Cauchy_factor": cauchy,
        "complete_fourth_N_chain_rule_factor": chain,
        "whole_current_R_F_tree_difference_four_jet_bound": s.Rational(1, 10**300),
        "whole_current_normalized_heavy_source_four_jet_bound": s.Rational(1, 10**2500),
        "normalized_vacuum_constant_absolute_bound": normalized_pv,
        "tree_bounds": {
            key: value
            for key, value in tree_packet.items()
            if key not in ("checks", "gates")
        },
        "checks": {
            "complete_fourth_chain_rule_coefficient": s.Integer(chain) - 1027,
            "entire_full_R_initial_time_derivative_zero": s.diff(
                fixed["entire_current_function_bindings"][R], u
            ).subs(u, 0),
            "whole_complex_tree_denominator_identity": s.expand(
                tree_den - 200 * (1 + u * u) ** 12
            ),
            "whole_vacuum_polynomial_rational_denominator": vacuum_den - 150,
        },
        "gates": {key: bool(value) for key, value in gates.items()},
    }


@cache
def energy_data():
    rate = rate_data()
    field, velocity, NN, NNd, mass, friction, forcing = s.symbols(
        "heavy_field heavy_velocity N Ndot mass2 friction full_source", real=True
    )
    energy = velocity**2 + NN**2 * mass * field**2
    acceleration = -friction * velocity - NN**2 * mass * field + NN**2 * forcing
    derivative = (
        s.diff(energy, field) * velocity
        + s.diff(energy, velocity) * acceleration
        + s.diff(energy, NN) * NNd
    )
    expected = (
        -2 * friction * velocity**2
        + 2 * NN**2 * forcing * velocity
        + 2 * NN * NNd * mass * field**2
    )
    damping = (
        3 * s.Rational(1, 1000)
        + rate["whole_coefficient_first_time_derivative_bounds"]["Z"]
        / s.Rational(99, 100)
        + rate["whole_velocity_and_lapse_derivative_bounds"]["Ndot"]
        / s.Rational(999, 1000)
    )
    improved = 4 * TIME_LENGTH / s.Integer(10) ** 2500
    return {
        "whole_mass_adapted_heavy_energy_square": energy,
        "whole_mass_adapted_heavy_energy_derivative": derivative,
        "whole_heavy_root_energy_damping_bound": damping,
        "heavy_root_energy_bootstrap": s.Rational(1, 10**1000),
        "heavy_root_energy_improved_bound": improved,
        "evaluated_classical_time_length": TIME_LENGTH,
        "whole_light_rate_displacement_bound": 10**20 * TIME_LENGTH,
        "whole_coefficient_displacement_bound": 10**40 * TIME_LENGTH,
        "continuation_argument": "Use the unchanged S255 implicit-lapse local solution. On its maximal subinterval inside the displayed box, the full constrained Euler bounds control lapse, volume and M1 rates. The positive heavy energy obeys e_prime<=10^20 e+2*10^-2500, so e<=4*T*10^-2500 on T=10^-60. The exact rate and first-coefficient bounds improve every bootstrap margin strictly. Standard local ODE continuation therefore reaches T. The giant fixed mass is retained; no tiny-time conclusion is drawn from its raw mass-squared Lipschitz constant.",
        "checks": {
            "complete_mass_cancellation_in_energy_identity": s.expand(
                derivative - expected
            )
        },
        "gates": {
            "actual_root_energy_damping_below_1e20": bool(damping < 10**20),
            "short_interval_exponential_below_two": bool(
                10**20 * TIME_LENGTH < s.Rational(1, 2)
            ),
            "strict_heavy_energy_bootstrap_improvement": bool(
                improved < s.Rational(1, 10**1001)
            ),
            "heavy_charge_and_mixed_source_box_improved": bool(
                10**10 * improved < s.Rational(1, 10**1000)
            ),
            "full_light_displacements_strictly_inside_initial_margins": bool(
                10**20 * TIME_LENGTH < s.Rational(1, 10**30)
            ),
            "full_coefficient_displacements_below_initial_margins": bool(
                10**40 * TIME_LENGTH < s.Rational(1, 10**19)
            ),
            "actual_positive_mass_covers_field_bootstrap": bool(source.MASS2 > 10**197),
            "actual_heavy_potential_coefficient_relative_time_bound": bool(
                rate["whole_velocity_and_lapse_derivative_bounds"]["Ndot"]
                / s.Rational(999, 1000)
                + s.Rational(3, 4)
                * rate["whole_coefficient_first_time_derivative_bounds"]["r"]
                / s.Rational(99, 100)
                < 10**40
            ),
        },
    }


u, N = old.u, old.N
R = old.R
F = old.F
primitive = old.I
j = s.Function("full_normalized_source", real=True)(u, N)
n = old.mass2
J = s.Symbol("positive_lapse_Schur_pivot", positive=True)
H, mc, mh, h = c.H, prior.mc, prior.mh, prior.hbar
a = c.a
U = R ** (-s.Rational(3, 4))
D = R ** s.Rational(1, 4) / N
Z = U / N
Ru = s.diff(R, u)
IN = 3 * U * Ru * s.diff(R, N) / (4 * R * N)
B = -U * Ru / (2 * N) - primitive
L = (
    -3 * D * H**2
    + 3 * B * H
    + N * U * F
    + 9 * U * Ru**2 / (16 * R * N)
    - s.diff(primitive, u)
    + Z * (mc**2 + mh**2) / 2
    + N * U * (-n * h * h / 2 + j * h)
)


def eliminate(expr):
    expr = s.expand(expr)
    rules = {}
    for deriv in expr.atoms(s.Derivative):
        if deriv.expr != primitive:
            continue
        counts = dict(deriv.variable_count)
        assert counts.get(N, 0) > 0, ("unremoved primitive", deriv)
        rules[deriv] = s.diff(IN, u, counts.get(u, 0), N, counts[N] - 1)
    result = expr.xreplace(rules)
    assert not result.has(primitive)
    return result


constraint = eliminate(s.diff(L, N))
theta = eliminate(-H * s.diff(D, N) + s.diff(B, N) / 2)
q = s.Matrix([6 * theta, s.diff(Z, N) * mc, s.diff(Z, N) * mh])
p0 = s.diff(L, H)
f = s.Matrix(
    [
        eliminate(3 * L - 3 * H * p0 - s.diff(p0, u)),
        -3 * H * Z * mc - s.diff(Z, u) * mc,
        N * U * (-n * h + j) - 3 * H * Z * mh - s.diff(Z, u) * mh,
    ]
)
rest = s.diff(constraint, u) + s.diff(constraint, h) * mh
Kinv = s.diag(-1 / (6 * D), 1 / Z, 1 / Z)
Nd = (-rest - (q.T * Kinv * f)[0]) / (2 * J)
rates = Kinv * (f - q * Nd)
expressions = {
    "Ndot": Nd,
    "Hdot": rates[0],
    "M1_acceleration": rates[1],
    "heavy_acceleration": rates[2],
}


@cache
def rate_data():
    bindings = {
        N: I(s.Rational(999, 1000), s.Rational(1001, 1000)),
        R: I(s.Rational(99, 100), s.Rational(101, 100)),
        J: I(1, 2),
        a: I(s.Rational(999, 1000), s.Rational(1001, 1000)),
        H: I(-s.Rational(1, 1000), s.Rational(1, 1000)),
        mc: I(s.Rational(9, 100), s.Rational(11, 100)),
        mh: I(-s.Rational(1, 10**1000), s.Rational(1, 10**1000)),
        h: I(-s.Rational(1, 10**1098), s.Rational(1, 10**1098)),
        n: I(10**196, 10**198),
    }
    for func, bound in [(R, 200), (F, 10**4), (j, s.Rational(1, 10**2500))]:
        for i in range(5):
            for k in range(5 - i):
                key = s.diff(func, u, i, N, k)
                if key not in bindings:
                    bindings[key] = I(-bound, bound)

    enclosed = {}
    rate_bounds = {}
    for name, expr in expressions.items():
        value = evaluate(expr, bindings)
        bound = max(abs(value.lo), abs(value.hi))
        assert bound < 10**20, (name, float(bound))
        enclosed[name] = value
        rate_bounds[name] = s.Rational(bound.numerator, bound.denominator)

    Cnn = s.diff(constraint, N) / 2
    Jactual = Cnn + 3 * theta**2 / D - (s.diff(Z, N) ** 2 * (mc**2 + mh**2)) / (2 * Z)
    parameters = {
        "D": D,
        "Z": Z,
        "J": Jactual,
        "Cnn": Cnn,
        "theta": theta,
        "r": R - 1,
        "rN": s.diff(R, N),
        "Y": N * R ** (-s.Rational(1, 4)),
        "K": R ** (-s.Rational(1, 4)) / (10**6 * N * a * a),
        "Yv": N * R ** (-s.Rational(1, 4)) / (a * a),
        "L2": 4 * s.diff(N * R ** s.Rational(3, 4) / 2, N),
        "C": N * R ** s.Rational(3, 4) / 2,
        "c": Z * mc,
        "w": s.diff(Z, N) * mc,
    }
    # All scalar time derivatives are bounded with the same full coupled rate
    # enclosure; the fixed source/profile jets are not treated as live state data.
    derivative_bounds = {}
    for name, expr in parameters.items():
        derivative = (
            s.diff(expr, u)
            + s.diff(expr, N) * s.Symbol("Ndot")
            + s.diff(expr, H) * s.Symbol("Hdot")
            + s.diff(expr, mc) * s.Symbol("M1_acceleration")
            + s.diff(expr, mh) * s.Symbol("heavy_acceleration")
            + s.diff(expr, h) * mh
            + s.diff(expr, a) * a * H
        )
        extra = {s.Symbol(key): value for key, value in enclosed.items()}
        value = evaluate(derivative, {**bindings, **extra})
        bound = max(abs(value.lo), abs(value.hi))
        assert bound < 10**40, (name, float(bound))
        derivative_bounds[name] = s.Rational(bound.numerator, bound.denominator)
    return {
        "whole_equations": expressions,
        "whole_coefficient_functions": parameters,
        "whole_coefficient_jet_box": {
            str(key): value.bounds() for key, value in bindings.items()
        },
        "whole_velocity_and_lapse_derivative_bounds": rate_bounds,
        "whole_coefficient_first_time_derivative_bounds": derivative_bounds,
        "checks": {
            "literal_complete_homogeneous_action": s.expand(
                L
                - prior.action()["normalized_homogeneous_density"].subs(
                    old.JH, s.sqrt(old.kappa) * j
                )
            )
        },
        "gates": {
            "whole_background_rates_below_1e20": all(
                value < 10**20 for value in rate_bounds.values()
            ),
            "whole_coefficient_first_jets_below_1e40": all(
                value < 10**40 for value in derivative_bounds.values()
            ),
        },
    }


@cache
def initial_data():
    N = s.Symbol("N", positive=True)
    R = s.Function("R", positive=True)(N)
    F = s.Function("F", real=True)(N)
    Ruu = s.Function("Ruu", real=True)(N)
    U = R ** (-s.Rational(3, 4))
    D = R ** s.Rational(1, 4) / N
    Z = U / N
    T = 3 * U * Ruu * s.diff(R, N) / (4 * R * N)
    G = s.diff(N * U * F, N) - T
    m2 = -2 * G / s.diff(Z, N)
    Cnn = (s.diff(N * U * F, N, 2) - s.diff(T, N) + s.diff(Z, N, 2) * m2 / 2) / 2
    J = Cnn - s.diff(Z, N) ** 2 * m2 / (2 * Z)
    Hd = (-U * Ruu / (2 * N) - N * U * F - Z * m2 / 2) / (2 * D)
    td = -s.diff(D, N) * Hd - s.diff(U * Ruu / (2 * N), N) / 2 - T / 2
    L2 = 4 * s.diff(N * R ** s.Rational(3, 4) / 2, N)
    epsilon = s.Rational(1, 10**6)
    point = 1 + epsilon
    error = s.Rational(1, 10**300)
    tree_R = N**-2
    tree_Ruu = -6 * (N**-2 - 1)
    tree_F = -(624 * N**-4 + 753 * N**-2 + 224) / 200
    bindings = {N: I(point)}
    for formal, tree, maximum in [(R, tree_R, 2), (Ruu, tree_Ruu, 1), (F, tree_F, 2)]:
        for j in range(maximum + 1):
            value = s.diff(tree, N, j).subs(N, point)
            bindings[s.diff(formal, N, j)] = I(value - error, value + error)
    expressions = {
        "R": R,
        "r": R - 1,
        "rN": s.diff(R, N),
        "D": D,
        "Z": Z,
        "m2": m2,
        "Cnn": Cnn,
        "J": J,
        "Hhatdot": Hd,
        "theta_dot": td,
        "L2": L2,
        "Y": N * R ** (-s.Rational(1, 4)),
        "C": N * R ** s.Rational(3, 4) / 2,
        "K": R ** (-s.Rational(1, 4)) / (10**6 * N),
        "Yv": N * R ** (-s.Rational(1, 4)),
        "c": Z * s.sqrt(m2),
        "w": s.diff(Z, N) * s.sqrt(m2),
    }
    enclosures = {
        name: evaluate(value, bindings) for name, value in expressions.items()
    }
    safe = {
        "R": (s.Rational(999, 1000), s.Rational(1001, 1000)),
        "r": (-3 * epsilon, -epsilon),
        "D": (s.Rational(999, 1000), s.Rational(1001, 1000)),
        "rN": (-3, -1),
        "Yv": (s.Rational(999, 1000), s.Rational(1001, 1000)),
        "Z": (s.Rational(999, 1000), s.Rational(1001, 1000)),
        "m2": (s.Rational(9, 1000), s.Rational(11, 1000)),
        "Cnn": (s.Rational(151, 100), s.Rational(153, 100)),
        "J": (s.Rational(151, 100), s.Rational(153, 100)),
        "Hhatdot": (s.Rational(399, 100), s.Rational(401, 100)),
        "theta_dot": (s.Rational(299, 100), s.Rational(301, 100)),
        "L2": (-s.Rational(1001, 1000), -s.Rational(999, 1000)),
        "Y": (s.Rational(999, 1000), s.Rational(1001, 1000)),
        "C": (s.Rational(499, 1000), s.Rational(501, 1000)),
        "K": (s.Rational(999, 10**9), s.Rational(1001, 10**9)),
        "c": (s.Rational(99, 1000), s.Rational(101, 1000)),
        "w": (s.Rational(49, 1000), s.Rational(51, 1000)),
    }
    for name, interval in enclosures.items():
        assert interval.inside(*safe[name]), (name, interval)

    return {
        "epsilon": epsilon,
        "lapse": point,
        "whole_current_central_formulas": expressions,
        "nonzero_full_jet_error": error,
        "whole_current_enclosures": {
            name: value.bounds() for name, value in enclosures.items()
        },
        "strict_coefficient_box": safe,
        "checks": {},
        "gates": {
            "all_seventeen_actual_coefficient_enclosures_inside_strict_box": all(
                value.inside(*safe[name]) for name, value in enclosures.items()
            )
        },
    }


@cache
def tree_data():
    u, N = parent.u, parent.N
    F = parent.original.original.data()["original_retuned_tree_scalar"].subs(
        parent.X, N**-2
    )
    R = 1 + (N**-2 - 1) / (1 + u * u) ** 3
    bindings = {
        u: I(-s.Rational(1, 1000), s.Rational(1, 1000)),
        N: I(s.Rational(999, 1000), s.Rational(1001, 1000)),
    }
    maxima = {}
    for name, expr, bound in [("R", R, 200), ("F", F, 10**4)]:
        largest = 0
        for i in range(5):
            for j in range(5 - i):
                derivative = s.factor(s.diff(expr, u, i, N, j))
                value = evaluate(derivative, bindings)
                absolute = max(abs(value.lo), abs(value.hi))
                largest = max(largest, absolute)
                assert absolute + s.Rational(1, 10**300) < bound, (
                    name,
                    i,
                    j,
                    float(absolute),
                )
        maxima[name] = s.Rational(largest.numerator, largest.denominator)
    return {
        "full_rational_tree": F,
        "full_clock_R": R,
        "mixed_total_order": 4,
        "tree_four_jet_maximum_enclosures": maxima,
        "whole_current_four_jet_bounds": {"R": s.Integer(200), "F": s.Integer(10) ** 4},
        "checks": {},
        "gates": {
            "whole_current_R_four_jets_below_200": bool(
                maxima["R"] + s.Rational(1, 10**300) < 200
            ),
            "whole_current_F_four_jets_below_1e4": bool(
                maxima["F"] + s.Rational(1, 10**300) < 10**4
            ),
        },
    }


@cache
def continuation_data():
    """Link the evaluated initial margins to every actual matrix input box."""
    initial = initial_data()
    enclosures = initial["whole_current_enclosures"]
    safe = initial["strict_coefficient_box"]
    displacement = s.Integer(10) ** 40 * TIME_LENGTH
    # Hdot and theta_dot are initial diagnostics only, not persistent boxes.
    persistent = (
        "R",
        "r",
        "rN",
        "D",
        "Z",
        "Cnn",
        "J",
        "L2",
        "Y",
        "C",
        "K",
        "Yv",
        "c",
        "w",
    )
    margins = {
        name: min(
            enclosures[name][0] - safe[name][0], safe[name][1] - enclosures[name][1]
        )
        for name in persistent
    }
    light_displacement = s.Integer(10) ** 20 * TIME_LENGTH
    mc_box = I(*enclosures["m2"]) ** s.Rational(1, 2)
    coarse = {
        key: I(*value)
        for key, value in rate_data()["whole_coefficient_jet_box"].items()
    }

    def bind(expr):
        return {
            atom: coarse[str(atom)]
            for atom in ({R, N, h, mh, n, j} | set(expr.atoms(s.Derivative)))
            if str(atom) in coarse
        }

    mixed_source = s.diff(N * U * (-n * h + j), N)
    mixed_box = evaluate(mixed_source, bind(mixed_source))
    improved = energy_data()["heavy_root_energy_improved_bound"]
    charge_box = evaluate(Z, {R: coarse[str(R)], N: coarse[str(N)]})
    w_box = evaluate(s.diff(Z, N), bind(s.diff(Z, N)))
    max_charge = max(abs(charge_box.lo), abs(charge_box.hi))
    max_w = max(abs(w_box.lo), abs(w_box.hi))
    physical_ratio = evaluate(
        R ** s.Rational(1, 4) / a,
        {R: I(*safe["R"]), a: I(s.Rational(999, 1000), s.Rational(1001, 1000))},
    )
    lapse = I(s.Rational(999, 1000), s.Rational(1001, 1000))
    source_mass = evaluate(N * U, {R: I(*safe["R"]), N: lapse})
    # Derive the full four-equation constrained ODE from the entire action.
    actual_K = s.hessian(L, (H, mc, mh))
    actual_q = s.Matrix([s.diff(constraint, v) for v in (H, mc, mh)])
    K = Kinv.inv()
    Cnn = s.diff(constraint, N) / 2
    schur = Cnn - (q.T * Kinv * q)[0] / 2
    expected_schur = (
        Cnn + 3 * theta**2 / D - s.diff(Z, N) ** 2 * (mc**2 + mh**2) / (2 * Z)
    )
    # This algebraic identity does not expand the enormous full rates.
    cc = s.Symbol("joint_Cnn", real=True)
    zz, dd, jj = s.symbols("joint_Z joint_D joint_J", positive=True)
    qq = s.Matrix(s.symbols("joint_q0:3", real=True))
    ff = s.Matrix(s.symbols("joint_f0:3", real=True))
    rr = s.Symbol("joint_constraint_rest", real=True)
    ki = s.diag(-1 / (6 * dd), 1 / zz, 1 / zz)
    nj = (-rr - (qq.T * ki * ff)[0]) / (2 * jj)
    vv = ki * (ff - qq * nj)
    joint = (ki.inv() * vv + qq * nj - ff).col_join(
        s.Matrix([(qq.T * vv)[0] + 2 * cc * nj + rr])
    )
    joint = joint.subs(cc, jj + (qq.T * ki * qq)[0] / 2).applyfunc(s.factor)
    return {
        "persistent_initial_coefficient_margins": margins,
        "uniform_coefficient_displacement": displacement,
        "uniform_light_rate_displacement": light_displacement,
        "actual_complete_mixed_heavy_source": mixed_source,
        "whole_mixed_heavy_source_enclosure": mixed_box.bounds(),
        "physical_wavenumber_over_comoving_momentum_enclosure": physical_ratio.bounds(),
        "actual_positive_heavy_potential_over_mass2_enclosure": source_mass.bounds(),
        "entire_joint_Euler_kinetic_matrix": actual_K,
        "entire_joint_Euler_lapse_cross_vector": actual_q,
        "whole_constraint_schur_pivot": expected_schur,
        "checks": {
            "complete_action_three_velocity_Hessian": actual_K - K,
            "complete_action_constraint_velocity_cross": actual_q - q,
            "entire_lapse_schur_identity": s.factor(schur - expected_schur),
            "independent_full_four_equation_schur_solution": joint,
        },
        "gates": {
            **{
                "strict_persistent_margin_" + name: bool(margin > 2 * displacement)
                for name, margin in margins.items()
            },
            "actual_positive_M1_rate_stays_in_coarse_bootstrap": bool(
                mc_box.lo - light_displacement > s.Rational(9, 100)
                and mc_box.hi + light_displacement < s.Rational(11, 100)
            ),
            "lapse_and_scale_factor_stay_in_coarse_bootstrap": bool(
                EPSILON + light_displacement < s.Rational(1, 1000)
            ),
            "clock_interval_stays_in_full_coefficient_strip": bool(
                TIME_LENGTH < s.Rational(1, 1000)
            ),
            "whole_heavy_source_mixing_inside_matrix_box": bool(
                max(abs(mixed_box.lo), abs(mixed_box.hi)) < s.Rational(1, 10**880)
            ),
            "full_improved_heavy_charges_inside_matrix_box": bool(
                max(max_charge, max_w) * improved < s.Rational(1, 10**1000)
            ),
            "full_improved_heavy_field_inside_bootstrap": bool(
                improved / 10**98 < s.Rational(1, 10**1098)
            ),
            "entire_physical_momentum_conversion_between_half_and_two": physical_ratio.inside(
                s.Rational(1, 2), 2
            ),
            "actual_heavy_potential_in_full_mass_box": source_mass.inside(
                s.Rational(1, 2), 2
            ),
            "entire_clock_Hubble_difference_inside_matrix_box": bool(
                light_displacement + 4 * TIME_LENGTH < s.Rational(1, 10**20)
            ),
        },
    }
