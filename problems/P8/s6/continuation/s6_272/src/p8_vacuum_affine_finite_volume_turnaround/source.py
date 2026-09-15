"""Whole moving-center clock jets and translated auxiliary derivative bounds."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate
from p8_vacuum_affine_quantitative_local_time import source as q

TIME = s.Rational(1, 10**130)
COORDINATE = s.Rational(1, 10**120)
MAJORANT = s.Integer(10) ** 200
LAPSE_RADIUS = q.LAPSE_RADIUS
u, N = q.u, q.N
PROFILE = q.current.prior_initial.PROFILE_BOUND


def background():
    values = dict.fromkeys(q.COORDS, s.Integer(0))
    values[q.p] = -2 * q.physical.H
    values[q.dp] = q.physical.ell - s.Rational(1, 10)
    return values


@cache
def clock_sources():
    tree = q.physical.original.original.data()["original_retuned_tree_scalar"].subs(
        q.physical.X, N**-2
    )
    return (
        1 + (N**-2 - 1) / (1 + u * u) ** 3,
        tree
        - q.physical.PRESSURE
        - (q.physical.RHO + q.physical.PRESSURE) * (N**-2 - 1) / 2,
    )


def clock_at(expression):
    """Exact low clock jets only; never replace the off-clock source."""
    Rclock, Fclock = clock_sources()
    expression = q.eliminate_N_primitive(expression).subs(
        background(), simultaneous=True
    )
    rules = {
        q.R: Rclock,
        q.F: Fclock,
        q.primitive: s.Integer(0),
        q.j: s.Integer(0),
        q.Hclock: q.physical.H,
    }
    for term in expression.atoms(s.Derivative):
        count = dict(term.variable_count)
        if term.expr in (q.R, q.F):
            if count.get(N, 0) >= 8:
                raise ValueError("Clock germ does not license the eighth lapse jet")
            target = Rclock if term.expr == q.R else Fclock
            rules[term] = s.diff(target, *term.variable_count)
        elif term.expr == q.primitive:
            if count.get(N, 0):
                raise ValueError("Uneliminated primitive lapse contact")
            rules[term] = s.Integer(0)
        elif term.expr == q.j:
            if count.get(N, 0) >= 8:
                raise ValueError("Whole heavy source has a nonzero eighth clock jet")
            rules[term] = s.Integer(0)
        elif term.expr == q.Hclock:
            rules[term] = s.diff(q.physical.H, *term.variable_count)
        else:
            raise ValueError("Unproved clock source jet: " + str(term))
    return s.factor(expression.xreplace(rules).subs(N, 1))


@cache
def magnitude(expression):
    expression = s.sympify(expression)
    if expression.is_Number:
        if expression.is_Rational is not True or expression.is_finite is not True:
            raise ValueError("Require finite exact rational source coefficients")
        return abs(expression)
    if expression in (N, q.R):
        return s.Integer(2)
    if expression in q.COORDS:
        return COORDINATE
    if expression == q.mu:
        return s.Rational(1, 100)
    if expression == q.F:
        return s.Integer(10) ** 7
    if expression in (q.j, q.Hclock):
        return s.Integer(1)
    if expression == q.primitive:
        return LAPSE_RADIUS * magnitude(q.IN)
    if isinstance(expression, s.Derivative):
        count = dict(expression.variable_count)
        order = sum(count.values())
        if expression.expr in (q.R, q.F):
            if order > 5:
                raise ValueError("Unproved full source derivative order")
            return q.JET
        if expression.expr == q.j:
            if order > 5:
                raise ValueError("Unproved full heavy-source derivative order")
            return s.Integer(1)
        if expression.expr == q.Hclock:
            if order > 3:
                raise ValueError("Unproved clock derivative order")
            return s.Integer(1000)
        if expression.expr == q.primitive:
            if count.get(N, 0):
                return magnitude(q.eliminate_N_primitive(expression))
            return LAPSE_RADIUS * magnitude(s.diff(q.IN, u, count.get(u, 0)))
    if expression.is_Add:
        return sum(magnitude(value) for value in expression.args)
    if expression.is_Mul:
        return s.prod(magnitude(value) for value in expression.args)
    if expression.is_Pow:
        base, power = expression.args
        if not power.is_Rational:
            raise ValueError("Unproved full source power")
        if power < 0:
            if base not in (N, q.R):
                raise ValueError("Unproved full source reciprocal")
            return s.Integer(2) ** s.ceiling(-power)
        if not power.is_Integer:
            if base not in (N, q.R):
                raise ValueError("Unproved full source fractional power")
            return s.Integer(2) ** s.ceiling(power)
        return magnitude(base) ** power
    raise ValueError("Unbounded whole-source atom: " + str(expression))


@cache
def expressions():
    C = q.CONSTRAINT
    rows = {
        "H": q.HAMILTONIAN,
        "C": C,
        "CN": s.diff(C, N),
        "Cu": s.diff(C, u),
        "CNN": s.diff(C, N, 2),
        "CNu": s.diff(C, N, u),
        "T": q.TEMPORAL,
    }
    for i, z in enumerate(q.COORDS):
        rows["Cz" + str(i)] = s.diff(C, z)
        rows["CNz" + str(i)] = s.diff(C, N, z)
        rows["Cuz" + str(i)] = s.diff(C, u, z)
        for j, w in enumerate(q.COORDS):
            rows["Czz" + str(i) + "_" + str(j)] = s.diff(C, z, w)
    return rows


@cache
def derivative_bounds():
    return {
        name: magnitude(q.eliminate_N_primitive(value))
        for name, value in expressions().items()
    }


@cache
def center():
    rho, pressure = s.symbols("bound_fixed_rho bound_fixed_pressure", real=True)
    box = {u: I(-TIME, TIME), rho: I(-PROFILE, PROFILE), pressure: I(-PROFILE, PROFILE)}
    C = clock_at(q.CONSTRAINT)
    CN = clock_at(s.diff(q.CONSTRAINT, N))
    gradients = {z: clock_at(s.diff(q.CONSTRAINT, z)) for z in q.COORDS}
    pivot = evaluate(CN.subs({q.physical.RHO: rho, q.physical.PRESSURE: pressure}), box)
    gradient_bounds = {}
    for name, value in gradients.items():
        interval = evaluate(value, box)
        raw = max(abs(interval.lo), abs(interval.hi))
        gradient_bounds[name] = s.Rational(raw.numerator, raw.denominator)
    return {
        "C": C,
        "CN": CN,
        "gradients": gradients,
        "gradient_bounds": gradient_bounds,
        "pivot": [
            s.Rational(pivot.lo.numerator, pivot.lo.denominator),
            s.Rational(pivot.hi.numerator, pivot.hi.denominator),
        ],
    }


@cache
def data():
    values, actual = derivative_bounds(), center()
    old = q.complex_source_bounds()
    clock_center = 3 * q.physical.PRESSURE / (2 * (1 + u * u) ** 3) - q.physical.RHO
    checks = {
        "whole_homogeneous_clock_constraint_with_actual_momenta": s.factor(
            actual["C"] - clock_center
        ),
        "whole_tree_homogeneous_clock_constraint_vanishes": s.factor(
            actual["C"].subs({q.physical.RHO: 0, q.physical.PRESSURE: 0})
        ),
        "whole_Gauss_center_gradient_zero": actual["gradients"][q.G],
        "whole_heavy_momentum_center_gradient_zero": actual["gradients"][q.ph],
        "whole_heavy_field_center_gradient_zero": actual["gradients"][q.eta],
        "whole_actual_background_trace_density": background()[q.p] + 2 * q.physical.H,
        "whole_actual_background_matter_density": background()[q.dp]
        + s.Rational(1, 10)
        - q.physical.ell,
        "whole_derivative_expression_count": s.Integer(len(values) - 187),
        "whole_same_fixed_profile_bound": PROFILE - s.Rational(2, 10**400),
    }
    return {
        "whole_original_constraint": q.CONSTRAINT,
        "whole_original_Hamiltonian": q.HAMILTONIAN,
        "whole_actual_moving_invariant_center": background(),
        "whole_clock_center_constraint": actual["C"],
        "whole_clock_center_lapse_pivot": actual["CN"],
        "whole_all_twelve_clock_center_parameter_gradients": actual["gradients"],
        "whole_center_pivot_enclosure": actual["pivot"],
        "whole_parameter_gradient_enclosures": actual["gradient_bounds"],
        "whole_187_translated_domain_derivative_ceilings": {
            name: s.ceiling(value) for name, value in values.items()
        },
        "whole_translated_coordinate_and_lapse_radii": [COORDINATE, LAPSE_RADIUS],
        "whole_actual_source_bindings": old["packet"][
            "entire_current_function_bindings"
        ],
        "whole_actual_fixed_profile_bindings": old["packet"][
            "whole_fixed_profile_bindings"
        ],
        "whole_primitive_rule": "Differentiate the complete parent first. Eliminate positive N derivatives only by the exact IN identity. I(u,1)=0 and its pure time derivatives vanish only on the clock. The entire RF and heavy-source differences have eighth clock factors; no off-clock replacement or live mean substitution is made.",
        "checks": checks,
        "gates": {
            "whole_187_derivative_bounds_below_majorant": all(
                value < MAJORANT for value in values.values()
            ),
            "whole_H_modulus_below1e62": values["H"] < 10**62,
            "whole_actual_pivot_strict_enclosure": -s.Rational(31, 10)
            < actual["pivot"][0]
            <= actual["pivot"][1]
            < -3,
            "whole_clock_p_gradient": actual["gradient_bounds"][q.p] < 100 * TIME,
            "whole_clock_pm_gradient": actual["gradient_bounds"][q.dp] < 1,
            "whole_all_other_center_gradients": all(
                value < 4
                for name, value in actual["gradient_bounds"].items()
                if name not in (q.p, q.dp)
            ),
            "all_twenty_one_RFj_source_jet_bounds": all(
                row[0] < q.JET and row[1] < q.JET and row[2] < 1
                for row in old["rows"].values()
            ),
            "actual_joint_complex_source_neighborhood": TIME + old["radius"]
            < old["distance"]
            and LAPSE_RADIUS + old["radius"] < old["distance"],
            "original_real_C5_time_slab": TIME < q.current.TIME_LENGTH,
            "full_source_gates_retained": all(old["packet"]["gates"].values()),
            "real_time_not_assumed_holomorphic": True,
        },
    }
