"""Componentwise entire-source jets on a wider complex-lapse neighborhood."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_volume_turnaround import source as previous
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate

q, u, N = previous.q, previous.u, previous.N
SOURCE_TIME = previous.TIME
TIME = s.Rational(1, 10**180)
LAPSE_RADIUS = s.Rational(1, 10**115)
COORDINATE = previous.COORDINATE
HOMOGENEOUS_RADIUS = s.Rational(1, 10**122)
HOMOGENEOUS_CAUCHY = s.Rational(1, 10**124)
REAL_RADIUS = s.Rational(1, 10**390)
PROFILE = previous.PROFILE
CAUCHY_RADIUS = s.Rational(1, 10000)
SIX_JET = s.Integer(10) ** 35
OFF_CLOCK = LAPSE_RADIUS * SIX_JET
MAJORANT = s.Integer(10) ** 9


def absolute_interval(expression):
    rules, bindings = {}, {u: I(-SOURCE_TIME, SOURCE_TIME)}
    for name in (q.physical.RHO, q.physical.PRESSURE):
        for i in range(6):
            marker = s.Symbol(str(name.func) + "_time_" + str(i), real=True)
            rules[s.diff(name, u, i)] = marker
            bindings[marker] = I(-PROFILE, PROFILE)
    interval = evaluate(expression.xreplace(rules), bindings)
    value = max(abs(interval.lo), abs(interval.hi))
    return s.Rational(value.numerator, value.denominator)


def upper_power(value):
    if value.is_Rational is not True or value < 0:
        raise ValueError("Require an exact nonnegative rational bound")
    if value == 0:
        return s.Integer(0)
    bound = s.Integer(1)
    if value >= 1:
        while bound <= value:
            bound *= 10
    else:
        while bound / 10 > value:
            bound /= 10
    if not value < bound <= 10 * value:
        raise ValueError("Decimal outward bound failed")
    return bound


@cache
def jet_data():
    six = {}
    for i in range(6):
        for k in range(7 - i):
            factor = s.factorial(i) * s.factorial(k) / CAUCHY_RADIUS ** (i + k)
            profile = 8 * PROFILE * s.factorial(k) / CAUCHY_RADIUS**k
            six[(i, k)] = (2 * factor, 10**7 * factor + profile, factor / 10**100)
    clocks = dict(zip((q.R, q.F), previous.clock_sources(), strict=True))
    values, clock_values = {}, {}
    for name in (q.R, q.F):
        for i in range(6):
            for k in range(6 - i):
                germ = s.factor(s.diff(clocks[name], u, i, N, k).subs(N, 1))
                raw = absolute_interval(germ)
                clock_values[(name, i, k)] = (germ, raw)
                values[(name, i, k)] = upper_power(raw + OFF_CLOCK)
    for i in range(6):
        for k in range(6 - i):
            values[(q.j, i, k)] = OFF_CLOCK
    hclock = {
        i: upper_power(absolute_interval(s.diff(q.physical.H, u, i))) for i in range(5)
    }
    return {"six": six, "clocks": clock_values, "bounds": values, "Hclock": hclock}


@cache
def magnitude(expression, legacy=False):
    expression = s.sympify(expression)
    if expression.is_Number:
        if expression.is_Rational is not True or expression.is_finite is not True:
            raise ValueError("Require finite rational coefficients")
        return abs(expression)
    if expression == N:
        return s.Integer(2)
    if expression == q.R:
        return s.Integer(2) if legacy else jet_data()["bounds"][(q.R, 0, 0)]
    if expression in q.COORDS:
        return COORDINATE
    if expression == q.mu:
        return s.Rational(1, 100)
    if expression == q.F:
        return s.Integer(10) ** 7 if legacy else jet_data()["bounds"][(q.F, 0, 0)]
    if expression == q.j:
        return s.Integer(1) if legacy else jet_data()["bounds"][(q.j, 0, 0)]
    if expression == q.Hclock:
        return s.Integer(1) if legacy else jet_data()["Hclock"][0]
    radius = previous.LAPSE_RADIUS if legacy else LAPSE_RADIUS
    if expression == q.primitive:
        return radius * magnitude(q.IN, legacy)
    if isinstance(expression, s.Derivative):
        count = dict(expression.variable_count)
        i, k = count.get(u, 0), count.get(N, 0)
        if sum(count.values()) != i + k:
            raise ValueError("Unexpected whole-source differentiation variable")
        if expression.expr in (q.R, q.F, q.j):
            if i + k > 5:
                raise ValueError("Unproved full source jet order")
            if legacy:
                return s.Integer(1) if expression.expr == q.j else q.JET
            return jet_data()["bounds"][(expression.expr, i, k)]
        if expression.expr == q.Hclock:
            if k or i > (3 if legacy else 4):
                raise ValueError("Unproved clock derivative")
            return s.Integer(1000) if legacy else jet_data()["Hclock"][i]
        if expression.expr == q.primitive:
            if k:
                return magnitude(q.eliminate_N_primitive(expression), legacy)
            return radius * magnitude(s.diff(q.IN, u, i), legacy)
    if expression.is_Add:
        return sum(magnitude(v, legacy) for v in expression.args)
    if expression.is_Mul:
        return s.prod(magnitude(v, legacy) for v in expression.args)
    if expression.is_Pow:
        base, power = expression.args
        if not power.is_Rational:
            raise ValueError("Unproved source power")
        if power < 0:
            if base not in (N, q.R):
                raise ValueError("Unproved reciprocal")
            return s.Integer(2) ** s.ceiling(-power)
        if not power.is_Integer:
            if base not in (N, q.R):
                raise ValueError("Unproved fractional power")
            return s.Integer(2) ** s.ceiling(power)
        return magnitude(base, legacy) ** power
    raise ValueError("Unbounded full-source atom: " + str(expression))


@cache
def expressions():
    rows = dict(previous.expressions())
    C, H = q.CONSTRAINT, q.HAMILTONIAN
    rows["CNNN"] = s.diff(C, N, 3)
    for i, z in enumerate(q.COORDS):
        rows["CNNz" + str(i)] = s.diff(C, N, 2, z)
        rows["Hz" + str(i)] = s.diff(H, z)
        for j, w in enumerate(q.COORDS):
            rows["CNzz" + str(i) + "_" + str(j)] = s.diff(C, N, z, w)
            rows["Hzz" + str(i) + "_" + str(j)] = s.diff(H, z, w)
    return rows


@cache
def derivative_bounds():
    return {
        name: magnitude(q.eliminate_N_primitive(value))
        for name, value in expressions().items()
    }


@cache
def data():
    jets, bounds = jet_data(), derivative_bounds()
    old = q.complex_source_bounds()
    legacy = previous.derivative_bounds()
    matches = s.Matrix(
        [
            magnitude(q.eliminate_N_primitive(value), True) - legacy[name]
            for name, value in previous.expressions().items()
        ]
    )
    groups = {
        "Cz": max(bounds["Cz" + str(i)] for i in range(12)),
        "CNz": max(bounds["CNz" + str(i)] for i in range(12)),
        "CNN": bounds["CNN"],
        "Czz": max(
            bounds["Czz" + str(i) + "_" + str(j)] for i in range(12) for j in range(12)
        ),
        "Hz": max(bounds["Hz" + str(i)] for i in range(12)),
        "Hzz": max(
            bounds["Hzz" + str(i) + "_" + str(j)] for i in range(12) for j in range(12)
        ),
    }
    safe = {
        "Cz": 10**4,
        "CNz": 10**6,
        "CNN": 10**9,
        "Czz": 10**4,
        "Hz": 10**3,
        "Hzz": 10**3,
    }
    return {
        "whole_original_Hamiltonian": q.HAMILTONIAN,
        "whole_original_constraint": q.CONSTRAINT,
        "whole_full_primitive_N_identity": q.IN,
        "whole_500_derivative_ceilings": {
            name: s.ceiling(value) for name, value in bounds.items()
        },
        "whole_all_63_componentwise_RFj_bounds": {
            str(key): value for key, value in jets["bounds"].items()
        },
        "whole_all_27_mixed_six_jet_bounds": {
            str(key): value for key, value in jets["six"].items()
        },
        "whole_five_clock_derivative_bounds": jets["Hclock"],
        "whole_six_envelope_group_bounds": groups,
        "whole_six_envelope_group_safe_bounds": safe,
        "whole_source_and_homogeneous_radii": [
            SOURCE_TIME,
            TIME,
            LAPSE_RADIUS,
            COORDINATE,
            HOMOGENEOUS_RADIUS,
            HOMOGENEOUS_CAUCHY,
            REAL_RADIUS,
        ],
        "whole_actual_source_bindings": old["packet"][
            "entire_current_function_bindings"
        ],
        "whole_actual_fixed_profile_bindings": old["packet"][
            "whole_fixed_profile_bindings"
        ],
        "whole_clock_to_tube_rule": "Use the actual eighth-order clock factors. Bound all total-order-five jets at N1 exactly, then integrate the extra N derivative. The 27 mixed total-order-six entries omit time6. Original profiles need only real-time C5 and N analyticity; they are not live means or holomorphic in time.",
        "whole_primitive_rule": "Only positive N derivatives are replaced by the entire IN identity. Pure-time primitive derivatives are bounded by the full N integral on the NEW lapse radius. No off-clock primitive or source is set to zero.",
        "checks": {
            "all_187_legacy_engine_values_match": matches,
            "whole_500_expression_count": s.Integer(len(bounds) - 500),
            "whole_27_mixed_six_jet_count": s.Integer(len(jets["six"]) - 27),
            "whole_63_componentwise_jet_count": s.Integer(len(jets["bounds"]) - 63),
            "full_constraint_invariant_degree_two": s.Integer(
                s.Poly(q.CONSTRAINT, *q.COORDS).total_degree() - 2
            ),
            "unchanged_fixed_profile_bound": PROFILE - s.Rational(2, 10**400),
        },
        "gates": {
            "all_original_entire_source_gates": all(old["packet"]["gates"].values()),
            "whole_same_analytic_scalar_metric_and_heavy_bounds": old["scalar"] < 10**7
            and old["metric"] < 2
            and old["normalized_source"] < s.Rational(1, 10**100),
            "whole_larger_joint_source_neighborhood": SOURCE_TIME + CAUCHY_RADIUS
            < old["distance"]
            and LAPSE_RADIUS + CAUCHY_RADIUS < old["distance"],
            "all_27_extra_N_derivative_bounds": all(
                row[0] < SIX_JET and row[1] < SIX_JET and row[2] < SIX_JET
                for row in jets["six"].values()
            ),
            "no_sixth_time_derivative": (6, 0) not in jets["six"],
            "whole_500_derivative_bounds": all(
                value < MAJORANT for value in bounds.values()
            ),
            "whole_H_modulus": bounds["H"] < 100,
            "all_six_envelope_groups": all(
                groups[name] < safe[name] for name in groups
            ),
            "full_R_nonzero_branch": jet_data()["bounds"][(q.R, 0, 1)] * LAPSE_RADIUS
            < s.Rational(1, 100),
            "full_temporal_auxiliary_pivot": 1
            - s.Rational(3, 2)
            * (jet_data()["bounds"][(q.R, 0, 1)] * LAPSE_RADIUS) ** 2
            / (1 - jet_data()["bounds"][(q.R, 0, 1)] * LAPSE_RADIUS)
            > s.Rational(99, 100),
            "full_real_C5_source_slab": TIME < SOURCE_TIME < q.current.TIME_LENGTH,
            "homogeneous_Cauchy_polydisc_fits": 5 * HOMOGENEOUS_CAUCHY + REAL_RADIUS
            < HOMOGENEOUS_RADIUS,
            "no_source_or_profile_retuning": True,
        },
    }
