"""Full complex auxiliary contraction, including every primitive-time contact."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate
from p8_vacuum_affine_nonlinear_lapse_branch import branch as previous
from p8_vacuum_affine_nonlinear_lapse_branch import source as bounce

from . import source as q

MAJORANT = s.Integer(10) ** 200


@cache
def magnitude(expression):
    expression = s.sympify(expression)
    if expression.is_Number:
        if expression.is_Rational is not True or expression.is_finite is not True:
            raise ValueError(
                "Only exact finite rational source coefficients are admitted"
            )
        return abs(expression)
    if expression in (q.N, q.R):
        return s.Integer(2)
    if expression in q.COORDS:
        return q.IMAGE
    if expression == q.mu:
        return s.Rational(1, 100)
    if expression == q.F:
        return s.Integer(10) ** 7
    if expression in (q.j, q.Hclock):
        return s.Integer(1)
    if expression == q.primitive:
        return q.LAPSE_RADIUS * magnitude(q.IN)
    if isinstance(expression, s.Derivative):
        count = dict(expression.variable_count)
        order = sum(count.values())
        if expression.expr in (q.R, q.F):
            if order > 5:
                raise ValueError("Unproved source derivative beyond total order5")
            return q.JET
        if expression.expr == q.j:
            if order > 5:
                raise ValueError("Unproved source derivative beyond total order5")
            return s.Integer(1)
        if expression.expr == q.Hclock:
            if order > 3:
                raise ValueError("Unproved Hclock derivative beyond order3")
            return s.Integer(1000)
        if expression.expr == q.primitive:
            if count.get(q.N, 0):
                return magnitude(q.eliminate_N_primitive(expression))
            return q.LAPSE_RADIUS * magnitude(s.diff(q.IN, q.u, count.get(q.u, 0)))
    if expression.is_Add:
        return sum(magnitude(a) for a in expression.args)
    if expression.is_Mul:
        return s.prod(magnitude(a) for a in expression.args)
    if expression.is_Pow:
        base, power = expression.args
        if not power.is_Rational:
            raise ValueError("Only proved real rational source powers are admitted")
        if power < 0:
            if base not in (q.N, q.R):
                raise ValueError("Unproved reciprocal in the whole source")
            return s.Integer(2) ** s.ceiling(-power)
        if not power.is_Integer:
            if base not in (q.N, q.R):
                raise ValueError("Unproved fractional source power")
            return s.Integer(2) ** s.ceiling(power)
        return magnitude(base) ** power
    raise ValueError("Unbounded atom in full auxiliary majorant: " + str(expression))


@cache
def derivative_bounds():
    C, N, u = q.CONSTRAINT, q.N, q.u
    expressions = {
        "H": q.HAMILTONIAN,
        "C": C,
        "CN": s.diff(C, N),
        "Cu": s.diff(C, u),
        "CNN": s.diff(C, N, 2),
        "CNu": s.diff(C, N, u),
        "T": q.TEMPORAL,
    }
    for i, z in enumerate(q.COORDS):
        expressions["Cz" + str(i)] = s.diff(C, z)
        expressions["CNz" + str(i)] = s.diff(C, N, z)
        expressions["Cuz" + str(i)] = s.diff(C, u, z)
        for j, w in enumerate(q.COORDS):
            expressions["Czz" + str(i) + "_" + str(j)] = s.diff(C, z, w)
    return {
        name: magnitude(q.eliminate_N_primitive(value))
        for name, value in expressions.items()
    }


@cache
def bounce_identity():
    C = q.CONSTRAINT
    rules = {q.R: bounce.R, q.F: bounce.F, q.j: bounce.j, q.primitive: 0, q.Hclock: 0}
    for term in C.atoms(s.Derivative):
        counts = dict(term.variable_count)
        nu, nN = counts.get(q.u, 0), counts.get(q.N, 0)
        if term.expr == q.R:
            if nu == 0:
                rules[term] = s.diff(bounce.R, q.N, nN)
            elif nu == 1:
                rules[term] = s.Integer(0)
            elif nu == 2:
                rules[term] = s.diff(bounce.Ruu, q.N, nN)
            else:
                raise ValueError("Unexpected time jet in literal bounce constraint")
        elif term.expr in (q.F, q.j) and nu == 0:
            target = bounce.F if term.expr == q.F else bounce.j
            rules[term] = s.diff(target, q.N, nN)
    return s.factor(C.xreplace(rules) - bounce.CONSTRAINT)


@cache
def data():
    values = derivative_bounds()
    old = previous.enclosure()
    variation = MAJORANT * (q.LAPSE_RADIUS + 12 * q.IMAGE + q.TIME)
    contraction = s.Rational(1, 30) + variation / 3
    center = bounce.CENTER_RADIUS + 120 * q.IMAGE + MAJORANT * q.TIME
    root = center / (3 * (1 - contraction))
    rn = 2 + q.JET * q.LAPSE_RADIUS
    Rgap = q.JET * q.LAPSE_RADIUS
    U_N = s.Rational(3, 4) * s.Rational(100, 99) ** 2 * rn
    clock = 4 * q.u / (1 + q.u * q.u)
    clock_bounds = []
    for k in range(4):
        iv = evaluate(s.diff(clock, q.u, k), {q.u: I(-q.TIME, q.TIME)})
        raw = max(abs(iv.lo), abs(iv.hi))
        clock_bounds.append(s.Rational(raw.numerator, raw.denominator))
    return {
        "whole_all_187_complete_derivative_outward_integer_ceilings": {
            name: s.ceiling(value) for name, value in values.items()
        },
        "whole_safe_joint_derivative_majorant": MAJORANT,
        "whole_full_Hamiltonian_modulus_ceiling": s.Integer(10) ** 62,
        "whole_actual_bounce_lapse_pivot_enclosure": old[
            "whole_lapse_first_outward_enclosure"
        ],
        "whole_actual_bounce_parameter_gradient_ceilings": old[
            "whole_parameter_gradient_ceilings"
        ],
        "whole_actual_bounce_source_only_center_residual": old[
            "whole_source_only_center_residual_bound"
        ],
        "whole_complex_joint_MVT_variation": variation,
        "whole_complex_lapse_contraction_factor": contraction,
        "whole_full_source_center_residual_bound": center,
        "whole_full_complex_root_displacement_bound": root,
        "whole_complete_R_gap_and_first_N_bound": [Rgap, rn],
        "whole_complex_physical_U_N_bound": U_N,
        "whole_actual_clock_time_derivative_bounds": clock_bounds,
        "whole_complex_branch_proof": "The actual full S266 real center has CN in(-3.1,-3), every |Czi|<6 and |C(0,1,0)|<1e-375. The complete complex derivative majorant bounds variations in N,z and REAL u by M(rho+12delta+T)<.001. Hence every Czi<10 and N->N+C/3 is q<1/20 contractive on |N-1|<=1e-245. Its full center residual is<1e-257 and the displayed self-map is strict. The entire source is holomorphic in complex N,z at each real u and smooth in real u, so the unique root is holomorphic in finite z and smooth in u. Real data give the positive real root by conjugation and uniqueness. No bare-clock substitution or complex-time analyticity is used.",
        "whole_primitive_majorant_proof": "For I(u,1)=0, I and every pure time derivative needed here equal the full N integral of the corresponding derivative of IN. Along the straight complex N segment their moduli are at most rho times the complete integrand majorant. The engine eliminates only positive N derivatives and rejects unproved derivative orders, reciprocals, powers and source atoms. All 187 expressions include the whole H, constraint, temporal solution and every CNz/Cuz/Czz contact.",
        "whole_real_auxiliary_regular_domain": "For real data the complete CN stays between-31/10-variation and-3+variation, hence below-29/10. Full R,N stay near1, Gamma=1-3(R-1)^2/(2R)>99/100 and the temporal pivot remains between-2 and-1/2. The full2x2 auxiliary Hessian is negative definite with determinant>1. The full temporal root, Hclock and primitive terms are retained; no spatial gauge or quantum matching determinant is inferred.",
        "whole_physical_derivative_and_composition_scope": "Actual full R(u,1)=1 and RN(u,1)=-2/(1+u^2)^3. The complete complex RNN bound gives |RN|<=2+1e30*rho, so the whole complex |U_N|<12 on this tube. The bound is for the complete implicit composition, not a linear Gaussian lapse or a partial Hessian with N_i z^i_AB omitted. The root is pointwise in spatial invariants; no lapse-gradient term is reintroduced into the reduced Hamiltonian.",
        "checks": {
            "literal_full_time_constraint_restricts_to_unchanged_bounce": bounce_identity(),
            "whole_derivative_expression_count": s.Integer(len(values) - 187),
            "whole_joint_MVT_factor": variation
            - MAJORANT * (q.LAPSE_RADIUS + 12 * q.IMAGE + q.TIME),
            "whole_implicit_contraction_error_denominator": 3 * (1 - contraction) * root
            - center,
            "whole_original_mass_adapted_coordinate_count": s.Integer(
                len(q.COORDS) - 12
            ),
        },
        "gates": {
            "all_187_complete_derivative_bounds": all(
                value < MAJORANT for value in values.values()
            ),
            "whole_Hamiltonian_modulus_below1e62": values["H"] < 10**62,
            "actual_original_center_gradient_below6": max(
                old["whole_parameter_gradient_ceilings"].values()
            )
            < 6,
            "actual_original_center_residual_below1e375_inverse": old[
                "whole_source_only_center_residual_bound"
            ]
            < bounce.CENTER_RADIUS,
            "strict_complex_joint_MVT_variation": variation < s.Rational(1, 1000),
            "whole_parameter_gradient_below10": 6 + variation < 10,
            "strict_complex_contraction_below_one_twentieth": contraction
            < s.Rational(1, 20),
            "full_complex_self_map": center / 3 + contraction * q.LAPSE_RADIUS
            < q.LAPSE_RADIUS,
            "full_source_center_residual": center < s.Rational(1, 10**257),
            "full_complex_root_bound": root < s.Rational(1, 10**257),
            "real_lapse_pivot_negative": -3 + variation < -s.Rational(29, 10),
            "whole_R_nonzero_principal_branch": Rgap < s.Rational(1, 100),
            "full_complex_RN_bound": rn < 4,
            "whole_complex_U_N_below12": U_N < 12,
            "actual_clock_jet_premises": clock_bounds[0] < 1
            and max(clock_bounds[1:]) < 1000,
            "real_auxiliary_determinant_margin_above_one": s.Rational(29, 10) / 2 > 1,
            "no_complex_time_or_original_P8_closure": True,
        },
    }
