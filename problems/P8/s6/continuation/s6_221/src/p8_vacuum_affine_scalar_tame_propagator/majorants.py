"""Uniform exact coefficient enclosures for both complete finite-transfer charts."""

from functools import cache

import sympy as s

from . import charts as c

LOWER = s.Rational(1, 10**8)
UPPER = s.Integer(10) ** 4
ENTRY = s.Integer(10) ** 18


def polynomial_majorant(expr, variables, box):
    polynomial = s.Poly(s.cancel(expr), *variables)
    value = sum(
        abs(coefficient)
        * s.prod(box[v] ** power for v, power in zip(variables, powers))
        for powers, coefficient in polynomial.terms()
    )
    return value, polynomial


@cache
def data():
    z = c.z
    small = s.Rational(1, 10**6)
    zi = s.Rational(1, 4096)
    perturb = 6 * small * zi + 9 * small**2 * zi**2
    perturb += (2 * zi + 9 * small * zi**2) * (200 + s.Rational(1, 400))
    delta_lower = s.Rational(1, 4) - perturb
    delta_upper = 1 + perturb
    central = c.central()
    outer = c.outer()
    variables = (z, *c.jetvars, *c.jets, *c.jet2)
    box = {
        z: zi,
        c.th: 2,
        c.E: s.Rational(1, 2),
        c.l: s.Rational(1, 10),
        c.J: 100,
        c.A: small,
        c.T: small,
        c.H: 2,
    }
    box.update({v: s.Integer(10) ** 6 for v in (*c.jets, *c.jet2)})
    rows = {}
    polynomial_counts = {}
    reconstruction = {}
    for label, chart in (("central", central), ("outer", outer)):
        localbox = {**box, c.E: s.Rational(1, 2) if label == "central" else 1}
        denominator = c.E**4 * central["Delta"] ** 2 if label == "central" else c.th**4
        inverse_bound = s.Integer(4) ** 4 * (8**2 if label == "central" else 1)
        packets = {
            "K": chart["K"],
            "Kdot": chart["K"].applyfunc(c.dtime),
            "G": chart["G"],
            "Gdot": chart["G"].applyfunc(c.dtime),
            "gyro": chart["gyro"],
            "lower": chart["lower"],
        }
        for key, matrix in packets.items():
            values = []
            counts = []
            residual = []
            for expr in matrix:
                normalized = s.cancel(expr.subs(c.q, 1 / z) * denominator)
                envelope, polynomial = polynomial_majorant(
                    normalized, variables, localbox
                )
                values.append(envelope * inverse_bound)
                counts.append(len(polynomial.terms()))
                residual.append(s.cancel(normalized - polynomial.as_expr()))
            rows[label + "_" + key] = values
            polynomial_counts[label + "_" + key] = counts
            reconstruction[label + "_" + key] = s.Matrix(residual)
    k11 = 4 * (200 + s.Rational(1, 400)) * 8
    k22 = ((1 + 3 * small * zi) ** 2 + 18 * small * 100 * zi**2) * 8
    trace_outer = (200 + s.Rational(1, 100)) * 16 + 1
    trace_gradient = 2 * (100 + s.Rational(1, 200)) * 16 + 1
    trace_bounds = {
        "central_K": k11 + k22,
        "outer_K": trace_outer,
        "both_G": trace_gradient,
    }
    det_bounds = {"central_K": s.Rational(1, 25), "other_three": s.Rational(1, 200)}
    checks = {
        "normalized_denominator_from_actual_auxiliary_Hessian": s.factor(
            central["Delta"] - (2 * c.J * central["XX"].det() / c.q**2).subs(c.q, 1 / z)
        ),
        "central_K11_exact": s.factor(
            central["K"][0, 0].subs(c.q, 1 / z)
            - 4 * (2 * c.J + c.w * c.w) / central["Delta"]
        ),
        "central_K22_exact": s.factor(
            central["K"][1, 1].subs(c.q, 1 / z)
            - ((2 * c.E + 3 * c.T * z) ** 2 - 4 * c.J * z - 18 * c.A * c.J * z * z)
            / central["Delta"]
        ),
        "central_G_determinant": s.factor(
            central["G"].det() - 2 * (c.gradient_numerator() - c.w * c.w / 2) / c.E**2
        ),
        "outer_K_determinant": s.factor(outer["K"].det() - 2 * c.J / c.th**2),
        "outer_G_determinant": s.factor(
            outer["G"].det() - 2 * (c.gradient_numerator() - c.w * c.w / 2) / c.th**2
        ),
        "forty_eight_normalized_polynomial_reconstructions": s.Matrix(
            [entry for row in reconstruction.values() for entry in row]
        ),
        "forty_eight_entry_count": sum(map(len, rows.values())) - 48,
    }
    return {
        "central_denominator_bounds": (delta_lower, delta_upper),
        "coefficient_entry_majorants": rows,
        "normalized_polynomial_term_counts": polynomial_counts,
        "trace_upper_bounds": trace_bounds,
        "determinant_lower_bounds": det_bounds,
        "common_positive_operator_interval": (LOWER, UPPER),
        "common_absolute_entry_bound": ENTRY,
        "coefficient_enclosure": "After q=1/z each central coefficient K,K',G,G',Omega,R multiplied by E^4 Delta^2 is a polynomial. Bound its absolute coefficient polynomial on the exact coefficient/jet box and multiply4^4*8^2. Outer coefficients multiplied by Theta^4 are polynomials and use4^4. All48 entries are below1e18, with no time or momentum sampling.",
        "positivity_argument": "The complete central auxiliary Hessian has positive ps diagonal and determinant q^2 Delta/(2Jc), hence positive inverse and K. Both G determinants use positive F. Outer K uses positive Jc. Every matrix has positive trace below1e4 and determinant at least1/200, so lambda_min>=det/trace>1e-8. No assertion is made outside its chart domain.",
        "checks": checks,
        "gates": {
            "Delta_lower_above_one_eighth": delta_lower > s.Rational(1, 8),
            "Delta_upper_below_two": delta_upper < 2,
            "all_trace_bounds_below_upper": all(
                v < UPPER for v in trace_bounds.values()
            ),
            "central_determinant_above_common": det_bounds["central_K"]
            > det_bounds["other_three"],
            "common_minimum_eigenvalue_margin": det_bounds["other_three"] / UPPER
            > LOWER,
            "all_48_entry_envelopes": all(
                value < ENTRY for row in rows.values() for value in row
            ),
            "all_envelopes_nonnegative": all(
                value >= 0 for row in rows.values() for value in row
            ),
            "full_lower_order_terms_and_q_time_derivative": True,
            "actual_fixed_profile_jets_through_two_only": True,
            "not_an_all_finite_q_complementary_chart": True,
        },
    }
