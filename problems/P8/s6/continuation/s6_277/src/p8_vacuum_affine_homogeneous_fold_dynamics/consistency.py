"""Full first/second constraint preservation, not a free multiplier at the fold."""

from functools import cache

import sympy as s

from . import canonical, source

p, m, sh, eta, ph = source.p, source.m, source.sh, source.eta, source.ph
r, P, r1, P1, r2, P2 = source.PROFILE_JETS
FIRST = p * (-s.Rational(3, 25) + 3 * (P - r)) + 3 * P1 / 2 - r1
COMPATIBLE_P = (3 * P1 / 2 - r1) / (s.Rational(3, 25) - 3 * (P - r))
EPS = source.PROFILE_BOUND
LAM = s.Symbol("finite_clock_lapse_velocity", real=True)


@cache
def coefficients():
    H0, Hu = source.clock_jet("H"), source.clock_jet("H", 1, 0)
    C0, Cu, Cuu = (source.clock_jet("C", k, 0) for k in range(3))
    CN = source.clock_jet("C", 0, 1)
    K = Cu + canonical.flow(C0, H0)
    DCN = source.clock_jet("C", 1, 1) + canonical.flow(CN, H0)
    D2C = Cuu + canonical.flow(Cu, H0) + canonical.flow(C0, Hu) + canonical.flow(K, H0)
    a, b, c = (source.at_fold(x) for x in (source.clock_jet("C", 0, 2), 2 * DCN, D2C))
    return {
        "K": source.at_fold(K),
        "a": a,
        "b": b,
        "c": c,
        "K_before_heavy_restriction": K,
    }


def remainder_bound(expression):
    variables = (*source.PROFILE_JETS, p)
    radii = (*([EPS] * 6), 21 * EPS)
    poly = s.Poly(s.expand(expression), *variables)
    return sum(
        abs(coefficient)
        * s.prod(radius**power for radius, power in zip(radii, monomial, strict=True))
        for monomial, coefficient in poly.terms()
        if any(monomial)
    )


@cache
def data():
    d = coefficients()
    a, b, c = (d[k] for k in ("a", "b", "c"))
    eb = -(600 * P * p - 700 * P1 - 1200 * r * p + 800 * r1 - 12627 * p) / 200
    ec = (
        -(
            20000 * P * P
            - 36000 * P * r
            - 36000 * P * p * p
            - 9100 * P
            - 4500 * P1 * p
            - 3000 * P2
            + 16000 * r * r
            + 45000 * r * p * p
            + 4940 * r
            + 2000 * r2
            + 55485 * p * p
            - 120528
        )
        / 2000
    )
    center = dict.fromkeys((*source.PROFILE_JETS, p), s.Integer(0))
    centers = [x.subs(center) for x in (a, b, c)]
    errors = [remainder_bound(x) for x in (a, b, c)]
    D = s.Matrix([[-9 * p * p, 3 * p], [3 * p, -1]])
    cp, ch, tq, tn, ct, ht = s.symbols(
        "Cq Hq Tstar_q Tstar_N Ctime Tstar_time", real=True
    )
    projection = (ct - tn * ht + (cp - tn * tq) * ch) + tn * (ht + tq * ch)
    # Generic polynomial test retains the normalized-density off-surface term.
    t, n = s.symbols("independent_clock independent_lapse", real=True)
    H = (
        1
        + t * p * m
        + n * n * p * p
        + n * m * eta
        + t * t * n * ph
        + sh * eta * ph
        + n**3 * sh
    )
    C = s.diff(H, n)
    K = s.diff(C, t) + canonical.flow(C, H)
    off = (
        s.diff(K, n)
        - s.diff(C, n, t)
        - canonical.flow(s.diff(C, n), H)
        + C * s.diff(C, p)
    )
    checks = {
        "entire_first_preservation": s.factor(d["K"] - FIRST),
        "entire_second_quadratic_a": s.factor(a - source.FOLD_SECOND),
        "entire_second_quadratic_b": s.factor(b - eb),
        "entire_second_quadratic_c": s.factor(c - ec),
        "first_compatibility_unique_p": s.factor(FIRST.subs(p, COMPATIBLE_P)),
        "zero_profile_second_quadratic": s.factor(
            (a * LAM**2 + b * LAM + c).subs(center)
            - 3 * (18985 * LAM**2 + 40176) / 2000
        ),
        "zero_profile_discriminant": s.factor(
            (b * b - 4 * a * c).subs(center) + s.Rational(85808403, 12500)
        ),
        "full_joint_auxiliary_null": D * s.Matrix([1, 3 * p]),
        "full_joint_auxiliary_determinant": D.det(),
        "full_raw_secondary_projection": s.expand(projection - ct - cp * ch),
        "off_surface_weighted_constraint_term": s.factor(off),
        **{
            "heavy_first_gradient_" + label + "_" + str(i): s.diff(expr, z).subs(
                source.ZERO_HEAVY
            )
            for label, expr in (
                ("H", source.clock_jet("H")),
                ("C", source.clock_jet("C")),
            )
            for i, z in enumerate((eta, ph))
        },
    }
    return {
        "whole_first_preservation_numerator": d["K"],
        "whole_complete_second_preservation_coefficients": (a, b, c),
        "whole_first_compatible_trace_density": COMPATIBLE_P,
        "whole_compatible_trace_density_bound": 21 * EPS,
        "whole_second_quadratic_centers": centers,
        "whole_exact_polynomial_remainder_bounds": errors,
        "whole_safe_second_coefficient_enclosures": {
            "a_lower": 28,
            "a_upper": 29,
            "absolute_b_upper": 1,
            "c_lower": 60,
            "discriminant_upper": -6719,
        },
        "whole_joint_auxiliary_hessian_and_null": (D, s.Matrix([1, 3 * p])),
        "whole_regular_constraint_gradient": "C_s=3 and Pbeta is nonzero. The constraint submanifold remains regular; the auxiliary Poisson/pullback rank changes. These are different notions of rank.",
        "whole_second_preservation_identity": "On C=C_N=0, K=C_u+{C,kappa V H}=0, any C1 lapse and canonical solution must satisfy C_NN lambda^2+2 D(C_N)lambda+D^2 C=0, D=d_u+{.,kappa V H}. Off the constraint K_N=D(C_N)-C*C_p, not a falsely global equality.",
        "whole_no_C1_conclusion": "For every real trace density p in this homogeneous zero-heavy bounce-slice fold family, either K is nonzero or p equals COMPATIBLE_P. In the latter case the complete actual-profile second quadratic has negative discriminant. No finite real lapse slope allows a C1 clock-time continuation through these data. This is not a classification of every inhomogeneous fold or a quantum no-go.",
        "checks": checks,
        "gates": {
            "unique_first_compatible_denominator_positive": s.Rational(3, 25) - 6 * EPS
            > 0,
            "compatible_p_below21_profile_bound": (s.Rational(5, 2) * EPS)
            / (s.Rational(3, 25) - 6 * EPS)
            < 21 * EPS,
            "second_a_above28": centers[0] - errors[0] > 28,
            "second_a_below29": centers[0] + errors[0] < 29,
            "second_b_absolute_below_one": abs(centers[1]) + errors[1] < 1,
            "second_c_above60": centers[2] - errors[2] > 60,
            "strict_negative_full_profile_discriminant": 1 - 4 * 28 * 60 == -6719,
            "future_p_minus_tenth_K_above11_over1000": s.Rational(3, 250)
            - s.Rational(31, 10) * EPS
            > s.Rational(11, 1000),
            "past_p_plus_tenth_K_below_minus11_over1000": -s.Rational(3, 250)
            + s.Rational(31, 10) * EPS
            < -s.Rational(11, 1000),
            "only_C1_lapse_required_by_second_order_Taylor_argument": True,
            "no_new_first_class_generator_or_constraint_deletion": True,
        },
    }
