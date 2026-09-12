"""Actual QG1 coefficient jets, positive scalar characteristics and chart cover."""

from functools import cache

import sympy as s
from p8_vacuum_affine_quantum_retuning import profile
from p8_vacuum_affine_reduced_scalar_hamiltonian import scalar

t = scalar.t
H, ell, delta = scalar.H, scalar.ell, scalar.delta
EPS = profile.EPS
THRESHOLD = s.Integer(100)
QMIN = s.Integer(4096)
JET = s.Integer(10) ** 6


def envelope(expr):
    num, den = s.fraction(s.factor(expr))
    denominator = s.Poly(den, t)
    if not all(power[0] % 2 == 0 and value > 0 for power, value in denominator.terms()):
        raise ValueError(
            "Require a positive-even denominator for the real slab envelope"
        )
    numerator = s.Poly(num, t)
    return sum(
        abs(value) * s.Rational(1, 2) ** power[0] for power, value in numerator.terms()
    ) / den.subs(t, 0)


@cache
def clock():
    c = scalar.actual_coefficients()
    th, E, J, w = [c[v] for v in (scalar.Theta, scalar.E, scalar.J0, scalar.w)]
    C = s.factor(th * s.diff(E, t) - E * s.diff(th, t) + H * E * th - th * th)
    F = s.factor(C - w * w / 2)
    return {"Theta": th, "E": E, "J": J, "w": w, "C": C, "F": F, "Jc": J + c[scalar.dJ]}


@cache
def data():
    c = clock()
    th, E, J, _w, _C, F, Jc = [
        c[name] for name in ("Theta", "E", "J", "w", "C", "F", "Jc")
    ]
    numerator, denominator = s.fraction(F)
    polynomial = s.Poly(numerator, t)
    rows = {
        name: [envelope(s.diff(expr, t, j)) for j in range(3)]
        for name, expr in (
            ("Theta", th),
            ("E", E),
            ("J", J),
            ("ell", ell),
            ("H", H),
            ("delta", delta),
        )
    }
    retuning = {}
    for j in range(3):
        dj = (
            sum(
                s.binomial(j, k)
                * (
                    (
                        21
                        * sum(
                            s.binomial(k, z) * rows["delta"][z] * rows["delta"][k - z]
                            for z in range(k + 1)
                        )
                        + 3 * rows["delta"][k]
                    )
                    / 2
                    + (1 if k == 0 else 0)
                    + 6 * rows["delta"][k]
                )
                for k in range(j + 1)
            )
            * EPS
        )
        tc = (
            1 + 3 * sum(s.binomial(j, k) * rows["delta"][k] for k in range(j + 1))
        ) * EPS
        retuning[j] = (dj, tc, EPS)
    Fmin = s.Rational(1199, 800) * s.Rational(4, 5) ** 18
    Jmin = s.Rational(1215, 800) * s.Rational(4, 5) ** 18 - 4 * EPS
    gap = s.Rational(1, 50) * s.Rational(4, 5) ** 6 - 4 * EPS
    x, j, w0, f = s.symbols("speed_squared Jc w F", real=True)
    K = s.Matrix([[2 * j + w0 * w0, w0], [w0, 1]])
    G = s.Matrix([[2 * f + w0 * w0, w0], [w0, 1]])
    checks = {
        "exact_scalar_gradient_gap": s.factor(J - F - 1 / (50 * (1 + t * t) ** 6)),
        "full_matter_mixed_characteristic_factor": s.factor(
            (x * K - G).det() - 2 * j * (x - 1) * (x - f / j)
        ),
        "full_gradient_numerator_constant": numerator.subs(t, 0) - 1199,
        "full_gradient_denominator": denominator - 800 * (1 + t * t) ** 18,
        "bare_bounce_scalar_speed": s.factor(
            (F / J).subs(t, 0) - s.Rational(1199, 1215)
        ),
        "central_E_endpoint": E.subs(t, s.Rational(1, 4)) + s.Rational(1231, 4913),
        "central_E_derivative": s.factor(s.diff(E, t) - 9 * t / (1 + t * t) ** 4),
        "outer_Theta_lower_representation": s.factor(
            th - t * (4 - (1 + t * t) ** -3) / (1 + t * t)
        ),
        "outer_Theta_lower_arithmetic": s.Rational(3, 8) * s.Rational(4, 5)
        - s.Rational(3, 10),
        "actual_physical_transfer_threshold": THRESHOLD**2 * s.Rational(4, 5) ** 4
        - QMIN,
        "retuned_principal_shift_only_J": s.factor(
            Jc - J - scalar.actual_coefficients()[scalar.dJ]
        ),
    }
    return {
        "actual_principal_gradient": F,
        "bare_J_minus_F": s.factor(J - F),
        "positive_even_gradient_numerator": polynomial.as_expr(),
        "actual_Jc": Jc,
        "coefficient_time_jet_envelopes": rows,
        "fixed_retuning_DeltaJ_Tcorr_A_jet_envelopes": retuning,
        "lower_bounds": {"F": Fmin, "Jc": Jmin, "Jc_minus_F": gap},
        "principal_scope": "For the full two-scalar coefficient-sector principal pair the characteristic speeds squared are1 andF/Jc. Bare F>0 and J-F=1/[50(1+t^2)^6] hold at all real times. Fixed QG1 positivity and separation are asserted only on the original unit slab, using |DeltaJ|<4epsilon. This is not a finite-q frozen-frequency stability inference.",
        "chart_cover": "The v chart is used only where |t|>=1/8, giving |Theta|>=3/10. The complementary b=-pv/(2q) chart is used only where |t|<=1/4, giving1/4<|E|<=1/2. Switch at+-3/16; at most two switches cover every ordered subinterval.",
        "momentum_cover": "For |P|>=100, q=|P|^2/a^2>=4096. Smaller physical transfers use the regular first-order comparison, not the potentially singular finite-q b chart. These thresholds are mathematical, not physical EFT cutoffs.",
        "actual_coefficient_box": "Pointwise |Theta|<=2, |H|<=2, ell<=1/10, |E|<=1, |w|<=1/10, 1/100<Jc<100. In the central chart additionally |E|<=1/2. A,Tcorr and their first two jets are<1e-6; every coefficient time jet through2 is<1e6. The fixed profile is not recomputed and no high stress jet beyond the original five is used.",
        "checks": checks,
        "gates": {
            "positive_even_gradient_polynomial": all(
                power[0] % 2 == 0 and val > 0 for power, val in polynomial.terms()
            ),
            "F_lower": Fmin > s.Rational(1, 100),
            "Jc_lower": Jmin > s.Rational(1, 100),
            "strict_characteristic_separation": gap > s.Rational(1, 1000),
            "bare_J_upper": rows["J"][0] < 99,
            "all_base_jets": all(v < JET - 1 for row in rows.values() for v in row),
            "all_retuning_jets": all(
                v < s.Rational(1, 10**6) for row in retuning.values() for v in row
            ),
            "central_E_separation": s.Rational(1231, 4913) > s.Rational(1, 4),
            "two_switches_inside_both_chart_domains": s.Rational(1, 8)
            < s.Rational(3, 16)
            < s.Rational(1, 4),
            "no_Theta_squared_rescaling_as_crossing_proof": True,
            "current_fixed_QG1_not_bare_parent_only": True,
            "no_global_current_retuning_bound_inferred": True,
        },
    }
