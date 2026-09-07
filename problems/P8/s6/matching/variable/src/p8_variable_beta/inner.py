"""Exact shrinking-inner-window tensor operator and nonadiabaticity audit.

Set c=2+epsilon², u=epsilon*x, epsilon>0. All limits hold on each fixed
compact x interval at fixed finite kbar. This is neither a regular finite-
coefficient parent at epsilon=0 nor a fixed physical-time-window theorem.
"""

from functools import cache

import sympy as sp

from . import canonical

EPS = sp.Symbol("epsilon", positive=True)
X = sp.Symbol("x", real=True)


def at_inner_limit(expression):
    d = canonical.derive()
    value = expression.subs({d["c"]: 2+EPS**2, d["u"]: EPS*X}, simultaneous=True)
    # Rational expressions have removable epsilon factors. Radical rotation
    # factors sqrt(2+epsilon²) are analytic on the selected positive branch.
    return sp.simplify(sp.limit(sp.factor(value), EPS, 0, dir="+"))


@cache
def checks():
    d = canonical.derive()
    q = d["mass_squared"]
    u, c = d["u"], d["c"]
    ratio = sp.factor((sp.diff(q, u, 2)/q**2).subs(u, 0))
    expected = -c*(7*c**2+146*c-240)/(8*(c+8)**2)
    cross = d["q_squared"]*d["D"]
    # Derivatives are taken in u first; only then transform the entire
    # second-order equation to x and multiply it by epsilon².
    coefficients = {
        "inner_heavy_potential": EPS**2*d["V_HH"],
        "inner_light_potential": EPS**2*d["V_LL"],
        "inner_heavy_first_derivative_cross": 2*EPS*d["omega"],
        "inner_light_first_derivative_cross": -2*EPS*d["omega"],
        "inner_heavy_zero_order_cross": EPS**2*(cross-2*d["omega"]*d["theta_sum"]),
        "inner_light_zero_order_cross": EPS**2*(cross-2*sp.diff(d["omega"], u)-2*d["omega"]*d["theta_sum"]),
    }
    values = {key: at_inner_limit(value) for key, value in coefficients.items()}
    values["inner_heavy_potential"] -= 80/(1+8*X**2)
    values.update({"exact_mass_squared_variation_ratio": sp.factor(ratio-expected),
                   "nonvanishing_mass_squared_variation_limit": sp.limit(ratio, c, 2, dir="+")+sp.Rational(1, 5),
                   "nonvanishing_frequency_variation_limit": sp.limit(ratio/2, c, 2, dir="+")+sp.Rational(1, 10),
                   "divergent_center_mass": sp.limit((c-2)*q.subs(u, 0), c, 2, dir="+")-80,
                   "locked_cone_excess_same_order": sp.limit((d["c_light_squared"].subs(u, 0)-1)/(c-2), c, 2, dir="+")-sp.Rational(4, 5)})
    return {key: sp.simplify(value) for key, value in values.items()}


def canonical_adiabatic_checks():
    d = canonical.derive()
    v, u, c = d["V_HH"], d["u"], d["c"]
    # The exact canonical corrections are analytic at c=2,u=0. Verify
    # their finite value and finite second derivative before using the
    # pole orders of the algebraic mass to transfer the ratio limits.
    correction = v-d["mass_squared"]
    finite = {"correction_at_center_limit": sp.simplify(correction.subs({u: 0, c: 2})),
              "correction_second_at_center_limit": sp.simplify(sp.diff(correction, u, 2).subs({u: 0, c: 2}))}
    if any(value.has(sp.oo, sp.zoo, sp.nan) for value in finite.values()):
        raise ValueError("A canonical normalization correction is singular in the inner chart")
    return finite


def controls():
    mass = 80/(1+8*X**2)
    exact = (1+8*X**2)/96
    algebraic = (1+8*X**2)/80
    residual = sp.simplify(sp.diff(algebraic, X, 2)+mass*algebraic-1)
    if sp.simplify(sp.diff(exact, X, 2)+mass*exact-1) != 0 or residual != sp.Rational(1, 5):
        raise ValueError("The limiting forced-operator omission control failed")
    return {"limiting_exact_particular": exact, "limiting_algebraic_representative": algebraic,
            "omitted_derivative_residual": residual, "particular_to_algebraic_ratio": sp.cancel(exact/algebraic),
            "source_scope": "Limiting ODE constant forcing only; not a specified physical-metric source or retarded matching exclusion"}
