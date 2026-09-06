"""Exact causal principal-log resolvent including the stiff Einstein term.

All three poles and the positive cut are retained. The associated constant-
coefficient linear calculation does not assert a runaway of the actual SEE.
"""

import sympy as sp
from p8a_existence.mode_lipschitz import rational


def positive(value):
    value = rational(value)
    if value <= 0:
        raise ValueError("an exact positive rational is required")
    return value


def named_parameters(scale="5/2", delta="1/100000000000000"):
    """The freezing scale is a preconditioner choice, not a new physical model."""
    a, delta = map(positive, (scale, delta))
    beta = sp.EulerGamma-sp.Rational(19, 30)-sp.log(a/2)
    c = a**2/(30*delta)
    return {"scale": a, "delta": delta, "beta": beta, "c": c,
            "Lambert_argument": 4*sp.exp(2*sp.EulerGamma-sp.Rational(19, 15))/(15*delta)}


def transfer(laplace, beta, coefficient):
    """Symbolic transform 2/(log(s)+beta-c/s^2), principal log, c>0."""
    s, beta, c = map(sp.sympify, (laplace, beta, coefficient))
    return 2*s**2/(s**2*(sp.log(s)+beta)-c)


def pole_expressions(beta, coefficient):
    """Three principal-s-sheet poles, with beta real and coefficient positive.

The explicit sheet/count proof is in notes/resolvent.md; no floating-point
root search is used to infer completeness of this list.
    """
    beta, c = map(sp.sympify, (beta, coefficient))
    argument = 2*c*sp.exp(2*beta)
    return {branch: sp.exp(sp.LambertW(argument, branch)/2-beta) for branch in (0, 1, -1)}


def residue(pole, beta):
    p, beta = map(sp.sympify, (pole, beta))
    return 2*p/(1+2*(sp.log(p)+beta))


def cut_density(radius, beta, coefficient):
    r, beta, c = map(sp.sympify, (radius, beta, coefficient))
    return 2/((sp.log(r)+beta-c/r**2)**2+sp.pi**2)


def identities():
    s, c, af, delta, t0 = sp.symbols("s c af delta T0", positive=True)
    beta, gamma, df, aa = sp.symbols("beta gamma df aa", real=True)
    theta = sp.Symbol("theta", real=True)
    w, z = sp.symbols("w z", complex=True)
    den = s**2*(sp.log(s)+beta)-c
    lhs_bank, rhs_bank = 2/(aa-sp.I*sp.pi), 2/(aa+sp.I*sp.pi)
    angular_derivative = 1/theta-4*sp.cot(2*theta)+4*theta/sp.sin(2*theta)**2
    angular_positive = ((sp.sin(2*theta)-2*theta*sp.cos(2*theta))**2
                        +4*theta**2*sp.sin(2*theta)**2)/(theta*sp.sin(2*theta)**2)
    named_beta = sp.EulerGamma-sp.Rational(19, 30)-sp.log(af/2)
    lambert_argument = 4*sp.exp(2*sp.EulerGamma-sp.Rational(19, 15))/(15*delta)
    formal_pole_residual = sp.exp(w-2*beta)*w/2-c
    formal_lambert_residual = sp.exp(-2*beta)*(w*sp.exp(w)-z)/2
    return {
        "full_normalized_log_and_Einstein_multiplier": sp.simplify(
            transfer(s, gamma+2*df, c)-1/((sp.log(s)+gamma)/2+df-c/(2*s**2))),
        "log_only_limit": sp.simplify(transfer(s, beta, 0)-2/(sp.log(s)+beta)),
        "simple_pole_residue": sp.simplify(2*s**2/sp.diff(den, s)-residue(s, beta)),
        "positive_cut_jump": sp.simplify(lhs_bank-rhs_bank-4*sp.I*sp.pi/(aa**2+sp.pi**2)),
        "Lambert_root_relation": sp.simplify(
            (formal_pole_residual-formal_lambert_residual).subs(z, 2*c*sp.exp(2*beta))),
        "angular_count_strict_monotonicity": sp.trigsimp(angular_derivative-angular_positive),
        "nonreal_root_real_equation": sp.trigsimp(
            (-theta*sp.cot(2*theta))*sp.cos(2*theta)-theta*sp.sin(2*theta)
            +theta/sp.sin(2*theta)),
        "named_argument_scale_cancellation": sp.simplify(
            2*af**2/(30*delta)*sp.exp(2*named_beta)-lambert_argument),
        "proper_growth_rate_units": sp.simplify(
            (2*(af**2/(30*delta))/w)/(af**2*t0**2)
            -192*sp.pi**2/(2880*sp.pi**2*delta*t0**2*w)),
        "positive_Neumann_resummation": sp.simplify(
            (2/(sp.log(s)+beta))/(1-c/(s**2*(sp.log(s)+beta)))-transfer(s, beta, c)),
        "Einstein_sign_at_log_pole": sp.simplify(
            (s**2*(sp.log(s)+beta)-c).subs(s, sp.exp(-beta))+c),
    }


def controls():
    s, p, weight = sp.symbols("s p weight", positive=True)
    return {"omitted_positive_pole_Laplace_error": weight/(s-p),
            "omitting_the_two_damped_poles_is_also_wrong": True,
            "full_rolling_Frechet_derivative_computed": False,
            "causal_zero_past_homogeneous_solution_is_only_zero": True,
            "forcing_projection_not_fixed_by_the_energy_constraint_alone": True}
