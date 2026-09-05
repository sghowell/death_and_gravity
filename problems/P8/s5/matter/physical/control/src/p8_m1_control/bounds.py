"""Exact coefficient majorants on overlapping compact-time high-q charts.

These are deliberately conservative sufficient bounds, not optimized scales,
an interacting cutoff, or a global fixed-comoving-mode assertion.
"""

from functools import cache

import sympy as sp

from . import model as m
from . import oscillator

ZMAX = sp.Rational(1, 1000)


def denominator_bases(chart):
    bases = [(m.J, sp.Rational(1, 10), "J"), (m.J0, sp.Rational(1, 10), "J0")]
    if chart == "unitary":
        bases += [(4-(1-m.x*m.x)**3, sp.Integer(3), "Theta/x"),
                  (m.x, sp.Rational(1, 9), "abs(x)")]
    elif chart == "gamma":
        bases += [(m.LAMBDA**2-m.J*m.z, sp.Rational(4, 125), "R=Lambda^2-J*z"),
                  (m.LAMBDA**2-m.J0*m.z, sp.Rational(31, 1000), "D=Lambda^2-J0*z"),
                  (m.LAMBDA, sp.Rational(1, 5), "abs(Lambda)")]
    else:
        raise ValueError("Use unitary or gamma chart")
    return bases


def coefficient_bound(value, chart):
    even, odd = m.parity(value)
    common = even.denom.lcm(odd.denom)
    numerator = ((even.numer*common.exquo(even.denom)).as_expr()
                 +m.l*(odd.numer*common.exquo(odd.denom)).as_expr())
    denominator = common.as_expr()
    remaining = sp.Poly(denominator, m.x, m.z, domain=sp.QQ)
    lower, factors = sp.Integer(1), []
    # Strip ONLY named, already-proved nonvanishing factors. No generic
    # factorization or numerical root tolerance is used to invent a domain.
    for base, absolute_lower, name in denominator_bases(chart):
        base_poly = sp.Poly(base, m.x, m.z, domain=sp.QQ)
        power = 0
        while remaining.total_degree() >= base_poly.total_degree():
            quotient, remainder = sp.div(remaining, base_poly)
            if not remainder.is_zero:
                break
            remaining = quotient
            power += 1
        if power:
            lower *= absolute_lower**power
            factors.append({"factor": name, "power": power, "absolute_lower": str(absolute_lower)})
    if remaining.total_degree() != 0:
        raise ValueError(f"Unproved control denominator in {chart}: {remaining.as_expr()}")
    lower *= abs(remaining.as_expr())
    if lower == 0:
        raise ValueError("Zero control denominator")
    xmax = sp.Integer(1) if chart == "unitary" else sp.Rational(1, 4)
    polynomial = sp.Poly(numerator, m.x, m.z, m.l, domain=sp.QQ)
    l1 = sum(abs(value)*xmax**powers[0]*ZMAX**powers[1]*sp.Rational(1, 10)**powers[2]
             for powers, value in polynomial.terms())
    return {"absolute_bound": str(l1/lower), "weighted_numerator_L1": str(l1),
            "denominator_absolute_lower": str(lower), "denominator_factors": factors}


def covering_checks():
    # Arithmetic implications of the elementary interval proofs in the notes.
    checks = {
        "gamma_lambda_lower": sp.Rational(1933, 8192)-sp.Rational(1, 5),
        "R_lower": sp.Rational(1, 25)-8*ZMAX-sp.Rational(4, 125),
        "D_lower": sp.Rational(1, 25)-9*ZMAX-sp.Rational(31, 1000),
        "R_relative_lower": sp.Rational(1, 5)*sp.Rational(1, 25)-8*ZMAX,
        "D_relative_lower": sp.Rational(1, 4)*sp.Rational(1, 25)-9*ZMAX,
        "chart_overlap": sp.Rational(1, 4)-sp.Rational(1, 9),
        "unitary_alpha_upper": 512-(144+1+sp.Rational(3, 10)**2),
        "gamma_alpha_upper": 512-(500+sp.Rational(4, 3)*(1+sp.Rational(5, 8)**2)),
        "unitary_alpha_lower": sp.Rational(1, 80)/(sp.Rational(1, 80)+1+sp.Rational(3, 10)**2)-sp.Rational(1, 100),
        "gamma_alpha_lower": sp.Rational(4, 5)/(sp.Rational(4, 5)+sp.Rational(4, 3)*(1+sp.Rational(5, 8)**2))-sp.Rational(1, 100),
        "root_ratio_unitary": 100-80,
        "root_ratio_gamma": 100-sp.Rational(5, 3),
        "window_gamma_margin": sp.Rational(1, 4)-(sp.Rational(9, 50)+sp.Rational(1, 99)),
        "window_unitary_margin": sp.Rational(9, 50)-sp.Rational(1, 99)-sp.Rational(1, 9),
        "window_q_lower": 2*sp.Rational(19, 20)**2-1,
        "window_q_upper": 2-(sp.Rational(101, 100)*sp.Rational(99, 95))**2,
        "energy_exponent_below_one": 1-sp.Rational(36, 99),
        "energy_growth_below_two": 2-sp.Rational(11, 7),
    }
    if any(value < 0 for value in checks.values()):
        raise ValueError("A covering or norm implication failed")
    return {key: str(value) for key, value in checks.items()}


@cache
def build_bounds():
    entries, constraints = {}, [sp.Integer(1000)]
    for chart in ("unitary", "gamma"):
        data = oscillator.derive(chart)
        names = ("mass11", "mass22", "mass12_factor", "connection_factor",
                 "covariant_mass11", "covariant_mass22", "covariant_mass12_factor")
        row = {name: coefficient_bound(data[name], chart) for name in names}
        b = {key: sp.Rational(value["absolute_bound"]) for key, value in row.items()}
        C0 = max(b["mass11"], b["mass22"])+10*b["mass12_factor"]
        C1 = max(b["covariant_mass11"], b["covariant_mass22"])+10*b["covariant_mass12_factor"]
        Comega = 10*b["connection_factor"]
        entries[chart] = {"coefficient_bounds": row, "mass_operator_bound": str(C0),
                          "covariant_mass_operator_bound": str(C1), "connection_operator_bound": str(Comega)}
        constraints.extend((2*C0, C1))
    entries["tensor"] = {"mass_operator_bound": "30", "covariant_mass_operator_bound": "132",
                         "connection_operator_bound": "0"}
    constraints.extend((sp.Integer(60), sp.Integer(132)))
    threshold = sp.Integer(1000)
    while threshold < max(constraints):
        threshold *= 10
    return {"charts": entries, "q_threshold": str(threshold),
            "threshold_comparison_margins": [str(threshold-value) for value in constraints],
            "potential_eigenvalue_bounds": ["q/(2*ell^2)", "3*q/(2*ell^2)"],
            "covariant_potential_derivative_bound": "9*q/ell^3",
            "relative_free_energy_rate_bound": "18/ell",
            "local_window_half_length": "ell0/100",
            "centre_q_sufficient": f"q0 >= {2*threshold}",
            "full_window_energy_ratio": ["7/11", "11/7"],
            "not_an_interacting_cutoff_or_optimized_threshold": True}
