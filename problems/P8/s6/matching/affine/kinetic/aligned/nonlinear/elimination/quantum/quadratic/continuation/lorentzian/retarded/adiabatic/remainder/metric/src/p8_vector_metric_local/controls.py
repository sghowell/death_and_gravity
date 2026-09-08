"""Exact failure controls for premature limits and missing vertices."""
from functools import cache

import sympy as sp

from . import consistency, counterterms, jets, radial


@cache
def omitted_counterterm():
    result = {}
    for order in (0, 1, 2):
        table = {}
        for output in ("N", "Z"):
            raw = sp.expand(counterterms.mass**(4-2*order)*jets.actual(radial.laurent(output, order)["finite_MSbar_mu_m"]))
            for source, row in (("N", jets.n), ("Z", jets.v)):
                table[output+source] = tuple(sp.factor(raw.coeff(field)) for field in row[:5])
        reverse = consistency.adjoint(table["ZN"])
        defect = tuple(sp.factor(left-right) for left, right in zip(table["NZ"], reverse))
        result[order] = {"mixed_adjoint_defect_coefficients": defect,
                         "bounce_unit_mass_fixtures": tuple(sp.factor(value.subs({jets.u: 0, counterterms.mass: 1})) for value in defect)}
    return result


@cache
def second_mass_vertex():
    value = radial.laurent("N", 0)["finite_MSbar_mu_m"]
    zero = {field: 0 for row in (jets.alpha2, jets.beta2) for field in row}
    omitted = jets.actual(value-value.subs(zero))
    return {"retained_minus_omitted": sp.factor(omitted),
            "bounce_unit_source_fixture": sp.factor(omitted.subs({jets.u: 0, jets.n[0]: 1}))}


@cache
def top_derivative():
    u, r = jets.u, 1+jets.u**2
    computed = sp.ImmutableMatrix([[consistency.coefficients(output, 2, source)[4] for source in ("N", "Z")]
                                   for output in ("N", "Z")])
    expected = sp.ImmutableMatrix([[sp.Rational(1256, 6561)/r**6, -sp.Rational(76, 243)/r**3],
                                   [-sp.Rational(76, 243)/r**3, -sp.Integer(4)]])
    return {"matrix": computed, "closed_matrix_residual": sp.ImmutableMatrix(computed-expected),
            "determinant": sp.factor(computed.det()),
            "closed_determinant_residual": sp.factor(computed.det()+sp.Rational(50992, 59049)/(1+u**2)**6),
            "scope": "This finite local fourth-derivative symbol is indefinite. It is not a resummed particle spectrum or a quantum stability verdict."}


@cache
def checks():
    return {"full_metric_fourth_derivative_matrix": top_derivative()["closed_matrix_residual"],
            "full_metric_fourth_derivative_determinant": top_derivative()["closed_determinant_residual"],
            "second_mass_contact_nonzero_fixture": sp.factor(second_mass_vertex()["bounce_unit_source_fixture"]+sp.Rational(826, 243))}
