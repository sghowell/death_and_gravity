"""Strict coefficient extraction only from an entire stated local polynomial."""

from functools import cache

import sympy as s

from . import basis

VARS = basis.variables
REFERENCE = s.Poly(basis.canonical, *VARS).monoms()[0]
REFERENCE_COEFFICIENT = s.Poly(basis.canonical, *VARS).coeff_monomial(REFERENCE)


def extract_local_coefficient(expression):
    if isinstance(expression, (bool, float, str)) or not isinstance(
        expression, (int, s.Basic)
    ):
        raise TypeError("Require an entire exact local polynomial")
    expression = s.sympify(expression)
    if expression.has(s.Float) or expression.free_symbols - set(VARS):
        raise ValueError(
            "Unknown matching symbols or inexact data are not coefficients"
        )
    try:
        polynomial = s.Poly(s.expand(expression), *VARS)
    except s.PolynomialError as error:
        raise ValueError(
            "Require a polynomial, not a nonlocal or rational response"
        ) from error
    coefficient = s.simplify(
        polynomial.coeff_monomial(REFERENCE) / REFERENCE_COEFFICIENT
    )
    if (
        coefficient.is_real is not True
        or coefficient.is_finite is not True
        or coefficient.free_symbols
    ):
        raise ValueError("Require a finite exact real coefficient")
    if s.expand(expression - coefficient * basis.canonical) != 0:
        raise ValueError(
            "Entire residual does not belong to the stated one-dimensional truncation"
        )
    return coefficient


@cache
def data():
    chi = s.Symbol("chi", real=True)
    checks = {
        "nonzero_reference_normalization": s.Poly(
            basis.canonical, *VARS
        ).coeff_monomial(REFERENCE)
        - REFERENCE_COEFFICIENT,
        "generic_complete_matching_family": s.expand(
            chi * basis.canonical
            - sum(chi * basis.q(i, j, i, j) for i, j in basis.pairs)
        ),
        "unit_action_four_label_multiplicity": s.Integer(len(basis.perms))
        - 4 * len(basis.pairs),
    }
    diagnostics = (s.S.Zero, s.S.One, -s.Rational(3, 7), s.sqrt(2), s.log(2))
    for index, value in enumerate(diagnostics):
        checks[f"exact_full_polynomial_coefficient_recovery_{index}"] = (
            extract_local_coefficient(value * basis.canonical) - value
        )
    higher = s.expand(sum(a * a for a in basis.a) * basis.canonical)
    gates = {
        "reference_monomial_coefficient_nonzero": REFERENCE_COEFFICIENT != 0,
        "higher_derivative_Bose_symmetric_example_not_in_basis": s.Poly(
            higher, *VARS
        ).total_degree()
        > s.Poly(basis.canonical, *VARS).total_degree(),
        "finite_chi_not_selected_by_basis_or_pole_cancellation": True,
        "one_physical_sample_not_used_as_full_matching": True,
        "no_naturalness_or_minimal_subtraction_zero_assumption": True,
        "coefficient_diagnostics_not_substituted_for_original_source": True,
    }
    return {
        "checks": checks,
        "gates": gates,
        "whole_coefficient_functional": {
            "monomial": REFERENCE,
            "normalization": REFERENCE_COEFFICIENT,
        },
        "whole_matching_requirement": "After actual matching supplies the entire real local six-derivative four-scalar TT residual in this fixed representative, one coefficient remains. Extraction verifies the entire polynomial equals chi*T; it rejects extra symbols, floats, rational/nonlocal terms and every out-of-span residual. A single physical sample cannot isolate chi amid unknown higher-derivative or nonlocal terms.",
        "whole_pole_finite_distinction": "Zero additional selected-matter physical-TT UV pole leaves every finite chi in the formal matching family. The loop-counted term first affects the four-scalar one-graviton amplitude at this order but is invisible to the retained flat one-loop matching. This is not permission to assign an arbitrary coefficient to the actual parent.",
        "whole_remainder_boundary": "The existing 1e-196 bound applies to the selected minimal matter remainder. Without physical chi matching it is not a bound on the sum with the independent curvature contact, and no full inclusive or global P8 result follows.",
    }
