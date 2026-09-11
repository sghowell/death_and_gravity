"""The exact beta(8,8) polynomial interpolation gate."""

from functools import cache

import sympy as s

X = s.Symbol("X", real=True)


@cache
def data():
    t = s.Symbol("t", real=True)
    derivative = 51480 * t**7 * (1 - t) ** 7
    T = s.expand(s.integrate(derivative, (t, 0, X)))
    q = 1 - T
    expected = (
        6435 * X**8
        - 40040 * X**9
        + 108108 * X**10
        - 163800 * X**11
        + 150150 * X**12
        - 83160 * X**13
        + 25740 * X**14
        - 3432 * X**15
    )
    quotient, remainder = s.div(q, (1 - X) ** 8, X)
    vacuum_quotient, vacuum_remainder = s.div(q - 1, X**8, X)
    checks = {
        "literal_normalized_beta_polynomial": s.expand(T - expected),
        "unit_total_beta_integral": T.subs(X, 1) - 1,
        "reflection_identity": s.expand(T + T.subs(X, 1 - X) - 1),
        "clock_eighth_factor": remainder,
        "vacuum_eighth_factor": vacuum_remainder,
        "clock_leading_eighth_coefficient": quotient.subs(X, 1) - 6435,
        "vacuum_leading_eighth_coefficient": vacuum_quotient.subs(X, 0) + 6435,
        "finite_polynomial_degree": s.degree(q, X) - 15,
    }
    for j in range(8):
        checks[f"vacuum_gate_derivative_{j}"] = s.diff(q, X, j).subs(X, 0) - int(j == 0)
        checks[f"clock_gate_derivative_{j}"] = s.diff(q, X, j).subs(X, 1)
    return {
        "step": T,
        "gate": q,
        "clock_eighth_factor_quotient": quotient,
        "vacuum_eighth_factor_quotient": vacuum_quotient,
        "step_nonzero_coefficients": {j: s.expand(T).coeff(X, j) for j in range(8, 16)},
        "definition": "Fhat(Psi)=Psi+q8(X)R_old[Psi], X=(partial Psi)^2/kappa, with the SAME literal cubic R and numerical canonical reference coefficients. q8 is a finite degree-fifteen polynomial, not the rational target's separate switching function.",
        "domain": "A finite polynomial in real or complex X. No new denominator, branch point or infinite scalar-field series is introduced by this gate. Its values outside the two stated domains are not bounded or used to assert a global inverse.",
        "checks": checks,
    }
