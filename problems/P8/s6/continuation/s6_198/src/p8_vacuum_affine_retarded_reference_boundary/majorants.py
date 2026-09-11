"""Exact Banach-valued Cauchy and full Leibniz majorants through order six."""

from functools import cache

import sympy as s
from p8_vacuum_affine_cd_metric_noise import reference

from . import ibp

RADIUS = reference.RADIUS
PAIR = reference.PAIR_REF


@cache
def recurrence():
    rows = [(s.Integer(1),)]
    for n in range(6):
        old = rows[-1]
        rows.append(
            tuple(
                (2 * n + 2 - j) * (old[j] if j < len(old) else 0)
                + (old[j - 1] if j else 0)
                for j in range(n + 2)
            )
        )
    return tuple(rows)


@cache
def data():
    C, G, r = s.symbols("C G r", positive=True)
    replacements = {ibp.A[j]: s.factorial(j) * C / r**j for j in range(8)}
    replacements.update({ibp.G[j]: s.factorial(j) * G / r**j for j in range(8)})
    checks = {}
    for n, expr in enumerate(ibp.iterates()):
        for j, value in enumerate(recurrence()[n]):
            raw = expr.coeff(ibp.H[j]).subs(replacements, simultaneous=True)
            checks[f"full_Leibniz_Cauchy_order_{n}_test_jet_{j}"] = s.cancel(
                raw / (C * G**n / r ** (n - j)) - value
            )
    return {
        "actual_analytic_input": "The complete phase-stripped reference pair amplitude has norm<=1e9 sqrt(nu mu) on the S186 radius1/20000 complex-time discs; g=1/(W_k+W_l) has modulus<=2/(nu+mu). Both constant actual-alpha factors and the full constrained ten-field readout are retained.",
        "coefficient_recurrence": "c[0,0]=1; c[n+1,j]=(2n+2-j)c[n,j]+c[n,j-1]. A term has n+1 inverse-phase factors after multiplication by g. Differentiating these and the amplitude contributes2n+2-j to the Cauchy factorial bound; differentiating the test contributes the previous coefficient.",
        "complete_coefficients": recurrence(),
        "sixth_sum": sum(recurrence()[6]),
        "fifth_sum": sum(recurrence()[5]),
        "bound": "||L^n(a h)|| <= Cr sqrt(nu mu) [2/(nu+mu)]^n sum_(j=0)^n c[n,j] r^(-(n-j)) ||partial_t^j h||, for n<=6. This is the complete product rule, not a slow-mixing derivative assumption.",
        "checks": checks,
        "gates": {
            "third_order_agrees_with_complete_old_coefficients": recurrence()[3]
            == (48, 33, 9, 1),
            "complete_sixth_coefficients": recurrence()[6]
            == (46080, 35685, 12645, 2640, 345, 27, 1),
            "complete_sixth_sum": sum(recurrence()[6]) == 97423,
            "complete_fifth_sum": sum(recurrence()[5]) == 7916,
            "actual_fixed_complex_disc": RADIUS == s.Rational(1, 20000),
        },
    }
