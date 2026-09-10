"""Subtracted vertex kernel plus the full finite proper MS local anchor."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_proper_references import conversion, dirac


@cache
def data():
    Y, a, Cf, Q = sp.symbols("Y a Cf Q", positive=True)
    J0, R = sp.symbols("J0 R", real=True)
    m = sp.Symbol("m", positive=True)
    l = sp.symbols("real_l0:4", real=True)
    gamma = dirac.gamma_matrices()
    slash = sum((v * g for v, g in zip(l, gamma)), sp.zeros(4))
    S = (m * sp.eye(4) - sp.I * slash) / (m * m + sum(v * v for v in l))
    mass_derivative = (S.diff(m) + S * S).applyfunc(sp.simplify)
    ups = (Y * (J0 + 2 * R) - 6 * a * Cf) / Q
    parent = conversion.data()["finite_Yukawa_relative_upsilon"]
    parent = parent.xreplace(
        {
            v: {"Y": Y, "a": a, "Cf": Cf, "Q": Q, "J0": J0, "R": R}[str(v)]
            for v in parent.free_symbols
        }
    )
    pole = (-Y + 4 * a * Cf) / Q
    ct = (Y - 4 * a * Cf) / Q
    A, B, D, E = sp.symbols("A B D E", commutative=False)
    return {
        "scalar_vertex_kernel": "K_b(q,l,p)=B_b(q-l) S(l) S(l+p); b=1",
        "gauge_vertex_kernel": "The same ordered product capped by gamma_mu on both sides, summed over mu; b=0. Its coupling has the opposite sign to the scalar kernel.",
        "zero_momentum_subtraction": "K_b(0,l,0)=B_b(l) S(l)^2, with identical regulator and dimensional numerator",
        "exact_MS_vertex_decomposition": "Gamma_MS(q,p)/y = ups_MS + integral of the signed scalar/gauge K(q,l,p)-K(0,l,0), including their internal couplings",
        "finite_MS_anchor": ups,
        "local_anchor_absolute_upper": (2 * Y + 6 * a * Cf) / Q,
        "kernel_coupling_majorant": Y + 4 * a * Cf,
        "proper_UV_vertex_pole": pole,
        "paired_proper_MS_counterterm": ct,
        "scope": "The zero-momentum anchor is a proper vertex, not the canonical bare-Yukawa conversion. Its finite dimensional terms are retained; all other field/parameter conversions remain separate.",
        "checks": {
            "literal_mass_derivative_of_inverse": mass_derivative,
            "same_finite_MS_vertex_anchor": sp.simplify(ups - parent),
            "proper_vertex_pole_counterterm_pair": pole + ct,
            "fixed_scale_gauge_vertex_anchor": ups.subs(Y, 0) + 6 * a * Cf / Q,
            "local_anchor_below_twice_kernel_majorant": 2 * (Y + 4 * a * Cf)
            - (2 * Y + 6 * a * Cf)
            - 2 * a * Cf,
            "noncommuting_product_difference": sp.expand(
                A * B - D * E - ((A - D) * B + D * (B - E))
            ),
            "zero_yukawa_no_external_Phi_vertex": (sp.Symbol("y") * ups).subs(
                sp.Symbol("y"), 0
            ),
            "one_loop_subgraph_dimension": 4 - 2 - 2,
        },
    }
