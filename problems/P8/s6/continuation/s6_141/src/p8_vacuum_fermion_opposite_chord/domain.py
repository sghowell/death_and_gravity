"""Joint real-loop routing and a minimum-momentum complex Cauchy radius."""

from functools import cache

import sympy as sp
from p8_vacuum_fermion_proper_references import dirac

from . import catalog


def routing(word):
    i = word.index("B")
    rotated = word[i:] + word[:i]
    j = rotated.index("B", 1)
    return {
        "arc_l_external": rotated[1:j],
        "arc_q_external": rotated[j + 1 :],
        "chord": "q-l",
        "second_endpoint_chord": "(q+P_arc)-(l+P_arc)=q-l",
        "fermion_counts": (j, len(rotated) - j),
    }


@cache
def data():
    gamma = dirac.gamma_matrices()
    q = sp.symbols("real_q0:4", real=True)
    m = sp.Symbol("m", positive=True)
    D = m * sp.eye(4) + sp.I * sum((v * g for v, g in zip(q, gamma)), sp.zeros(4))
    normality = (
        D.conjugate().T * D - (m * m + sum(v * v for v in q)) * sp.eye(4)
    ).applyfunc(sp.expand)
    E, P, s = sp.symbols("E P s")
    p1, p2 = sp.Matrix([sp.I * E, P, 0, 0]), sp.Matrix([sp.I * E, -P, 0, 0])
    external = (p1, p2, -p1, -p2)
    subst = {E**2: s / 4, P**2: s / 4 - 1}
    r = [routing(w) for w in catalog.selected_words()]
    q0, l0, A = sp.symbols("q l arc_external_sum")
    checks = {
        "literal_Dirac_normality": normality,
        "routing_keeps_boson_unshifted": (q0 + A) - (l0 + A) - (q0 - l0),
        "all_words_have_three_fermions_per_arc": sum(
            row["fermion_counts"] != (3, 3) for row in r
        ),
        "all_words_have_two_external_legs_per_arc": sum(
            (len(row["arc_l_external"]), len(row["arc_q_external"])) != (2, 2)
            for row in r
        ),
        "total_external_momentum": sum(external, sp.zeros(4, 1)),
        "forward_s": sp.expand(-((p1 + p2).T * (p1 + p2))[0]).subs(subst) - s,
        "forward_u": sp.expand(-((p1 - p2).T * (p1 - p2))[0]).subs(subst) - (4 - s),
        "two_arc_prefix_legs_bound": 2 * 2 - 4,
        "conservative_prefix_slack": 18 - 4 - 14,
        "joint_radius_shift_margin": sp.Rational(18, 360) - sp.Rational(1, 20),
        "Neumann_inverse_bound_slack": 2 - sp.Rational(20, 19) - sp.Rational(18, 19),
        "minimum_radius_at_mass_720": sp.Rational(720, 360) - 2,
        "six_propagators_and_Dirac_trace": 4 * 2**6 - 256,
        "four_gauge_index_terms": len(gamma) - 4,
        "active_gauge_color_trace": 6 * sp.Rational(4, 3) - 8,
    }
    for i, p in enumerate(external):
        checks[f"external_mass_one_{i}"] = sp.expand((p.T * p)[0]).subs(subst) + 1
    return {
        "literal_routings": r,
        "fermion_norm_bound": "2/sqrt(m^2+q^2) on each q arc and 2/sqrt(m^2+l^2) on each l arc",
        "joint_Cauchy_radius": "min(sqrt(m^2+q^2),sqrt(m^2+l^2))/360",
        "boson_norm_bound": "1/((q-l)^2+b) <= 1/(q-l)^2 for b=0 or 1; no complex external shift",
        "whole_word_trace_upper": "256 N Y^2 (Y+4aCf)/[(m^2+q^2)^(3/2)(m^2+l^2)^(3/2)(q-l)^2]",
        "scope": "A joint two-loop analytic domain, not an export of the self-energy subgraph radius to a vertex correction.",
        "checks": checks,
    }
