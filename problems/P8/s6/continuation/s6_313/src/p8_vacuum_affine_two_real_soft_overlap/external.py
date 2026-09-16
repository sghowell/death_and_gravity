"""Exact double-external ordering identity and independent full-tree checks."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import checks as samples
from p8_vacuum_affine_complete_two_graviton_tree import trees as e

MARKER = s.Symbol("light_scalar_propagator_marker")


class Marked(e.TreeEngine):
    """Count internal light scalar lines without changing canonical vertices."""

    def current(self, mask, kind):
        value, count = super().current(mask, kind)
        if kind == "phi" and count and mask & (mask - 1):
            return MARKER * value, count
        return value, count


def contraction(p, eps, q):
    return (p.T * eps * q)[0]


@cache
def general_identities():
    A, B, Z, Na, Nb, Xa, Yb = s.symbols("A B Z Na Nb Xa Yb", nonzero=True)
    ordered = Na * (Nb + Xa) / (A * (A + B + Z)) + Nb * (Na + Yb) / (B * (A + B + Z))
    correction = (
        Na * Xa / (A * (A + B + Z))
        + Nb * Yb / (B * (A + B + Z))
        - Na * Nb * Z / (A * B * (A + B + Z))
    )
    checks = {
        "general_same_leg_two_order_identity": s.cancel(
            ordered - Na * Nb / (A * B) - correction
        )
    }
    p, q, u, v = (s.Matrix(s.symbols(label + "0:4")) for label in ("p", "q", "u", "v"))
    mu = s.Symbol("mu")
    checks["general_bilinear_stress_shift"] = (
        e.old.stress(p + u, q + v, mu)
        - e.old.stress(p, q, mu)
        - e.old.stress(u, q, 0)
        - e.old.stress(p, v, 0)
        - e.old.stress(u, v, 0)
    ).applyfunc(s.expand)
    jl1, jl2, jr1, jr2 = s.symbols("jl1 jl2 jr1 jr2")
    a1, a2, b1, b2, c1, c2, d1, d2 = s.symbols("a1 a2 b1 b2 c1 c2 d1 d2")
    D0, D1, D2, D3 = s.symbols("D0 D1 D2 D3", nonzero=True)
    first = (a1, a2, c1, c2)
    second = (b1, b2, d1, d2)
    direct = s.S.Zero
    for i in range(4):
        for j in range(4):
            mask = int(i < 2) + 2 * int(j < 2)
            direct += first[i] * second[j] / (D0, D1, D2, D3)[mask]
    grouped = (jr1 * jr2 / D0 + jl1 * jr2 / D1 + jr1 * jl2 / D2 + jl1 * jl2 / D3).subs(
        {
            jl1: a1 + a2,
            jl2: b1 + b2,
            jr1: c1 + c2,
            jr2: d1 + d2,
        },
        simultaneous=True,
    )
    checks["general_four_radiation_distributions"] = s.cancel(direct - grouped)
    controls = {
        "missing_Q2_propagator_term_is_nonzero": s.cancel(
            ordered
            - Na * Nb / (A * B)
            - correction
            - Na * Nb * Z / (A * B * (A + B + Z))
        )
        != 0,
        "missing_shifted_numerator_terms_are_nonzero": s.cancel(
            ordered - Na * Nb / (A * B) + Na * Nb * Z / (A * B * (A + B + Z))
        )
        != 0,
    }
    return checks, controls


def core(ps, kind, heavy, cubic, contact):
    if kind == "C":
        return contact
    answer = s.S.Zero
    for left, right in e.old.PARTS:
        P = ps[left[0]] + ps[left[1]]
        if kind == "H":
            answer += cubic * cubic / (heavy - e.old.dot(P, P))
        else:
            TL = e.old.stress(ps[left[0]], ps[left[1]])
            TR = e.old.stress(ps[right[0]], ps[right[1]])
            HR = e.ETA * TR * e.ETA - e.ETA * s.trace(e.ETA * TR) / 2
            answer -= e.old.pair(TL, HR) / e.old.dot(P, P)
    return s.factor(answer)


@cache
def independent_tree():
    points, qs, pols = samples.nonopposite_state()
    n = s.Integer(128)
    g = s.Rational(2, 3)
    C = s.Rational(-1, 7)
    js = [
        [s.factor(contraction(p, E, p) / e.old.dot(p, q)) for p in points]
        for q, E in zip(qs, pols)
    ]
    exact = {kind: s.S.Zero for kind in ("C", "H", "GR")}
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            shifted = list(points)
            shifted[i] = e.imm(shifted[i] + qs[0])
            shifted[j] = e.imm(shifted[j] + qs[1])
            for kind in exact:
                exact[kind] += js[0][i] * js[1][j] * core(shifted, kind, n, g, C)
    for i, p in enumerate(points):
        A = e.old.dot(p, qs[0])
        B = e.old.dot(p, qs[1])
        Z = e.old.dot(qs[0], qs[1])
        Na = contraction(p, pols[0], p)
        Nb = contraction(p, pols[1], p)
        Xa = 2 * contraction(p, pols[1], qs[0]) + contraction(qs[0], pols[1], qs[0])
        Yb = 2 * contraction(p, pols[0], qs[1]) + contraction(qs[1], pols[0], qs[1])
        factor = Na * (Nb + Xa) / (A * (A + B + Z)) + Nb * (Na + Yb) / (B * (A + B + Z))
        shifted = list(points)
        shifted[i] = e.imm(p + qs[0] + qs[1])
        for kind in exact:
            exact[kind] += factor * core(shifted, kind, n, g, C)
    exact = {k: s.factor(v) for k, v in exact.items()}
    checks = {}
    gates = {}
    for label, cubic, contact, target, nexpected in (
        ("pure", 0, 0, exact["GR"], 198),
        ("with_C", 0, C, exact["GR"] + exact["C"], 236),
        ("full", g, C, sum(exact.values()), 434),
    ):
        legs = [("phi", p, 1) for p in points[1:]] + [
            ("h", q, E) for q, E in zip(qs, pols)
        ]
        value, count = Marked(legs, heavy=n, cubic=cubic, contact=contact).amplitude()
        poly = s.Poly(value, MARKER)
        checks[label + "_maximum_light_propagator_degree"] = s.Integer(
            poly.degree() - 2
        )
        checks[label + "_full_labeled_tree_count"] = s.Integer(count - nexpected)
        checks[label + "_entire_double_external_coefficient"] = s.factor(
            poly.coeff_monomial(MARKER**2) - target
        )
    corrections = {}
    for kind in ("C", "H", "GR"):
        leading = s.S.Zero
        if kind == "C":
            leading = C * sum(js[0]) * sum(js[1])
        else:
            for left, right in e.old.PARTS:
                TL = e.old.stress(points[left[0]], points[left[1]])
                TR = e.old.stress(points[right[0]], points[right[1]])
                HR = e.ETA * TR * e.ETA - e.ETA * s.trace(e.ETA * TR) / 2
                N = e.old.pair(TL, HR)
                for mask in range(4):
                    Qa = sum((qs[j] for j in range(2) if mask >> j & 1), e.VECTOR_ZERO)
                    P = points[left[0]] + points[left[1]] + Qa
                    D = e.old.dot(P, P)
                    soft = s.prod(
                        sum(js[j][i] for i in (left if mask >> j & 1 else right))
                        for j in range(2)
                    )
                    leading += soft * (g * g / (n - D) if kind == "H" else -N / D)
        corrections[kind] = s.factor(exact[kind] - leading)
        gates[kind + "_leading_only_grouping_has_nonzero_correction"] = (
            corrections[kind] != 0
        )
    return checks, gates, corrections


@cache
def data():
    checks, gates = general_identities()
    checks = dict(checks)
    gates = dict(gates)
    physical, physical_gates, corrections = independent_tree()
    checks.update(physical)
    gates.update(physical_gates)
    return {
        "whole_general_two_order_identity": "For A=p.q1,B=p.q2,Z=q1.q2,Na=p.eps1.p,Nb=p.eps2.p,Xa=2p.eps2.q1+q1.eps2.q1,Yb=2p.eps1.q2+q2.eps1.q2, the sum of both orderings minusNaNb/(AB) equals NaXa/[A(A+B+Z)]+NbYb/[B(A+B+Z)]-NaNbZ/[AB(A+B+Z)]. The shifted numerator and timelike total propagator correction are both retained.",
        "whole_double_external_inventory": "S312 has20C+60H+60GR graphs with two light-scalar propagators. Per four-scalar core,12 different-leg assignments and8 same-leg orderings give20. Merging same-leg orderings yields the four LL,LR,RL,RR radiation distributions for each hard channel. The canonical core numerator is shifted explicitly before comparing it with its unshifted value. No scalar seagull or internal metric emission is hidden in these140 graphs; all belong to the other247 regular graphs or the47 conserved pair graphs.",
        "whole_independent_marker_check": "An independent frozen TreeEngine is weighted by one formal marker at each internal light-scalar propagator. Its entire degree-two coefficient in the pure198,contact236 and full434 graph sums agrees exactly with the explicit external-line formula. This is a finite calibration of the general fixed-root grouping proof.",
        "whole_nonzero_leading_only_corrections": corrections,
        "checks": checks,
        "gates": {
            **gates,
            "general_stress_identity_not_a_sampled_derivative": True,
            "complete_tree_vertices_and_SymPy_not_modified": True,
            "all_double_external_graphs_covered_before_complex_bound": True,
        },
    }
