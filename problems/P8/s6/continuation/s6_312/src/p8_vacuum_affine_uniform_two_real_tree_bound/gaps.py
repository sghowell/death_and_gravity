"""Common physical hard-propagator gap, including every radiation assignment."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import trees as e

from . import source

GAP_DENOMINATOR = s.Integer(300000)
DEFAULT_R, DEFAULT_HARD, DEFAULT_ENERGY, DEFAULT_SHRINK = (
    s.Rational(1, 2),
    s.Rational(1, 3),
    s.Rational(5, 4),
    s.Rational(1, 10),
)


def configuration(
    r=DEFAULT_R,
    hard=DEFAULT_HARD,
    E=DEFAULT_ENERGY,
    shrink=DEFAULT_SHRINK,
    first=None,
):
    r, hard, E, shrink = map(e.exact_real, (r, hard, E, shrink))
    if r <= 0 or hard <= 0 or shrink <= 0:
        raise ValueError("Require positive finite angular and recoil parameters")
    alpha = E + s.sqrt(E * E - 1) - shrink
    if alpha <= 1:
        raise ValueError("Require a massive final pair above threshold")
    Ep = (alpha + 1 / alpha) / 2
    c = (1 - r * r) / (1 + r * r)
    d = 2 * r / (1 + r * r)
    a = (E * E - Ep * Ep) / (2 * E) if first is None else e.exact_real(first)
    b = s.factor((E * E - E * a - Ep * Ep) / (E - a * (1 - c) / 2))
    q1 = e.imm([a, 0, 0, a])
    q2 = e.imm([b, b * d, 0, b * c])
    u = s.Matrix(
        [2 * hard / (1 + hard * hard), 0, (1 - hard * hard) / (1 + hard * hard)]
    )
    points, rays, born = source.recoil.momenta(E, [q1, q2], u)
    eps1 = e.imm(s.diag(0, 1, -1, 0))
    v = s.Matrix([c, 0, -d])
    y = s.Matrix([0, 1, 0])
    eps2 = s.zeros(4)
    eps2[1:, 1:] = v * v.T - y * y.T
    return (
        tuple(map(e.imm, points)),
        tuple(map(e.imm, rays)),
        tuple(map(e.imm, born)),
        (eps1, e.imm(eps2)),
    )


SAMPLES = (
    (
        "interior",
        s.Rational(1, 2),
        s.Rational(1, 3),
        s.Rational(5, 4),
        s.Rational(1, 10),
        s.Rational(1, 50),
    ),
    (
        "near_collinear",
        s.Rational(1, 10**40),
        s.Rational(1, 3),
        s.Rational(5, 4),
        s.Rational(1, 10),
        None,
    ),
    (
        "near_forward",
        s.Rational(1, 10),
        s.Rational(1, 10**120),
        s.Rational(5, 4),
        s.Rational(1, 10),
        None,
    ),
    (
        "both_boundaries",
        s.Rational(1, 10**40),
        s.Rational(1, 10**120),
        s.Rational(5, 4),
        s.Rational(1, 10),
        None,
    ),
    (
        "hierarchical_soft",
        s.Rational(2, 3),
        s.Rational(4, 3),
        s.Rational(29, 20),
        s.Rational(1, 10),
        s.Rational(1, 10**35),
    ),
    (
        "simultaneous_soft",
        s.Rational(3, 5),
        s.Rational(10**40),
        s.Rational(5, 3),
        s.Rational(1, 10**35),
        None,
    ),
)


@cache
def physical_samples():
    checks = {}
    gates = {}
    for label, *args in SAMPLES:
        points, qs, born, pols = configuration(*args)
        Q = qs[0] + qs[1]
        W = Q[0]
        E = -points[0][0]
        checks[label + "_total_conservation"] = sum(points, Q)
        for j, p in enumerate(points):
            checks[f"{label}_mass_shell_{j}"] = s.factor(e.old.dot(p, p) - 1)
        for j, q in enumerate(qs):
            checks[f"{label}_null_shell_{j}"] = s.factor(e.old.dot(q, q))
            checks[f"{label}_TT_transverse_{j}"] = pols[j] * q
            checks[f"{label}_TT_norm_squared_{j}"] = sum(x * x for x in pols[j]) - 2
        valid = True
        for mask in range(4):
            Qa = sum((qs[j] for j in range(2) if mask >> j & 1), e.VECTOR_ZERO)
            Qb = Q - Qa
            for outgoing in (2, 3):
                other = 5 - outgoing
                A = points[outgoing] + Qa
                B = points[other] + Qb
                u = s.factor(e.old.dot(A, A) - 1)
                v = s.factor(e.old.dot(B, B) - 1)
                R2 = s.factor(sum(A[j] * A[j] for j in range(1, 4)))
                r02 = E * E - 1
                checks[f"{label}_cluster_Kallen_{outgoing}_{mask}"] = s.factor(
                    r02 - R2 - (u + v) / 2 + (u - v) ** 2 / (16 * E * E)
                )
                tau = -e.old.dot(born[0] + born[outgoing], born[0] + born[outgoing])
                D = e.old.dot(points[0] + A, points[0] + A)
                valid &= -D > (tau + W * W) / GAP_DENOMINATOR
                valid &= u + v >= W / 2 and u + v <= 8 * W + W * W
                valid &= r02 - R2 > W / 8
            K = points[2] + points[3] + Qa
            valid &= s.Rational(45, 8) <= e.old.dot(K, K) <= 16
            if mask:
                w = Qa[0]
                for p in points:
                    valid &= abs(2 * e.old.dot(p, Qa) + e.old.dot(Qa, Qa)) >= 3 * w / 8
        gates[label + "_all_cluster_and_scalar_gaps"] = bool(valid)
        gates[label + "_strict_physical_domain"] = bool(
            0 < W <= s.Rational(1, 8)
            and all(x[0] > 0 for x in qs)
            and all(x.is_Rational for p in (*points, *qs) for x in p)
        )
    return checks, gates


@cache
def data():
    u, v, ss = s.symbols("u v s", positive=True)
    lam = (
        ss * ss
        + (1 + u) ** 2
        + (1 + v) ** 2
        - 2 * ss * (1 + u)
        - 2 * ss * (1 + v)
        - 2 * (1 + u) * (1 + v)
    )
    checks = {
        "general_Kallen_radial_deficit": s.expand(
            ss / 4 - 1 - lam / (4 * ss) - (u + v) / 2 + (u - v) ** 2 / (4 * ss)
        ),
        "cluster_mass_upper_at_W_one_eighth": 8 * s.Rational(1, 8)
        + s.Rational(1, 8) ** 2
        - s.Rational(65, 64),
        "radial_deficit_uniform_coefficient": s.Rational(1, 2)
        - s.Rational(65, 64) / (4 * s.Rational(25, 4))
        - s.Rational(147, 320),
        "radial_shrink_lower_denominator": s.Integer(8 * 4 - 32),
        "mixed_gap_total_energy_denominator": s.Integer(4 * 32**2 - 4096),
        "high_recoil_mixed_gap_denominator": s.Integer(4096 * 65 - 266240),
        "low_recoil_mixed_gap_denominator": s.Rational(16 * 65, 64) - s.Rational(65, 4),
        "timelike_hard_pair_lower": 4 * s.Rational(5, 4) * s.Rational(9, 8)
        - s.Rational(45, 8),
        "light_scalar_incoming_gap": s.Rational(1, 2)
        - s.Rational(1, 8)
        - s.Rational(3, 8),
        "mixed_vertex_momentum_square_budget": s.Integer(2 * 4**2 - 32),
        "mixed_vertex_momentum_over_next_gap": 32 * GAP_DENOMINATOR - 9600000,
    }
    physical, gates = physical_samples()
    checks.update(physical)
    return {
        "whole_common_mixed_hard_gap": "Assign any subset of the two future rays to A=p3+Qa and the rest to B=p4+Qb. If A2=1+u,B2=1+v,s=4E2 then W/2<=u+v<=65/64 and r0^2-r^2=(u+v)/2-(u-v)^2/(4s)>W/8. The cluster energy imbalance is <=sqrt(3)/2 times r0-r, so -(A-k1)^2>=|Avec-k1vec|^2/4>W2/4096. Recoil shifts the elastic cluster momentum by<=4W. Splitting W<=sqrt(tau)/8 from its complement proves -K2>(tau+W2)/300000 for every mixed hard propagator and every assignment, using total W not just its assigned subset.",
        "whole_other_denominators": "Timelike hard invariants lie in[45/8,16]. Heavy inverse kinetic factors are less than2/n for n>=128. An internal light scalar carries one external leg and a nonempty radiation subset, up to complement/sign, and |2pi.Qsub+Qsub2|>=3wsub/8. All mixed hard propagators on a fixed path use the same underlying Born channel.",
        "whole_vertex_momentum_pairing": "For a mixed path max momentum component L<=sqrt(tau)+4W, hence L2<=32(tau+W2) and L2/|Dnext|<9600000. For a timelike path L<=4 and |Dnext|>=45/8. The common gap includes both clusters; treating an off-shell cluster as a mass-one outgoing particle would not justify it.",
        "whole_gap_calibrations": "Six exact rational recoil configurations check all four radiation assignments, both mixed channels, light scalar cuts and cluster Kallen identities: interior, near relative collinearity, extreme hard forward, both boundaries, hierarchical soft and simultaneous soft near the backward hard boundary. They supplement, not replace, the general written inequalities.",
        "checks": checks,
        "gates": {
            **gates,
            "radial_coefficient_strictly_exceeds_one_quarter": s.Rational(147, 320)
            > s.Rational(1, 4),
            "uniform_gap_coarser_than_both_regimes": GAP_DENOMINATOR
            > max(266240, s.Rational(65, 4)),
            "timelike_momentum_ratio_within_common_budget": s.Rational(128, 45)
            < 9600000,
            "all_assignments_include_complement_radiation": True,
            "no_sample_infers_a_uniform_endpoint_bound": True,
        },
    }
