"""Original current budget and exact47-graph factorization calibrations."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import checks as samples
from p8_vacuum_affine_complete_two_graviton_tree import trees as e
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree as matter

from . import source


def require_energy(value):
    value = e.exact_real(value)
    if not 0 < value <= s.Rational(1, 8):
        raise ValueError("Require a positive energy at most1/8")
    return value


def hard_current_relative_upper(total):
    total = require_energy(total)
    return 1 / (64 * total)


def selected_pair_upper(first, second):
    first, second = require_energy(first), require_energy(second)
    require_energy(first + second)
    return 9 * (first + second) / (s.sqrt(source.KAPPA) * first * second)


@cache
def factorization():
    points, qs, eps = samples.nonopposite_state()
    Q = qs[0] + qs[1]

    class WithoutPair(e.TreeEngine):
        def current(self, mask, kind):
            if mask == 24 and kind == "h":
                return e.ZERO, 0
            return super().current(mask, kind)

    checks = {}
    for label, pars, expected, paircount in (
        ("pure", {"heavy": 128, "cubic": 0, "contact": 0, "kappa": 1}, 198, 21),
        (
            "full",
            {
                "heavy": 128,
                "cubic": s.Rational(2, 3),
                "contact": s.Rational(-1, 7),
                "kappa": 1,
            },
            434,
            47,
        ),
        (
            "original",
            {
                "heavy": source.HEAVY_MASS2,
                "cubic": source.CUBIC,
                "contact": source.CONTACT,
                "kappa": source.KAPPA,
            },
            434,
            47,
        ),
    ):
        legs = [("phi", p, 1) for p in points[1:]] + [
            ("h", q, E) for q, E in zip(qs, eps)
        ]
        whole, count = e.TreeEngine(legs, **pars).amplitude()
        regular, nregular = WithoutPair(legs, **pars).amplitude()
        hard = e.TreeEngine([("phi", p, 1) for p in points], **pars)
        U, nhard = hard.amputated(hard.full, "h")
        H = e.imm(e.ETA * U * e.ETA - e.ETA * s.trace(e.ETA * U) / 2)
        paired = -e.cubic_gravity((eps[0], eps[1], H), (qs[0], qs[1], -Q)) / (
            s.sqrt(pars["kappa"]) * e.old.dot(Q, Q)
        )
        checks[label + "_whole_minus_regular_minus_pair"] = s.factor(
            whole - regular - paired
        )
        checks[label + "_whole_count"] = s.Integer(count - expected)
        checks[label + "_complementary_count"] = s.Integer(
            nregular - expected + paircount
        )
        checks[label + "_hard_current_count"] = s.Integer(nhard - paircount)
        checks[label + "_hard_current_divergence"] = (Q.T * e.ETA * U).applyfunc(
            s.factor
        )
    return checks


@cache
def original_samples():
    E = s.Rational(5, 4)
    alpha = s.Rational(19, 10)
    Ep = (alpha + 1 / alpha) / 2
    a = s.Rational(1, 50)
    pars = {
        "heavy": source.HEAVY_MASS2,
        "cubic": source.CUBIC,
        "contact": source.CONTACT,
        "kappa": source.KAPPA,
    }
    checks = {}
    gates = {}
    for label, r, hard in (
        ("interior", s.Rational(1, 2), None),
        ("near_collinear", s.Rational(1, 10**40), None),
        ("near_hard_forward", s.Rational(1, 10), s.Rational(1, 10**120)),
        ("both_boundaries", s.Rational(1, 10**40), s.Rational(1, 10**120)),
    ):
        c = (1 - r * r) / (1 + r * r)
        d = 2 * r / (1 + r * r)
        b = s.factor((E * E - E * a - Ep * Ep) / (E - a * (1 - c) / 2))
        q1 = e.imm([a, 0, 0, a])
        q2 = e.imm([b, b * d, 0, b * c])
        Q = q1 + q2
        W = a + b
        u = (
            s.Matrix([s.Rational(3, 5), s.Rational(4, 5), 0])
            if hard is None
            else s.Matrix(
                [2 * hard / (1 + hard * hard), 0, (1 - hard * hard) / (1 + hard * hard)]
            )
        )
        points, _, born = source.recoil.momenta(E, [q1, q2], u)
        A0 = s.factor(
            matter.born_continuation(
                born, source.HEAVY_MASS2, source.CUBIC, source.CONTACT
            )
            + e.old.born(born) / source.KAPPA
        )
        engine = e.TreeEngine([("phi", p, 1) for p in points], **pars)
        U, count = engine.amputated(engine.full, "h")
        normsq = sum(U[i, j] ** 2 for i in range(1, 4) for j in range(1, 4))
        eps1 = e.imm(s.diag(0, 1, -1, 0))
        v = s.Matrix([c, 0, -d])
        y = s.Matrix([0, 1, 0])
        eps2 = s.zeros(4)
        eps2[1:, 1:] = v * v.T - y * y.T
        eps2 = e.imm(eps2)
        H = e.imm(e.ETA * U * e.ETA - e.ETA * s.trace(e.ETA * U) / 2)
        pair = -e.cubic_gravity((eps1, eps2, H), (q1, q2, -Q)) / (
            2 * s.sqrt(source.KAPPA) * e.old.dot(Q, Q)
        )
        checks[label + "_current_count"] = s.Integer(count - 47)
        checks[label + "_exact_current_divergence"] = (Q.T * e.ETA * U).applyfunc(
            s.factor
        )
        checks[label + "_total_conservation"] = sum(points, Q)
        gates[label + "_exact_rational_recoil"] = all(
            all(x.is_Rational for x in p) for p in points
        )
        gates[label + "_original_current_bound"] = (
            normsq / A0**2 < hard_current_relative_upper(W) ** 2
        )
        gates[label + "_original_pair_bound"] = abs(pair) / A0 < selected_pair_upper(
            a, b
        )
    return checks, gates


@cache
def data():
    n, K = source.HEAVY_MASS2, source.KAPPA
    matter_external = 4 * 72 * 4 + 3 * 4 * 72 * 2
    matter_regular = 4 + 3 * 8 + 3 * 2 * 2
    matter_constant = matter_external + s.Rational(1, 8) * matter_regular
    matter_relative = 3072 * n * n / s.sqrt(K)
    seagull = 4 * (17 * 12 + 128 + 128 + 256 + 256)
    cubic = 4 * 6 * (512 * 3 * s.Rational(3, 2) ** 2 + 2 * 128 * 6 * s.Rational(3, 2))
    external_budget = 12 * 72 * 16 * 27 * 270000
    seagull_budget = 6 * seagull * 270000
    mixed_cubic_budget = 2 * cubic * 18 * 270000**2
    timelike_cubic_budget = cubic * 16 * 2**2 * s.Rational(65, 64)
    gravity_budget = (
        external_budget + seagull_budget + mixed_cubic_budget + timelike_cubic_budget
    )
    gravity_relative = s.Rational(3, 8) * 10**18 / s.sqrt(K)
    checks = {
        "matter_external_component_count": s.Integer(matter_external - 2880),
        "matter_regular_component_count": s.Integer(matter_regular - 40),
        "matter_component_budget": matter_constant - 2885,
        "general_root_seagull_bound": s.Integer(seagull - 3888),
        "general_root_Einstein_cubic_bound": cubic - 138240,
        "gravity_external_budget": s.Integer(external_budget - 100776960000),
        "gravity_seagull_budget": s.Integer(seagull_budget - 6298560000),
        "gravity_mixed_cubic_budget": mixed_cubic_budget - 362797056000000000,
        "gravity_timelike_cubic_budget": timelike_cubic_budget - 8985600,
        "gravity_total_component_budget": gravity_budget - 362797163084505600,
        "source_heavy_hierarchy": n - (s.Rational(10**200, 512) + 2),
        "recoil_low_high_gap_worst_denominator": s.Integer(768 * 37 - 28416),
        "outgoing_timelike_pair_lower": 4 * s.Rational(5, 4) * s.Rational(9, 8)
        - s.Rational(45, 8),
        "external_incoming_scalar_gap": s.Rational(1, 2)
        - s.Rational(1, 8)
        - s.Rational(3, 8),
        "selected_pair_coefficient_margin": s.Integer(9)
        - s.Rational(530, 64)
        - s.Rational(23, 32),
    }
    checks.update(factorization())
    physical, physical_gates = original_samples()
    checks.update(physical)
    return {
        "whole_exact_current_budgets": {
            "matter_component_coefficient": matter_constant,
            "matter_relative_Frobenius_coefficient": matter_relative,
            "gravity_component_coefficient": gravity_budget,
            "gravity_relative_Frobenius_coefficient": gravity_relative,
        },
        "whole_uniform_current_theorem": "For the original S300 recoil W<=1/8 and E in[5/4,2], the complete47-tree current obeys||Uspatial||F/(Am+AG)<1/(64W), at every nonforward elastic angle. The matter graph budget is below4096g2/(n sqrtK W) per component. The gravity budget is below1e18/[K^(3/2) W(delta+W2)] per component, using the physical mixed-transfer gap and mixed-cubic momentum suppression. Positive Born fractions combine both norm estimates without omitting any graph.",
        "whole_selected_pair_theorem": "For unit-Frobenius physical TT polarizations and a,b>0,a+b<=1/8, the exact47 pair-propagator graphs satisfy|Mpair|/(Am+AG)<9(a+b)/(sqrtK ab)=9e-400(1/a+1/b), uniformly in relative emitted angle and all nonforward hard angles. The full434 amplitude equals this class plus the complementary387 graphs, whose bound is not asserted here.",
        "whole_original_boundary_calibrations": "Exact original-parameter rational states check the complete current norm and unit-TT pair bound in the interior, at half-angle1e-40, at a hard-angle half-parameter1e-120, and at both boundaries together. Squared norms are compared exactly; no approximate square root or altered SymPy backend enters a scientific comparison.",
        "whole_integration_scope": "The retained1/a+1/b behavior is not an integrated IR-finite rate. Soft overlaps and virtual pairing must be subtracted consistently before integrating the full amplitude. Angular pole cancellation of the selected conserved class does not prove an all-N detector theorem, matching condition or full P8 closure.",
        "checks": checks,
        "gates": {
            **physical_gates,
            "original_heavy_mass_above128": n >= 128,
            "matter_component_bound_below4096": matter_constant < 4096,
            "gravity_component_bound_below1e18": gravity_budget < 10**18,
            "matter_relative_current_coefficient_below_one_over64": matter_relative
            < s.Rational(1, 64),
            "gravity_relative_current_coefficient_below_one_over64": gravity_relative
            < s.Rational(1, 64),
            "selected_pair_coefficient_below9": s.Rational(530, 64) < 9,
            "original_positive_Born_fractions_not_subtracted": True,
            "uniform_angle_bound_not_an_integrated_rate": True,
        },
    }
