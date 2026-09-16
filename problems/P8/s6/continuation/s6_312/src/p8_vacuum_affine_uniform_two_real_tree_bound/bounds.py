"""Canonical vertex budgets and full434-tree original-parameter upper bound."""

from functools import cache

import sympy as s
from p8_vacuum_affine_complete_two_graviton_tree import trees as e
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree as matter
from p8_vacuum_affine_two_real_collinear_current import bounds as pair_bound

from . import gaps, source

DET = (s.S.One, s.Integer(4), s.Integer(48), s.Integer(960))
INV = (s.S.One, s.Integer(2), s.Integer(32), s.Integer(768))
DENSITY = (s.S.One, s.Integer(6), s.Integer(96), s.Integer(2400))
C3 = 4 * 6 * (512 * 3 * s.Rational(3, 2) ** 2 + 2 * 128 * 6 * s.Rational(3, 2))
C4 = (
    8
    * 24
    * (
        512 * 12 * s.Rational(3, 2) ** 2
        + 2 * 512 * 3 * 6 * s.Rational(3, 2)
        + 2 * 128 * 24 * s.Rational(3, 2)
        + 128 * 6**2
    )
)
R3 = 3 * C3 * 9600000
R4 = 3 * C4 * 9600000
MATTER_BUDGET = (33 * 4 + 177 * 8) * 1024**2 * s.Rational(8, 3) ** 2
GRAVITY_BUDGET = 177 * 1024**4 * s.Rational(8, 3) ** 2 * 900000 * max(1, R3, R3**2, R4)


def energies(first, second):
    first, second = map(e.exact_real, (first, second))
    if first <= 0 or second <= 0 or first + second > s.Rational(1, 8):
        raise ValueError("Require a,b>0 with a+b<=1/8")
    return first, second


def regular_upper(first, second):
    first, second = energies(first, second)
    return 1 / (s.Integer(10) ** 395 * first * second)


def full_upper(first, second):
    return 2 * regular_upper(first, second)


class WithoutPair(e.TreeEngine):
    """The canonical gauge-fixed complementary387 graphs, not an observable."""

    def current(self, mask, kind):
        if mask == 24 and kind == "h":
            return e.ZERO, 0
        return super().current(mask, kind)


@cache
def original_samples():
    checks = {}
    gates = {}
    record = []
    pars = {
        "heavy": source.HEAVY_MASS2,
        "cubic": source.CUBIC,
        "contact": source.CONTACT,
        "kappa": source.KAPPA,
    }
    for label, *args in gaps.SAMPLES[:4]:
        points, qs, born, pols = gaps.configuration(*args)
        Am = s.factor(
            matter.born_continuation(
                born, source.HEAVY_MASS2, source.CUBIC, source.CONTACT
            )
        )
        AG = s.factor(e.old.born(born) / source.KAPPA)
        A0 = Am + AG
        legs = [("phi", p, 1) for p in points[1:]] + [
            ("h", q, E) for q, E in zip(qs, pols)
        ]
        total, ntotal = e.TreeEngine(legs, **pars).amplitude()
        regular, nregular = WithoutPair(legs, **pars).amplitude()
        # Each chosen plus tensor has squared Frobenius norm2.
        # Multilinearity divides the two-polarization amplitude by2.
        total = s.factor(total / 2)
        regular = s.factor(regular / 2)
        paired = s.factor(total - regular)
        a, b = qs[0][0], qs[1][0]
        checks[label + "_whole_graph_count"] = s.Integer(ntotal - 434)
        checks[label + "_regular_graph_count"] = s.Integer(nregular - 387)
        checks[label + "_exact_split"] = s.factor(total - regular - paired)
        gates[label + "_positive_original_Born_parts"] = bool(Am > 0 and AG > 0)
        gates[label + "_regular_original_bound"] = bool(
            abs(regular) / A0 < regular_upper(a, b)
        )
        gates[label + "_retained_S311_pair_bound"] = bool(
            abs(paired) / A0 < pair_bound.selected_pair_upper(a, b)
        )
        gates[label + "_complete434_original_bound"] = bool(
            abs(total) / A0 < full_upper(a, b)
        )
        record.append(
            {
                "label": label,
                "scaled_full_square": s.factor((a * b * total / A0) ** 2),
                "scaled_regular_square": s.factor((a * b * regular / A0) ** 2),
            }
        )
    return checks, gates, record


@cache
def data():
    n, K = source.HEAVY_MASS2, source.KAPPA
    checks = {
        "density_first_coefficient_budget": DET[1] + INV[1] - DENSITY[1],
        "density_second_coefficient_budget": DET[2]
        + 2 * DET[1] * INV[1]
        + INV[2]
        - DENSITY[2],
        "density_third_coefficient_budget": DET[3]
        + 3 * DET[2] * INV[1]
        + 3 * DET[1] * INV[2]
        + INV[3]
        - DENSITY[3],
        "scalar_first_vertex_budget": 144 * DENSITY[1] + DET[1] - 868,
        "scalar_second_vertex_budget": 144 * DENSITY[2] + DET[2] - 13872,
        "scalar_third_vertex_budget": 144 * DENSITY[3] + DET[3] - 346560,
        "Einstein_cubic_component_budget": C3 - 138240,
        "Einstein_quartic_component_budget": C4 - 10616832,
        "hard_cubic_propagator_step": R3 - 3981312000000,
        "hard_quartic_propagator_step": R4 - 305764761600000,
        "regular_matter_component_budget": MATTER_BUDGET - 11542724608,
        "regular_Einstein_component_budget": GRAVITY_BUDGET
        - s.Integer("19742652106045453359826285363200000000000000000"),
        "first_trace_reversed_hard_propagator": 3 * gaps.GAP_DENOMINATOR - 900000,
        "soft_product_common_numerator": s.Rational(8, 3) ** 2 - s.Rational(64, 9),
    }
    a, b = s.symbols("a b", positive=True)
    z = s.Symbol("z", positive=True)
    checks["nested_soft_propagator_order_sum"] = s.cancel(
        1 / (a * (a + b)) + 1 / (b * (a + b)) - 1 / (a * b)
    )
    checks["squared_envelope_times_radial_measure"] = s.cancel(
        (2 / (s.Integer(10) ** 395 * a * b)) ** 2 * a * b
        - 4 / (s.Integer(10) ** 790 * a * b)
    )
    checks["soft_regulator_logarithm_derivative"] = s.diff(s.log(z / a), a) + 1 / a
    physical, physical_gates, _ = original_samples()
    checks.update(physical)
    return {
        "whole_vertex_component_budgets": {
            "determinant": DET,
            "inverse_metric": INV,
            "density_inverse": DENSITY,
            "Einstein_cubic": C3,
            "Einstein_quartic": C4,
            "hard_cubic_step": R3,
            "hard_quartic_step": R4,
            "regular_matter": MATTER_BUDGET,
            "regular_Einstein": GRAVITY_BUDGET,
        },
        "whole_regular_matter_bound": "The33 regular contact and177 regular H graphs have no internal graviton and total metric degree2. Scalar vertices/potential densities are bounded by1024^r, H metric vertices by n*1024^r; j<=2 heavy insertions with j+1 propagators give at most8/n. At most two light-scalar propagators contribute(8/3)^2/(ab). |C|<4g2/n gives |Mregular_m|<11542724608*g2/(n*kappa*ab). Dividing by Am>4g2/n3 gives an original coefficient below1e-395/(ab).",
        "whole_regular_Einstein_bound": "The177 GR graphs have light metric degree<=4, at most two scalar propagators, one initial hard graviton with bound900000/(delta+W2), followed by zero, one or two cubic steps or one quartic step. The paired L2/D bounds give total coefficient19742652106045453359826285363200000000000000000<1e48 over kappa2*(delta+W2)*ab. AG>8/(kappa*delta) implies an original relative coefficient below1e-395/(ab), uniformly as either hard transfer tends to zero.",
        "whole_complete_two_real_tree_theorem": "Positive Am and AG combine the regular sector bounds without adding their coefficients: |Mregular|/(Am+AG)<1e-395/(ab). The unchanged S311 pair bound adds9(a+b)/(sqrt(kappa)*ab), whose coefficient is below1e-395 on W<=1/8. Therefore the full434-tree amplitude obeys |M6|/(Am+AG)<2e-395/(ab), for all stated physical emitted directions and nonforward hard angles. Component bounds allow complex TT polarizations of unit spatial Frobenius norm. The bound is deliberately conservative, not a leading-soft subtraction or a small integrated detector-error estimate.",
        "whole_original_extreme_angle_calibrations": "Direct exact-rational full434 and complementary387 sums use original n,g,C,kappa at interior, near-collinear, near-hard-forward and simultaneous angular boundaries. Each tests the positive Born parts, regular bound, retained selected-pair bound and full bound. Unit polarization normalization uses exact multilinearity, not approximate square roots.",
        "whole_infrared_limit": "The bare envelope squared times radial phase space a da b db is4e-790 da db/(ab). Its integrals have logarithmic soft divergences; no nonzero prefactor cures them. Gauge-consistent overlap subtraction and real/virtual pairing remain required. The two gauge-fixed computational graph classes are not separately asserted as physical detector observables.",
        "checks": checks,
        "gates": {
            **physical_gates,
            "all_light_scalar_vertices_below1024_power": all(
                144 * DENSITY[r] + DET[r] < 1024**r for r in (1, 2, 3)
            ),
            "all_heavy_scalar_vertices_below_n1024_power": all(
                400 * DENSITY[r] + 128 * DET[r] < 128 * 1024**r for r in (1, 2)
            ),
            "all_potential_density_vertices_below1024_power": all(
                DET[r] <= 1024**r for r in (0, 1, 2)
            ),
            "matter_coefficient_below2e10": MATTER_BUDGET < 2 * 10**10,
            "Einstein_coefficient_below1e48": GRAVITY_BUDGET < 10**48,
            "original_matter_relative_coefficient_below1e_minus395": MATTER_BUDGET
            * n
            * n
            / (4 * K)
            < s.Rational(1, 10**395),
            "original_Einstein_relative_coefficient_below1e_minus395": GRAVITY_BUDGET
            / (8 * K)
            < s.Rational(1, 10**395),
            "selected_pair_coefficient_below1e_minus395": s.Rational(9, 8) / s.sqrt(K)
            < s.Rational(1, 10**395),
            "original_heavy_mass_at_least128": n >= 128,
            "positive_Born_weights_not_absolute_interference_replacement": True,
            "IR_divergent_envelope_not_reported_as_finite_rate": True,
            "full_physical_tree_not_individual_class_is_Ward_invariant": True,
        },
    }
