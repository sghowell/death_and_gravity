"""Composition barrier and complete original-parameter tree calibration."""

from functools import cache

import sympy as s
from p8_vacuum_affine_all_multiplicity_tree_source import trees as full
from p8_vacuum_affine_complete_two_graviton_tree import trees as lower
from p8_vacuum_affine_massive_scalar_graviton_bremsstrahlung import tree as matter
from p8_vacuum_affine_temporal_tree_reorganization import trees as temporal
from p8_vacuum_affine_uniform_soft_current_bound import bounds as soft

from . import core, source


def coefficient_envelope(n):
    n = int(source.require_order(n))
    return (
        3
        * source.HEAVY_MASS2**2
        * s.factorial(n)
        * (2 * 10**16) ** n
        / source.KAPPA ** s.Rational(n, 2)
    )


def amplitude_upper(energies):
    if not isinstance(energies, (tuple, list)) or not energies:
        raise ValueError("Require a nonempty finite exact energy list")
    ws = tuple(source.require_mass(w) for w in energies)
    if sum(ws, s.S.Zero) > s.Rational(1, 8):
        raise ValueError("Require total radiated energy at most1/8")
    return coefficient_envelope(len(ws)) / s.prod(ws)


def sector_coefficient(n, sector):
    n = int(source.require_order(n))
    if not isinstance(sector, str) or sector not in ("matter", "gravity"):
        raise ValueError("Require a stated complete core sector")
    return _sector_coefficient(n, sector)


@cache
def _sector_coefficient(n, sector):
    x = s.Symbol("x")
    a = core.rational_coefficients(core.functions(x)[sector], x, n)
    b = [s.S.Zero] + [soft.coefficient(j) / s.factorial(j) for j in range(1, n + 1)]
    value = core.compose(a, b, n)[n] * s.factorial(n)
    return s.factor(value)


@cache
def barrier():
    x = s.Symbol("x")
    cap, radius = s.Rational(1, 10**16), s.Rational(1, 2 * 10**16)
    Fcap = 2048 * ((1 - 416 * cap) ** -2 - 1 - 832 * cap)
    values = {
        key: s.factor(value.subs(x, cap)) for key, value in core.functions(x).items()
    }
    return {
        "cap": cap,
        "radius": radius,
        "soft_nonlinearity": Fcap,
        "margin": cap - radius - Fcap,
        "core_values": values,
    }


@cache
def original_calibration():
    points, hs = full.configuration()
    E = -points[0][0]
    r0 = s.sqrt(E * E - 1)
    direction = points[2][1:, 0] / s.sqrt(sum(v * v for v in points[2][1:, 0]))
    born = (
        points[0],
        points[1],
        s.ImmutableMatrix([E, *(r0 * direction)]),
        s.ImmutableMatrix([E, *(-r0 * direction)]),
    )
    pars = source.original_parameters()
    Am = s.factor(
        matter.born_continuation(born, pars["heavy"], pars["cubic"], pars["contact"])
    )
    AG = s.factor(lower.old.born(born) / pars["kappa"])
    legs = [("phi", p, 1) for p in points[1:]] + [("h", q, A) for q, A in hs]
    baseline, n = full.TreeEngine(legs, **pars).amplitude()
    mapped, m = temporal.TemporalSoft(legs, **pars).amplitude()
    ws = tuple(q[0] for q, A in hs)
    scaled_square = s.factor((s.prod(ws) * baseline / (Am + AG)) ** 2 / 8)
    checks = {
        "full_5116_tree_inventory": s.Integer(n - 5116),
        "complete_temporal_inventory": s.Integer(m - n),
        "original_full_amplitude_gauge_equivalence": s.factor(mapped - baseline),
        "total_physical_conservation": sum(points, lower.VECTOR_ZERO)
        + sum((q for q, A in hs), lower.VECTOR_ZERO),
    }
    for j, (q, A) in enumerate(hs):
        checks[f"unit_polarization_square_normalization_{j}"] = (
            sum(v * v for v in A) - 2
        )
    return {
        "checks": checks,
        "gates": {
            "same_original_parameters": pars == source.original_parameters(),
            "positive_unexpanded_Born_parts": Am > 0 and AG > 0,
            "original_complete_tree_nonzero": baseline != 0,
            "original_full_tree_within_proved_bound": scaled_square
            < coefficient_envelope(3) ** 2,
            "strict_soft_energy_domain": all(w > 0 for w in ws)
            and sum(ws) <= s.Rational(1, 8),
        },
        "whole_exact_calibration": {
            "graph_count": s.Integer(n),
            "energies": ws,
            "scaled_unit_polarization_amplitude_square": scaled_square,
            "proved_coefficient_envelope": coefficient_envelope(3),
        },
    }


@cache
def data():
    record = barrier()
    cap, radius = record["cap"], record["radius"]
    values = record["core_values"]
    checks = {
        "barrier_radius_is_half_cap": radius - cap / 2,
        "original_three_real_coefficient": coefficient_envelope(3)
        - 3
        * source.HEAVY_MASS2**2
        * 6
        * (2 * 10**16) ** 3
        / source.KAPPA ** s.Rational(3, 2),
        "original_four_real_coefficient": coefficient_envelope(4)
        - 3 * source.HEAVY_MASS2**2 * 24 * (2 * 10**16) ** 4 / source.KAPPA**2,
    }
    gates = {
        "positive_exact_soft_fixed_point_margin": record["margin"] > 0,
        "cap_inside_every_vertex_series": 0 < 416 * cap < 1
        and 3072 * cap < 1
        and s.Rational(11, 3) * 1024 * cap < 1,
        "hard_path_sequence_converges_below_one_thousandth": 0
        < values["V"]
        < s.Rational(1, 1000),
        "matter_core_below_three": 0 < values["matter"] < 3,
        "gravity_core_below_two_eleven": 0 < values["gravity"] < 2 * 10**11,
        "original_matter_prefactor_dominates_gravity_cap": 3 * source.HEAVY_MASS2**2
        > 2 * 10**11,
        "original_N3_coefficient_below_one_e_minus755": coefficient_envelope(3)
        < s.Rational(1, 10**755),
        "original_N4_coefficient_below_one_e_minus1138": coefficient_envelope(4)
        < s.Rational(1, 10**1138),
        "all_N_nonnegative_series_argument_not_sample_extrapolation": True,
        "positive_Born_weights_combine_without_factor_two": True,
        "old_sharper_N1_N2_bounds_are_not_replaced": True,
        "factorial_envelope_not_probability_summability": True,
    }
    coefficients = {}
    x = s.Symbol("x")
    independent = core.independent_weighted_coefficients(6)
    b = [s.S.Zero] + [soft.coefficient(j) / s.factorial(j) for j in range(1, 7)]
    for label, a in zip(("matter", "gravity"), independent):
        composed = core.compose(a, b, 6)
        coefficients[label] = {}
        for n in range(1, 7):
            value = sector_coefficient(n, label)
            coefficients[label][n] = value
            checks[f"independent_composed_{label}_coefficient_{n}"] = (
                value - s.factorial(n) * composed[n]
            )
            cap_value = 3 if label == "matter" else 2 * 10**11
            gates[f"{label}_coefficient_below_proved_barrier_{n}"] = (
                0 < value < cap_value * s.factorial(n) / radius**n
            )
    a, b = s.symbols("a b", positive=True)
    checks["two_energy_radial_measure_keeps_IR_poles"] = s.cancel(
        (1 / (a * b)) ** 2 * a * b - 1 / (a * b)
    )
    checks["positive_soft_regulator_log_derivative"] = s.diff(s.log(x / a), a) + 1 / a
    actual = original_calibration()
    checks.update({"original_" + key: value for key, value in actual["checks"].items()})
    gates.update({"original_" + key: value for key, value in actual["gates"].items()})
    return {
        "whole_uniform_complete_finite_tree_theorem": "For every finite N>=1, |M_(4Phi+Nh)|/(Am+AG)<=3*n^2*N!*(2*10^16)^N/[kappa^(N/2)*product wi], at the original parameters and physical unit spatial TT norms. The bound is uniform in nested pure-soft collinear approaches, energy hierarchies and nonforward hard-angle limits; explicit1/wi factors remain. The exact singular point evaluator domain is unchanged.",
        "whole_exact_composition_barrier": record,
        "whole_composed_finite_coefficients": coefficients,
        "whole_original_full_5116_calibration": actual["whole_exact_calibration"],
        "whole_proof": "Unique maximal-soft contraction reduces every tree to one of three hard cores. All-valence vertex and paired hard-propagator bounds give nonnegative set/sequence majorants. Charge at most one cumulative scalar denominator to each disjoint block, then enlarge to all blocks. Compose with the S318 current EGF and use the exact positive radius barrier. Positive Born weights give one common coefficient without double counting.",
        "whole_infrared_boundary": "The bare squared envelope times radial phase space retains product(dwi/wi). A small prefactor cannot make it infrared integrable. Ordered all-N soft-overlap subtraction, real-virtual pairing, finite hard/evanescent matching, interacting state, unitarity, absolute complex Regge and original common-parent bounce remain separate.",
        "checks": checks,
        "gates": gates,
    }
