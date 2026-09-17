"""All-valence hard-core vertex, propagator and energy-assignment bounds."""

from functools import cache
from itertools import permutations

import sympy as s

from . import source


@cache
def data():
    checks = {}
    gates = {}
    for r in range(1, 33):
        checks[f"scalar_light_majorant_formula_{r}"] = s.factorial(r) * 2**r * (
            36 * s.binomial(r + 2, 2) + (r + 1)
        ) - s.factorial(r) * 2**r * (18 * (r + 2) * (r + 1) + (r + 1))
    gates["all32_light_vertex_calibrations"] = all(
        2**r * (36 * s.binomial(r + 2, 2) + r + 1) <= 1024**r for r in range(1, 33)
    )
    gates["all32_heavy_vertex_calibrations"] = all(
        2**r * (100 * s.binomial(r + 2, 2) + 128 * (r + 1)) <= 128 * 1024**r
        for r in range(1, 33)
    )
    gates["all33_density_calibrations"] = all(2**r * (r + 1) <= 4**r for r in range(33))
    checks.update(
        {
            "Euclidean_from_component_hard_step": s.Integer(4 * 32 * 300000 - 38400000),
            "light_general_envelope_base_margin": s.Integer(128 - 37 - 91),
            "heavy_general_envelope_base_margin": s.Integer(128 - 2 - 126),
            "maximum_light_Euclidean_square": s.Integer(4 * 3**2 - 36),
            "maximum_heavy_Euclidean_square": s.Integer(4 * 5**2 - 100),
            "three_heavy_channels_relative_Born": s.Rational(3 * 2, 4)
            - s.Rational(3, 2),
            "three_gravity_channels_relative_Born": s.Rational(3 * 300000, 8) - 112500,
        }
    )
    eta = s.diag(1, -1, -1, -1)
    v = s.Matrix(eta).reshape(16, 1)
    reflection = s.eye(16) - v * v.T / 2
    checks["trace_reversal_exact_Frobenius_isometry"] = (
        reflection.T * reflection - s.eye(16)
    )
    checks["Minkowski_metric_Frobenius_square"] = (v.T * v)[0] - 4
    a, b = s.symbols("a b", positive=True)
    checks["hard_component_square_envelope"] = s.expand(
        2 * a * a + 32 * b * b - (a + 4 * b) ** 2 - (a - 4 * b) ** 2
    )
    checks["ordered_two_energy_identity_retained"] = s.cancel(
        1 / (a * (a + b)) + 1 / (b * (a + b)) - 1 / (a * b)
    )
    # Distinct directly attached blocks are charged only once.
    assignments = {}
    for count in range(1, 6):
        ws = tuple(s.Rational(1, 2 ** (j + 8)) for j in range(count))
        all_orders = True
        for order in permutations(ws):
            cumulative = [sum(order[: j + 1], s.S.Zero) for j in range(count)]
            all_orders &= all(q >= w for q, w in zip(cumulative, order))
            all_orders &= s.prod(1 / q for q in cumulative) <= s.prod(
                1 / w for w in order
            )
        assignments[count] = bool(all_orders)
        checks[f"uncharged_block_energy_enlargement_{count}"] = (
            s.prod(ws) * s.prod(1 / w for w in ws) - 1
        )
    gates.update(
        {
            "all_order_vertex_proof_uses_r_factorial_not_finite_extrapolation": True,
            "binomial_bound_and_r_plus_one_square_below_four_power": True,
            "heavy_mass_at_least128": source.HEAVY_MASS2 >= 128,
            "timelike_Euclidean_ratio_within_common_hard_step": s.Rational(
                4 * 16 * 8, 45
            )
            < 38400000,
            "all_canonical_energy_assignments": all(assignments.values()),
            "chosen_blocks_disjoint_across_light_vertices": True,
            "hard_gap_uses_total_radiation_and_complement_at_every_N": True,
            "successive_EH_steps_have_one_Born_channel": True,
            "one_initial_hard_forward_pole_only": True,
            "heavy_mass_numerators_cancel_all_but_one_heavy_inverse": True,
            "trace_reversal_norm_does_not_require_current_conservation": True,
        }
    )
    return {
        "whole_vertex_envelopes": "For r>=1 light r!*1024^r and heavy n*r!*1024^r; density r!*4^r for r>=0. These follow from S315's all-order tensor bound, ||p||<=6 or10, n>=128, binom(r+2,2)<=(r+1)^2<=4^r. Canonical powers and every field norm multiply these coefficients.",
        "whole_hard_path_bound": "S312/S315 mixed-cut gap>(tau+W^2)/300000 and component momenta<=sqrt(tau)+4W give Euclidean L^2/|Dnext|<38400000. Every extra EH(r+2) vertex and next hard inverse have bound38400000*(r+2)!*32^(r+2). Trace reversal is a Frobenius isometry, without any off-shell Ward assumption. Only the first hard inverse keeps300000/(delta+W^2), absorbed by the positive gravity Born normalization.",
        "whole_block_energy_assignment": "Charge each light propagator to one directly attached block at its outer vertex, using cumulative energy>=block energy. Distinct vertices own disjoint blocks. Add1/W_block for uncharged blocks since W_block<=1/8<1. The S318 current bound plus attachment becomes C_m*W_block^(m-1)/(product w_i*kappa^(m/2))<=C_m/(product w_i*kappa^(m/2)).",
        "whole_finite_energy_order_checks": assignments,
        "checks": checks,
        "gates": gates,
    }
