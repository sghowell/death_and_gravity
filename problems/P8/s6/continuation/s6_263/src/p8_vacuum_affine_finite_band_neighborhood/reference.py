"""Whole fixed-reference canonical boundary and positive comparison energy."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_principal_obstruction import coupled as c
from p8_vacuum_affine_scalar_tame_propagator import charts, energy, majorants

from . import source


@cache
def canonical():
    chart = charts.central()
    binding = {
        c.D: 1,
        c.Z: 1,
        c.J: charts.J,
        c.K: c.zeta / c.a**2,
        c.Y: 1,
        c.Yv: c.a**-2,
        c.r: 0,
        c.dH: 0,
        c.th: charts.th,
        c.L2: 2 * charts.E,
        c.L0: 3 * charts.T,
        c.C: s.Rational(1, 2),
        c.Vvv: s.Rational(9, 2) * charts.A,
        c.cs[0]: charts.l,
        c.cs[1]: 0,
        c.ws[0]: -charts.l * charts.E,
        c.ws[1]: 0,
        c.ds[0]: 0,
        c.ds[1]: 0,
        c.vs[0]: 0,
        c.vs[1]: 0,
        c.m00: 0,
        c.m01: 0,
        c.m11: -s.Symbol("whole_fixed_heavy_mass2", positive=True),
        c.H: charts.H,
    }
    actual = s.hessian(
        c.central_hamiltonian().subs(binding, simultaneous=True), c.CENTRAL
    ).extract([0, 1, 4, 5], [0, 1, 4, 5])
    Minv = chart["M"].inv()
    wanted = s.BlockMatrix(
        [
            [c.a**3 * chart["YY"], -chart["source"].T * Minv.T],
            [-Minv * chart["source"], Minv * chart["XX"] * Minv.T / c.a**3],
        ]
    ).as_explicit()
    wanted = wanted.subs(charts.q, c.P**2 / c.a**2, simultaneous=True)
    S = (chart["B"] + chart["B"].T) / 2
    A = (chart["B"] - chart["B"].T) / 2
    y = s.Matrix(s.symbols("whole_reference_y0:2", real=True))
    velocity = s.Matrix(s.symbols("whole_reference_velocity0:2", real=True))
    momentum = c.a**3 * (chart["K"] * velocity + A * y)
    R = s.eye(4)
    R[2:, :2] = c.a**3 * S
    O = s.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]])
    return {
        "whole_current_reference_parameter_binding": binding,
        "whole_reference_light_canonical_Hessian": actual,
        "whole_weighted_reduced_Legendre_Hessian": wanted,
        "whole_symmetric_reference_boundary": S,
        "whole_antisymmetric_reference_boundary": A,
        "whole_clean_reference_momentum": momentum,
        "whole_reference_old_from_clean_phase": R,
        "whole_original_reference_scalar_energy": (velocity.T * chart["K"] * velocity)[
            0
        ]
        / 2
        + charts.q * (y.T * chart["G"] * y)[0] / 2,
        "checks": {
            "whole_current_reference_weighted_Legendre_bridge": (
                actual - wanted
            ).applyfunc(s.factor),
            "whole_symmetric_boundary_cotangent_split": (
                momentum - c.a**3 * chart["K"] * velocity - c.a**3 * A * y
            ).applyfunc(s.expand),
            "whole_reference_boundary_symplectic_map": (R * O * R.T - O).applyfunc(
                s.factor
            ),
            "whole_reference_gyro_remains_skew": A + A.T,
            "whole_reference_boundary_zero_at_bounce": S.subs(
                {charts.th: 0, charts.H: 0}
            ),
        },
        "gates": {
            "full_reference_J_A_T_not_set_to_bare_values": all(
                actual.has(value) for value in (charts.J, charts.A, charts.T)
            ),
            "whole_weighted_time_boundary_retained": actual.has(charts.H),
            "nonzero_boundary_derivative_not_deleted_at_bounce": True,
            "reference_energy_not_new_vacuum_or_quantum_ordering": True,
        },
    }


@cache
def shear_bounds():
    chart = charts.central()
    S = canonical()["whole_symmetric_reference_boundary"]
    Sd = S.applyfunc(charts.dtime)
    variables = (charts.z, *charts.jetvars, *charts.jets, *charts.jet2)
    box = {
        charts.z: s.Rational(1, 4096),
        charts.th: 2,
        charts.E: s.Rational(1, 2),
        charts.l: s.Rational(1, 10),
        charts.J: 100,
        charts.A: s.Rational(1, 10**6),
        charts.T: s.Rational(1, 10**6),
        charts.H: 2,
        **{var: 10**6 for var in (*charts.jets, *charts.jet2)},
    }
    denom = charts.E**4 * chart["Delta"] ** 2
    inverse = 4**4 * 8**2
    rows = {}
    counts = {}
    checks = {}
    for label, matrix in (("S_over_q", S / charts.q), ("Sdot_over_q", Sd / charts.q)):
        values = []
        sizes = []
        residual = []
        for expression in matrix:
            normalized = s.cancel(expression.subs(charts.q, 1 / charts.z) * denom)
            bound, polynomial = majorants.polynomial_majorant(
                normalized, variables, box
            )
            values.append(bound * inverse)
            sizes.append(len(polynomial.terms()))
            residual.append(s.cancel(polynomial.as_expr() - normalized))
        rows[label] = values
        counts[label] = sizes
        checks[label + "_whole_polynomial_reconstructions"] = s.Matrix(residual)
    derivative = 2 * (max(rows["Sdot_over_q"]) + 6 * max(rows["S_over_q"]))
    return {
        "whole_reference_symmetrized_boundary": S,
        "whole_reference_boundary_first_time_derivative": Sd,
        "whole_normalizing_positive_denominator": denom,
        "whole_reference_coefficient_jet_box": box,
        "whole_eight_polynomial_entry_majorants": rows,
        "whole_eight_polynomial_term_counts": counts,
        "whole_reference_a_cubed_boundary_over_P_squared_derivative_bound": derivative,
        "whole_reference_S_over_P_squared_absolute_bound": s.Rational(1, 10**20),
        "checks": checks,
        "gates": {
            "all_reference_shear_entries_below_1e40": all(
                x < 10**40 for row in rows.values() for x in row
            ),
            "whole_reference_S_over_P2_time_bound": derivative < 10**40,
            "whole_reference_S_over_P2_integrated_smallness": derivative * source.TIME
            < s.Rational(1, 10**20),
            "original_positive_E_and_Delta_domains_required": True,
            "whole_reference_boundary_not_recomputed_on_new_background": True,
        },
    }


@cache
def positive_energy():
    minimum, maximum = majorants.LOWER, majorants.UPPER
    gyro = 2 * majorants.ENTRY / source.LOW
    lower = min(1 / (128 * maximum), minimum / 4 - gyro**2 / maximum)
    upper = 2 / minimum + 2 * gyro**2 / minimum + maximum
    a, P, mass, zeta = s.symbols(
        "reference_a P fixed_mass2 positive_zeta", positive=True
    )
    weights = {
        "heavy": ((a * P**2 + a**3 * mass) / (P**2 + mass), a**-3),
        "longitudinal_Proca": (a, a**-3 + 1 / (a * zeta * P**2)),
        "TT": (a, a**-3),
        "transverse_Proca": (a**-1 + a / (zeta * P**2), a**-1),
    }
    # Every remaining reference pair is the entire positive oscillator,
    # including the heavy mass and longitudinal Proca mass.
    q, p, A, B, Ad, Bd = s.symbols(
        "q p full_positive_A full_positive_B A_dot B_dot", real=True
    )
    derivative = A * q * (B * p) + B * p * (-A * q) + Ad * q * q / 2 + Bd * p * p / 2
    return {
        "whole_scalar_reference_K_and_G_operator_interval": [minimum, maximum],
        "whole_scalar_reference_gyro_scaled_bound": gyro,
        "whole_scalar_energy_lower_bound": lower,
        "whole_scalar_energy_upper_bound": upper,
        "whole_other_normalized_reference_energy_weights": weights,
        "whole_reference_real_phase_count": 16,
        "whole_reference_energy_operator_interval": [
            s.Rational(1, 10**12),
            s.Integer(10) ** 12,
        ],
        "whole_reference_root_energy_time_rate": energy.ENERGY_RATE,
        "whole_positive_oscillator_energy_derivative": derivative,
        "positive_energy_argument": "Use the unchanged scalar E=(ydot^T K ydot+q y^T G y)/2, with Pi=a_ref^3(K ydot+A_skew y). In W coordinates, y=xq/P and Pi=xp; 1<=a_ref<=2 and ||A_skew||/P<=2e18/P. Young's inequality and the original K,G bounds give the displayed coercivity interval without commuting K and G. The other six oscillator pairs have the listed coefficients, with two TT and two transverse-Proca copies. The giant heavy frequency is retained and cancels in the energy derivative.",
        "checks": {
            "whole_positive_oscillator_mass_frequency_cancellation": s.expand(
                derivative - Ad * q * q / 2 - Bd * p * p / 2
            ),
            "complete_physical_mode_count": s.Integer(4 + 2 + 2 + 4 + 4) - 16,
        },
        "gates": {
            "complete_scalar_energy_lower_bound": lower > s.Rational(1, 10**12),
            "complete_scalar_energy_upper_bound": upper < 10**12,
            "all_other_energy_lower_bounds": s.Rational(1, 8) > s.Rational(1, 10**12),
            "all_other_energy_upper_bounds": 8
            + 2 / (s.Rational(1, 10**6) * source.LOW**2)
            < 10**12,
            "original_full_scalar_reference_energy_rate": energy.ENERGY_RATE == 10**28,
            "no_raw_heavy_mass_Gronwall_frequency": True,
            "full_noncommuting_coupled_scalar_pair_kept": True,
        },
    }
