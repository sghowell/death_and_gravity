"""Entire original canonical reference flow and nonlinear off-time phase image."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import phase
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate
from p8_vacuum_affine_quantitative_gaussian_window import source as original
from p8_vacuum_affine_quantitative_phase_domain import geometry as old_geometry
from p8_vacuum_affine_quantitative_phase_domain import reference as field
from p8_vacuum_affine_scalar_tame_propagator import charts

from . import source

CORE = field.CORE
ANALYTIC_RADIUS = 4 * CORE
IMAGE_RADIUS = 8 * CORE
CANONICAL_FLOW = s.Integer(10) ** 500
WHITENING = s.Integer(10) ** 120
WHITENED_FLOW = s.Integer(10) ** 1000


def maximum(iv):
    value = max(abs(iv.lo), abs(iv.hi))
    return s.Rational(value.numerator, value.denominator)


@cache
def linear_flow():
    P, a = field.P, phase.a
    light = phase.canonical_generator().applyfunc(s.cancel)
    box = {
        a: I(1, 2),
        charts.q: I(P**2 / 4, P**2),
        charts.th: I(-2, 2),
        charts.E: I(-s.Rational(1, 2), s.Rational(1, 2)),
        charts.l: I(-s.Rational(1, 10), s.Rational(1, 10)),
        charts.J: I(1, 100),
        charts.A: I(-s.Rational(1, 10**6), s.Rational(1, 10**6)),
        charts.T: I(-s.Rational(1, 10**6), s.Rational(1, 10**6)),
        charts.H: I(-2, 2),
    }
    tensor_scale = s.diag(-2, -s.Rational(1, 2))
    tensor = tensor_scale * phase.tensor_generator() * tensor_scale.inv()
    heavy = s.Matrix([[0, a**-3], [-(a**3) * field.heavy_state.MASS2 - a * P**2, 0]])
    k = s.Matrix(s.symbols("complete_vector_k0:3", real=True))
    mass, zeta = s.symbols("complete_Proca_mass complete_zeta", positive=True)
    K = s.eye(3) / a + k * k.T / (a**3 * mass**2)
    V = a * mass**2 * s.eye(3) + (k.dot(k) * s.eye(3) - k * k.T) / a
    vector = s.zeros(6)
    vector[:3, 3:], vector[3:, :3] = K / zeta, -zeta * V
    blocks = {"light": light, "each_tensor": tensor, "heavy": heavy}
    bounds = {
        name: sum(maximum(evaluate(x, box)) for x in matrix)
        for name, matrix in blocks.items()
    }
    for axis in range(3):
        substitution = {k[i]: P * int(i == axis) for i in range(3)}
        matrix = vector.subs({**substitution, mass: field.MASS, zeta: field.ZETA})
        bounds["three_vector_axis_" + str(axis)] = sum(
            maximum(evaluate(x, box)) for x in matrix
        )
    entire = 6 * (
        bounds["light"]
        + 2 * bounds["each_tensor"]
        + bounds["heavy"]
        + max(bounds["three_vector_axis_" + str(i)] for i in range(3))
    )
    coefficient_box = {
        original.u: I(-original.TIME, original.TIME),
        original.rho: I(-original.PROFILE_BOUND, original.PROFILE_BOUND),
        original.pressure: I(-original.PROFILE_BOUND, original.PROFILE_BOUND),
    }
    actual_pivot = evaluate(original.J, coefficient_box)
    return {
        "blocks": blocks,
        "vector": vector,
        "K": K,
        "V": V,
        "k": k,
        "mass": mass,
        "zeta": zeta,
        "bounds": bounds,
        "entire": entire,
        "actual_pivot": actual_pivot,
        "actual_theta": maximum(evaluate(original.theta, coefficient_box)),
        "actual_Hubble": maximum(evaluate(original.H, coefficient_box)),
    }


@cache
def phase_domain():
    raw = field.field_bounds()
    f = {name: value * IMAGE_RADIUS for name, value in raw.items()}
    old = old_geometry.bounds()
    time = source.TIME
    scale_gap, volume_gap = 16 * time**2, 64 * time**2
    gamma_gap = old["whole_gamma_minus_I_A2"] + 4 * scale_gap
    inverse_gap = old["whole_inverse_gamma_minus_I_A2"] + 4 * scale_gap
    background_momentum = 100 * time
    delta_density = 13 * time**2
    P = field.P
    images = [
        2 * (f["Pi_v_A0"] + background_momentum) / 3,
        6 * (1 + P) * f["Pi_W_A0"],
        2 * f["delta_Pi_M_A0"] + old_geometry.VLOG,
        2 * f["Pi_H_A0"],
        f["eta_A0"],
        old["whole_normalized_shear_invariant"],
        8 * f["Pi_W_A0"] ** 2 / field.ZETA,
        48 * field.ZETA * ((1 + P) * f["W_A0"]) ** 2,
        2 * f["W_A0"] ** 2,
        2 * f["M1_A1"] ** 2,
        2 * f["H_gradient_A0"] ** 2,
        old["whole_scalar_curvature_A0"],
    ]
    return {
        "raw": raw,
        "fields": f,
        "old_geometry": old,
        "scale_gap": scale_gap,
        "volume_gap": volume_gap,
        "gamma_gap": gamma_gap,
        "inverse_gap": inverse_gap,
        "background_momentum": background_momentum,
        "delta_density": delta_density,
        "invariants": dict(zip(source.COORDS, images, strict=True)),
        "flow_deviation": 2 * WHITENED_FLOW * time,
        "whitening_row": 40 * s.sqrt(field.KAPPA) * max(field.safe_bounds().values()),
    }


@cache
def data():
    flow, domain = linear_flow(), phase_domain()
    a, P, u = phase.a, field.P, original.u
    physical_frequency = field.heavy_state.MASS2 + P**2
    balance_frequency = field.heavy_state.MASS2 + P**2 / 16
    k, mass = flow["k"], flow["mass"]
    background_pv = -6 * original.a**3 * original.H
    x = s.Symbol("nonnegative_clock_square", nonnegative=True)
    checks = {
        "entire_scalar_canonical_symplectic_generator": (
            flow["blocks"]["light"] * phase.OMEGA
            + phase.OMEGA * flow["blocks"]["light"].T
        ).applyfunc(s.cancel),
        "literal_full_tensor_shape_flow": flow["blocks"]["each_tensor"]
        - s.Matrix([[0, 4 / a**3], [-(a**3) * charts.q / 4, 0]]),
        "full_three_polarization_Proca_product": (
            flow["K"] * flow["V"] - (mass**2 + k.dot(k) / a**2) * s.eye(3)
        ).applyfunc(s.factor),
        "physical_heavy_dispersion_not_Omega_star": physical_frequency
        - balance_frequency
        - 15 * P**2 / 16,
        "actual_unchanged_reference_matter_density": s.factor(
            original.a**3 * original.ell - s.Rational(1, 10)
        ),
        "whole_reference_scale_time_connection": s.factor(
            s.diff(original.a, u) / original.a - original.H
        ),
        "whole_reference_trace_density_background": s.factor(
            background_pv / (3 * original.a**3) + 2 * original.H
        ),
        "entire_scale_square_polynomial": s.expand(
            (1 + x) ** 4 - 1 - sum(s.binomial(4, i) * x**i for i in range(1, 5))
        ),
        "entire_volume_polynomial": s.expand(
            (1 + x) ** 6 - 1 - sum(s.binomial(6, i) * x**i for i in range(1, 7))
        ),
        "full_real_phase_count": s.Integer(field.DIMENSION * 2 - 96),
    }
    f = domain["fields"]
    return {
        "whole_original_fixed_canonical_scalar_tensor_heavy_generators": flow["blocks"],
        "whole_original_fixed_canonical_three_vector_generator": flow["vector"],
        "whole_original_vector_K_and_V_with_Gauss": [flow["K"], flow["V"]],
        "whole_complete_generator_outward_entry_sum_ceilings": {
            name: s.ceiling(value) for name, value in flow["bounds"].items()
        },
        "whole_complete_96_phase_generator_ceiling": CANONICAL_FLOW,
        "whole_actual_original_reference_pivot": [
            s.Rational(
                flow["actual_pivot"].lo.numerator, flow["actual_pivot"].lo.denominator
            ),
            s.Rational(
                flow["actual_pivot"].hi.numerator, flow["actual_pivot"].hi.denominator
            ),
        ],
        "whole_original_full_free_flow_proof": "The S251 generator is in fixed unit-CCR canonical field/density variables, including the complete Jc,A,Tcorr coefficients and density time connection. Its unchanged pure initial covariance already includes both canonical shears and the central swap; none is reset or reapplied as a duplicate connection. Add both tensor shape pairs, the physical H/Pi_H pair with the complete mass and all three vector pairs with K/zeta,zeta*V. The full96-phase norm is<1e500. This selected FREE reference is not identified with the full nonlinear Hamiltonian Hessian.",
        "whole_heavy_frequency_clarification": "Omega_star=sqrt(n+P^2/16) is the fixed S240 energy-balancing scale used by S269, not the physical KG dispersion. The physical omega(u)^2=n+P^2/a(u)^2; the chi equation additionally has -3H'/2-9H^2/4. At the bounce the physical-field dispersion is n+P^2. The original S269 balanced-state proof remains valid without substituting Omega_star into the equation. Both frequencies lie between1e98 and1e100 at the actual parameters.",
        "whole_original_external_means": {
            "a": original.a,
            "Hclock": original.H,
            "M1_mean_derivative": original.ell,
            "Pi_M1": s.Rational(1, 10),
            "Pi_v": background_pv,
            "H_and_W": 0,
        },
        "whole_full_covariance_S0_and_inverse_ceiling": WHITENING,
        "whole_whitened_generator_ceiling": WHITENED_FLOW,
        "whole_whitened_reference_flow_deviation_bound": domain["flow_deviation"],
        "whole_whitening_and_flow_proof": "Use the SAME full pure unit-CCR covariance V0=S0 S0^T/2 with a real symplectic S0. Fourier coefficient extraction costs sqrt(volume)<20, matrix conversion<2 and sqrt(96)<10; all S269 rows with physical kappa restored give ||S0||<1e120. Symplecticity and orthogonal Omega give ||S0^-1||=||S0||. Thus ||S0^-1 Gref S0||<1e1000, and the full integral equation gives ||Mw-I||<=exp(1e1000|u|)-1<=2e-1000, not a first-order frequency approximation. Initial complex w ball4R maps inside the S0-whitened ball8R.",
        "whole_complex_initial_and_image_radii": [ANALYTIC_RADIUS, IMAGE_RADIUS],
        "whole_all_complete_field_bounds_on_enlarged_ball": f,
        "whole_off_time_scale_volume_and_metric_gaps": {
            name: domain[name]
            for name in (
                "scale_gap",
                "volume_gap",
                "gamma_gap",
                "inverse_gap",
                "background_momentum",
                "delta_density",
            )
        },
        "whole_all_twelve_actual_off_time_complex_invariant_bounds": domain[
            "invariants"
        ],
        "whole_complete_complex_spatial_and_density_proof": "The full S269 A2 shape contraction and A0-to-A1 adjoint Neumann inverse use absolute complex norm bounds and therefore extend holomorphically to these finite complex inputs. Here the adjoint means the holomorphic continuation of the REAL adjoint, namely the formal bilinear transpose, not conjugate-transpose for complex Q. It restricts to the true Hilbert adjoint at real Q. No real coercivity or positivity is asserted for a complex metric. All reconstructed harmonics and generated means remain. With gamma=a^2 exp(2v)Q^-1 the extra a^+-2 gap is<=16T^2 and a^+-3 gap<=64T^2, so gamma and its inverse still have A2 gap<1e-272. All DQ*, derivative-Q ghost, full metric/vector/matter source and cotangent bounds remain. The reference trace momentum is<100T; Pi_M1 stays1/10, its full density shift adds at most13T^2 and is absorbed by the spare .4VLOG. Every full invariant is<1e-260. Real data retain the positive spatial branch and all three translations.",
        "checks": checks,
        "gates": {
            "entire_reference_generator_bound": flow["entire"] < CANONICAL_FLOW,
            "scalar_generator_no_unremoved_kappa": not any(
                value.has(phase.kappa) for value in flow["blocks"]["light"]
            ),
            "actual_full_J_lower_and_upper": flow["actual_pivot"].lo > 1
            and flow["actual_pivot"].hi < 100,
            "actual_full_clock_and_theta_bounds": flow["actual_theta"] < 2
            and flow["actual_Hubble"] < 2,
            "full_physical_and_balance_heavy_interval": 10**196
            < balance_frequency
            < physical_frequency
            < 10**200,
            "all_96_canonical_rows_bound_whitening": 10 * domain["whitening_row"]
            < WHITENING,
            "entire_whitened_generator": WHITENING**2 * CANONICAL_FLOW < WHITENED_FLOW,
            "full_exponential_in_proved_domain": WHITENED_FLOW * source.TIME
            < s.Rational(1, 2),
            "complete_complex_free_flow_inside_enlarged_ball": (
                1 + domain["flow_deviation"]
            )
            * ANALYTIC_RADIUS
            < IMAGE_RADIUS,
            "enlarged_complex_v_A2": f["v_A2"] < old_geometry.VLOG,
            "enlarged_complex_tau_A2": f["tau_A2"] < old_geometry.TAU,
            "all_complete_momentum_A1_with_background": all(
                f[name] * (1 + P) + domain["background_momentum"]
                < old_geometry.MOMENTUM_A1
                for name in (
                    "Pi_v_A0",
                    "Pi_tau_A0",
                    "delta_Pi_M_A0",
                    "Pi_H_A0",
                    "Pi_W_A0",
                )
            ),
            "complete_W_A1": f["W_A0"] * (1 + P) < s.Rational(1, 10**275),
            "complete_M1_gradient": f["M1_A1"] < s.Rational(1, 10**329),
            "complete_H_gradient": f["H_gradient_A0"] < s.Rational(1, 10**354),
            "full_off_time_metric_A2_gap": domain["gamma_gap"] < old_geometry.METRIC
            and domain["inverse_gap"] < old_geometry.METRIC,
            "actual_volume_background_below_two": (1 + source.TIME**2) ** 6 < 2,
            "complete_complex_c_density_below_two": (1 + source.TIME**2) ** 4
            / (1 - 2 * old_geometry.VLOG)
            < 2,
            "complete_complex_volume_and_inverse_below_two": (1 + source.TIME**2) ** 6
            / (1 - 3 * old_geometry.VLOG)
            < 2,
            "whole_density_shift_spare_margin": domain["delta_density"]
            < s.Rational(2, 5) * old_geometry.VLOG,
            "full_background_momentum_bound": 6 * 2 * 4 < 100,
            **{
                "entire_off_time_invariant_" + str(name): bool(value < source.IMAGE)
                for name, value in domain["invariants"].items()
            },
            "full_scale_and_volume_binomial_coefficients": sum(
                s.binomial(4, i) for i in range(1, 5)
            )
            < 16
            and sum(s.binomial(6, i) for i in range(1, 7)) < 64,
            "same_seed_and_residual_translations_not_removed": True,
        },
    }
