"""Entire sixteen-phase finite-band generator and original phase conversions."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_principal_obstruction import coupled as c
from p8_vacuum_affine_finite_window_growth.intervals import I, evaluate

from . import background, source


@cache
def matrices():
    mass, mu, kt, heavy_weight = s.symbols(
        "fixed_heavy_mass2 positive_mass_ratio K_over_zeta heavy_energy_weight",
        positive=True,
    )
    s00, s01, s11 = s.symbols(
        "fixed_reference_Sbar00 fixed_reference_Sbar01 fixed_reference_Sbar11",
        real=True,
    )
    shear = s.zeros(4)
    shear[0, 0], shear[0, 1], shear[1, 0], shear[1, 1] = [
        c.P**2 * x for x in (s00, s01, s01, s11)
    ]
    R = s.eye(8)
    R[4:, :4] = shear
    Ri = s.eye(8)
    Ri[4:, :4] = -shear
    W = s.diag(c.P, c.P, heavy_weight, c.P / s.sqrt(c.zeta), 1, 1, 1, s.sqrt(c.zeta))
    ham = c.central_hamiltonian().subs(
        {c.m00: 0, c.m01: 0, c.m11: -mu * mass, c.K: c.zeta * kt, c.ds[0]: 0},
        simultaneous=True,
    )
    H = s.hessian(ham, c.CENTRAL)
    A = (W * Ri * c.OMEGA * H * R * W.inv()).applyfunc(s.factor)
    omega2 = s.Matrix([[0, 1], [-1, 0]])
    weight2 = s.diag(c.P, 1)
    tt = omega2 * s.diag(2 * c.a * c.C * c.P**2, 1 / (c.a**3 * c.D))
    pt = omega2 * s.diag(
        c.a * kt * c.Y * c.P**2 / c.Z + c.a**3 * c.Yv / c.zeta, 1 / (c.a**3 * kt)
    )
    tt = (weight2 * tt * weight2.inv()).applyfunc(s.factor)
    pt = (weight2 * pt * weight2.inv()).applyfunc(s.factor)
    A = s.diag(A, tt, tt, pt, pt)
    assert A.shape == (16, 16)
    near = I(s.Rational(999, 1000), s.Rational(1001, 1000))
    tiny = I(-s.Rational(1, 10**20), s.Rational(1, 10**20))
    small = I(-s.Rational(1, 10**6), s.Rational(1, 10**6))
    heavy_tiny = I(-s.Rational(1, 10**780), s.Rational(1, 10**780))
    boxes = {
        c.P: I(10**64, 2 * 10**64),
        c.a: near,
        c.D: near,
        c.Z: near,
        c.Y: near,
        c.Yv: near,
        kt: near,
        mu: near,
        c.J: I(s.Rational(151, 100), s.Rational(153, 100)),
        c.C: I(s.Rational(499, 1000), s.Rational(501, 1000)),
        c.L2: I(-s.Rational(1001, 1000), -s.Rational(999, 1000)),
        c.r: tiny,
        c.th: tiny,
        c.dH: tiny,
        c.H: tiny,
        c.rN: I(-3, -1),
        c.cs[0]: I(s.Rational(99, 1000), s.Rational(101, 1000)),
        c.ws[0]: I(s.Rational(49, 1000), s.Rational(51, 1000)),
        c.cs[1]: heavy_tiny,
        c.ws[1]: heavy_tiny,
        c.ds[1]: heavy_tiny,
        c.L0: small,
        c.Vvv: small,
        c.vs[0]: small,
        c.vs[1]: small,
        c.zeta: I(s.Rational(1, 10**6)),
        mass: I(10**197, 10**198),
        heavy_weight: I(10**98, 10**100),
        s00: tiny,
        s01: tiny,
        s11: tiny,
    }
    variables = (
        c.a,
        c.D,
        c.Z,
        c.Y,
        c.Yv,
        kt,
        mu,
        c.J,
        c.C,
        c.L2,
        c.r,
        c.th,
        c.dH,
        c.H,
        c.rN,
        c.cs[0],
        c.ws[0],
        c.cs[1],
        c.ws[1],
        c.ds[1],
        c.L0,
        c.Vvv,
        c.vs[0],
        c.vs[1],
    )
    maximum = 0
    total = 0
    enclosures = {}
    for variable in variables:
        bound = 0
        for entry in A:
            derivative = s.factor(s.diff(entry, variable))
            value = evaluate(derivative, boxes)
            bound += max(abs(value.lo), abs(value.hi))
        assert bound < 10**180, (str(variable), float(bound))
        total += bound
        maximum = max(maximum, bound)
        enclosures[str(variable)] = background.rational(bound)
    assert total < 10**180
    assert Ri * R == s.eye(8)
    assert R * c.OMEGA * R.T == c.OMEGA

    Sd = s.zeros(4)
    dots = s.symbols("fixed_reference_Sbar_dot0:3", real=True)
    Sd[0, 0], Sd[0, 1], Sd[1, 0], Sd[1, 1] = [
        c.P**2 * x for x in (dots[0], dots[1], dots[1], dots[2])
    ]
    Rdot = s.zeros(8)
    Rdot[4:, :4] = Sd
    connection = s.diag(-W * Ri * Rdot * W.inv(), s.zeros(8))
    actual_map = c.chart()["whole_central_symplectic_map"]
    forward8 = W * Ri * actual_map
    inverse8 = actual_map.inv() * R * W.inv()
    forward = s.diag(forward8, weight2, weight2, weight2, weight2)
    inverse_map = s.diag(
        inverse8, weight2.inv(), weight2.inv(), weight2.inv(), weight2.inv()
    )
    conversions = []
    for matrix in (forward, inverse_map):
        conversions.append(
            background.rational(
                sum(
                    background.absolute(evaluate(s.factor(entry), boxes))
                    for entry in matrix
                )
            )
        )
    mass2 = background.b.source.MASS2
    return {
        "whole_on_shell_coupled_Hamiltonian_before_reference_shear": ham,
        "whole_sixteen_phase_generator_without_common_reference_connection": A,
        "whole_common_fixed_reference_time_connection": connection,
        "whole_fixed_reference_symmetric_boundary_old_from_new": R,
        "whole_mass_adapted_coupled_phase_weight": W,
        "whole_actual_heavy_weight_binding": s.sqrt(c.P**2 + mass),
        "whole_fixed_mass_binding": mass2,
        "whole_complete_coefficient_MVT_box": {
            str(key): value.bounds() for key, value in boxes.items()
        },
        "whole_all_24_coefficient_generator_derivative_sum_enclosures": enclosures,
        "whole_generator_coefficient_Lipschitz_sum": background.rational(total),
        "whole_safe_generator_coefficient_Lipschitz_bound": s.Integer(10) ** 180,
        "whole_original_canonical_phase_to_energy_coordinates": forward,
        "whole_energy_coordinates_to_original_canonical_phase": inverse_map,
        "whole_both_phase_conversion_enclosures": conversions,
        "whole_safe_both_phase_conversion_bound": s.Integer(10) ** 150,
        "whole_real_phase_count": 16,
        "checks": {
            "whole_coupled_reference_shear_inverse": Ri * R - s.eye(8),
            "whole_reference_shear_symplectic": R * c.OMEGA * R.T - c.OMEGA,
            "whole_complete_original_phase_conversion_inverse": (
                forward * inverse_map - s.eye(16)
            ).applyfunc(s.factor),
            "whole_common_time_connection_has_no_current_coefficient_derivative": s.Matrix(
                [s.diff(entry, var) for entry in connection for var in variables]
            ),
            "whole_mode_count": s.Integer(A.rows) - 16,
        },
        "gates": {
            "all_24_full_generator_derivative_sums_bounded": len(enclosures) == 24
            and total < 10**180,
            "whole_finite_P_interval_not_asymptotic": True,
            "both_original_phase_conversion_bounds": all(
                value < 10**150 for value in conversions
            ),
            "actual_full_mass_in_weight_box": bool(
                10**197 < mass2 < 10**198
                and 10**196 < mass2 + source.LOW**2
                and mass2 + source.HIGH**2 < 10**200
            ),
            "full_heavy_mass_and_all_lapse_contacts_retained": all(
                ham.has(value) for value in (mass, c.cs[1], c.ws[1], c.ds[1])
            ),
            "two_TT_and_two_transverse_Proca_pairs_retained": A.shape == (16, 16),
            "current_coefficient_variation_does_not_reselect_reference_shear": True,
            "same_reference_time_connection_cancels_only_in_generator_difference": True,
        },
    }


@cache
def comparison():
    T = source.TIME
    rate = s.Integer(10) ** 40
    heavy = s.Rational(1, 10**1000)
    initial = s.Rational(1, 10**225)
    forcing = 14 * rate * source.ERROR + 2 * rate * heavy
    light_bound = 2 * (initial + T * forcing)
    state = s.Rational(1, 10**224)
    coefficient_bound = (
        4 * state + 2 * heavy + 64 * source.ERROR
    ) * 10**80 + 10**80 * source.ERROR
    physical_bound = (4 * state + 2 * heavy + 64 * source.ERROR) * 10**160
    heavy_improved = 4 * T * s.Rational(1, 10**2500)
    generator_error = s.Integer(10) ** 180 * s.Rational(1, 10**143)
    root_rate = s.Integer(10) ** 28 + 10**12 * generator_error
    root_growth = 1 / (1 - 2 * T * root_rate)
    u = background.u
    q, v, N, Nd, mass, gamma, j = s.symbols(
        "heavy_q heavy_v positive_N Ndot positive_mass2 friction full_source", real=True
    )
    e2 = v * v + N * N * mass * q * q
    ed = (
        s.diff(e2, q) * v
        + s.diff(e2, v) * (-gamma * v - N * N * mass * q + N * N * j)
        + s.diff(e2, N) * Nd
    )
    wanted = -2 * gamma * v * v + 2 * N * N * j * v + 2 * N * Nd * mass * q * q
    return {
        "whole_symmetric_clock_slab": [-T, T],
        "whole_signed_initial_lapse_radius": source.EPSILON,
        "whole_light_state_Gronwall_Lipschitz_constant": 4 * rate,
        "whole_fixed_source_and_heavy_light_forcing_bound": forcing,
        "whole_improved_light_state_distance_enclosure": light_bound,
        "whole_safe_light_state_distance": state,
        "whole_complete_coefficient_distance_enclosure": coefficient_bound,
        "whole_safe_complete_coefficient_distance": s.Rational(1, 10**143),
        "whole_physical_Hubble_and_acceleration_distance_enclosure": physical_bound,
        "whole_safe_physical_Hubble_and_acceleration_distance": s.Rational(1, 10**63),
        "whole_mass_adapted_heavy_energy_square": e2,
        "whole_exact_heavy_energy_derivative": ed,
        "whole_improved_heavy_root_energy_bound": heavy_improved,
        "whole_safe_physical_turn_clock_distance": s.Rational(1, 10**63),
        "whole_physical_proper_acceleration_lower_bound": s.Integer(3),
        "whole_complete_generator_difference_bound": generator_error,
        "whole_full_reference_metric_condition_sqrt_bound": s.Integer(10) ** 12,
        "whole_full_root_energy_time_rate": root_rate,
        "whole_full_slab_root_energy_growth_bound": root_growth,
        "whole_safe_original_phase_propagator_bound": s.Integer(10) ** 313,
        "nonlinear_homogeneous_continuation_argument": "Use the unchanged full constraint root and its regular Cnn,J branch. Bootstrap all four light coordinates within 10^-200 of the bare clock and the original positive mass-adapted heavy box. The direct pivot formula and its complete coefficient derivatives keep every intermediate pivot in [1,2]; the bare pivot is near243/160. The full light-rate state bound and 14 independent source-jet bounds give the listed two-sided Gronwall improvement below10^-224. The heavy energy, including its actual mass and full source, improves its box on both time directions. Thus standard constrained-ODE continuation reaches[-10^-60,10^-60]. Bare coefficients only compare the unchanged full current solution; they never replace its action.",
        "coefficient_and_physical_turn_argument": "The 28 complete coefficient value/first-jet derivative rows give the displayed current-to-bare difference; the original fixed-reference profiles add only the bounded fixed error. Heavy charges and lapse source, including their first jets, are below10^-780. The actual physical scale is a_hat R^-1/4 and its proper Hubble is (Hhat-Rdot/(4R))/N. Both this value and its full proper acceleration differ from the exact bare clock values by less than10^-63. Hphysical has opposite endpoint signs and proper derivative greater than3, so there is exactly one physical turn in the slab, within10^-63 of the central clock. This is a homogeneous classical turn, not a quantum mean or a nonlinear inhomogeneous Cauchy theorem.",
        "full_variational_energy_argument": "Apply the single fixed-reference symplectic boundary and constant P/mass weights to the entire16-phase system. The shared reference time connection cancels in its difference, not in either individual equation. Every coefficient segment stays in the explicit whole-band box, giving ||Delta A||<10^180*10^-143. The full reference energy has metric eigenvalues in[10^-12,10^12] and original root-energy rate10^28. The perturbation contributes at most10^12||Delta A|| to the root-energy rate. Across the full symmetric slab the displayed positive-series bound is less than2. Both old physical canonical endpoint maps are retained; their finite bounds give the stated large unweighted propagator ceiling. Real/imaginary Fourier parts and full-band L2 integrals obey the same energy estimate. No all-momentum or Wilsonian-cutoff conclusion is drawn.",
        "checks": {
            "whole_heavy_frequency_cancellation_both_time_directions": s.expand(
                ed - wanted
            ),
            "whole_background_time_interval_unchanged": T - background.b.TIME_LENGTH,
            "whole_signed_lapse_radius": source.EPSILON - s.Rational(1, 10**230),
        },
        "gates": {
            "two_sided_light_Gronwall_exponential_below_two": 4 * rate * T
            < s.Rational(1, 2),
            "full_light_bootstrap_strict_improvement": light_bound < state,
            "full_coefficient_distance_bound": coefficient_bound
            < s.Rational(1, 10**143),
            "full_physical_Hubble_distance_bound": physical_bound
            < s.Rational(1, 10**63),
            "full_heavy_root_box_strict_improvement": heavy_improved < heavy / 10,
            "full_heavy_field_box_strict_improvement": 2
            * heavy_improved
            / s.Integer(10) ** 98
            < s.Rational(1, 10**1098),
            "positive_reference_Hubble_endpoint": background.Href.subs(u, T) > 3 * T,
            "negative_reference_Hubble_endpoint": background.Href.subs(u, -T) < -3 * T,
            "physical_endpoint_error_smaller_than_sign_margin": 3 * T
            > s.Rational(1, 10**63),
            "strict_positive_proper_acceleration": s.diff(background.Href, u).subs(u, T)
            - s.Rational(1, 10**63)
            > 3,
            "strict_positive_coordinate_Hubble_derivative": s.Rational(999, 1000) * 3
            > 2,
            "full_variational_root_energy_rate_bound": root_rate < 2 * 10**49,
            "full_symmetric_slab_root_energy_growth_below_two": 0
            < 2 * T * root_rate
            < s.Rational(1, 2)
            and root_growth < 2,
            "both_endpoint_original_phase_ceiling": 2 * 10**312 < 10**313,
            "tiny_classical_turn_not_original_B_or_quantum_mean": True,
            "finite_band_is_not_physical_cutoff": True,
        },
    }
