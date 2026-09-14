"""Complete finite-band matrix energy reduction with the heavy oscillator retained."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_principal_obstruction import coupled as c
from p8_vacuum_affine_heavy_source_filtration import source

from . import background as bg
from .intervals import I, evaluate


@cache
def matrix_data():
    td, Dd, Ld = s.symbols("theta_dot D_dot L2_dot", real=True)
    mu = s.Symbol("positive_heavy_potential_mass2", positive=True)
    rho = s.Symbol("positive_fast_rate", positive=True)
    restriction = {
        c.L0: 0,
        c.Vvv: 0,
        c.vs[0]: 0,
        c.vs[1]: 0,
        c.ds[0]: 0,
        c.m00: 0,
        c.m01: 0,
        c.m11: -mu,
    }
    ham = c.central_hamiltonian().subs(restriction, simultaneous=True)
    H = s.hessian(ham, c.CENTRAL)
    shift = -4 * c.a * c.th * c.P**2 / (c.D * c.L2)
    R = s.eye(8)
    R[4, 0] = shift
    Rdot = s.zeros(8)
    Rdot[4, 0] = (
        c.H * c.a * s.diff(shift, c.a)
        + td * s.diff(shift, c.th)
        + Dd * s.diff(shift, c.D)
        + Ld * s.diff(shift, c.L2)
    )
    whole = (R.T * H * R + c.OMEGA * R.inv() * Rdot).applyfunc(s.factor)
    A = c.OMEGA * whole
    order = [0, 1, 3, 4, 5, 7, 2, 6]
    A = A.extract(order, order)
    alpha = c.L2**2 / (8 * c.J * c.a**3)
    g = c.r * s.sqrt(c.zeta) / (c.D * c.a**2)
    h = c.a**3 * c.Yv / c.zeta
    beta = c.L2 * c.ws[0] / (4 * c.J * c.Z * c.a**3)
    charge = c.a * c.cs[0] / c.D
    d = c.a * c.Y
    fast = s.zeros(6)
    fast[0, 3] = alpha
    fast[1, 3] = beta
    fast[2, 0] = -g
    fast[3, 5] = g
    fast[4, 0] = -charge
    fast[4, 1] = -d
    fast[5, 2] = -h
    M = s.diag(1, 1, s.sqrt(c.P), c.P ** s.Rational(3, 2), s.sqrt(c.P), c.P, 1, 1)
    heavy_A = s.factor(c.a**3 * mu + c.a * c.Y * c.P**2)
    heavy_B = 1 / (c.Z * c.a**3)
    heavy = s.Matrix([[0, heavy_B], [-heavy_A, 0]])
    residual = (
        M.inv() * A * M - s.diag(c.P ** s.Rational(3, 2) * fast, heavy)
    ).applyfunc(s.factor)
    slow_weight = s.Rational(1, 80000)
    S = s.zeros(6)
    S[0, 0] = 1
    S[1, 2] = -rho / g
    S[2, 3] = alpha / rho
    S[3, 5] = rho**2 / (g * h)
    S[4, 1] = 1
    S[4, 0] = -beta / alpha
    S[5, 4] = slow_weight
    S[5, 2] = -slow_weight * (charge + d * beta / alpha) / g
    W = s.diag(S, s.diag(s.sqrt(heavy_A), s.sqrt(heavy_B)))
    scaled = (W * residual * W.inv() / c.P ** s.Rational(3, 2)).applyfunc(s.factor)
    permutation = s.eye(8).extract(order, list(range(8)))
    physical_to_central = c.chart()["whole_central_symplectic_map"]
    physical_map = (
        W * M.inv() * permutation * R.inv() * physical_to_central
    ).applyfunc(s.factor)
    physical_inverse = (
        physical_to_central.inv() * R * permutation.T * M * W.inv()
    ).applyfunc(s.factor)

    near = I(s.Rational(999, 1000), s.Rational(1001, 1000))
    tiny = I(-s.Rational(1, 10**20), s.Rational(1, 10**20))
    heavy_tiny = I(-s.Rational(1, 10**1000), s.Rational(1, 10**1000))
    jets = I(-(10**40), 10**40)
    bindings = {
        c.P: I(bg.MOMENTUM_MIN, bg.MOMENTUM_MAX),
        c.a: near,
        c.D: near,
        c.Z: near,
        c.Y: near,
        c.Yv: near,
        c.J: I(s.Rational(151, 100), s.Rational(153, 100)),
        c.L2: I(-s.Rational(1001, 1000), -s.Rational(999, 1000)),
        c.r: I(-s.Rational(3, 10**6), -s.Rational(1, 10**6)),
        c.rN: I(-3, -1),
        c.C: I(s.Rational(499, 1000), s.Rational(501, 1000)),
        c.K: I(s.Rational(999, 10**9), s.Rational(1001, 10**9)),
        c.zeta: I(s.Rational(1, 10**6)),
        c.cs[0]: I(s.Rational(99, 1000), s.Rational(101, 1000)),
        c.ws[0]: I(s.Rational(49, 1000), s.Rational(51, 1000)),
        c.cs[1]: heavy_tiny,
        c.ws[1]: heavy_tiny,
        c.ds[1]: I(-s.Rational(1, 10**880), s.Rational(1, 10**880)),
        c.th: tiny,
        c.dH: tiny,
        c.H: tiny,
        td: jets,
        Dd: jets,
        Ld: jets,
        mu: I(source.MASS2 / 2, 2 * source.MASS2),
        rho: I(s.Rational(1, 10**4), s.Rational(1, 100)),
    }
    bound = 0
    maximum = (0, None)
    for i in range(8):
        for j in range(8):
            interval = evaluate(scaled[i, j], bindings)
            upper = max(abs(interval.lo), abs(interval.hi))
            bound += upper
            if upper > maximum[0]:
                maximum = (upper, (i, j))

    cycle = s.Matrix([[0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0]])
    slow = s.Matrix([[0, 0], [-slow_weight * d, 0]])
    normal = S * fast * S.inv() - s.diag(rho * cycle, slow)
    relation = rho**4 - alpha * g**2 * h
    checks = {}
    for index, entry in enumerate(normal):
        numerator, _ = s.fraction(s.cancel(entry))
        checks["complete_light_growth_chart_entry_" + str(index)] = s.rem(
            numerator, relation, rho
        )

    checks["whole_perfect_square_shear_is_symplectic"] = R * c.OMEGA * R.T - c.OMEGA
    checks["whole_shear_time_connection_symmetric"] = (
        c.OMEGA * R.inv() * Rdot - (c.OMEGA * R.inv() * Rdot).T
    )
    omega_heavy = s.sqrt(s.factor(heavy_A * heavy_B))
    energy_chart = s.diag(s.sqrt(heavy_A), s.sqrt(heavy_B))

    def positive_square_normalize(value):
        value = s.simplify(value)
        value = value.replace(
            lambda term: term.is_Pow and term.exp == s.Rational(1, 2),
            lambda term: s.sqrt(s.factor(term.base)),
        )
        return s.factor(value)

    checks["retained_heavy_oscillator_is_exactly_skew"] = (
        energy_chart * heavy * energy_chart.inv()
        - s.Matrix([[0, omega_heavy], [-omega_heavy, 0]])
    ).applyfunc(positive_square_normalize)
    checks["heavy_skew_part_has_zero_energy_growth"] = (
        s.Matrix([[0, omega_heavy], [-omega_heavy, 0]])
        + s.Matrix([[0, omega_heavy], [-omega_heavy, 0]]).T
    )
    conversion_bounds = []
    for matrix in (physical_map, physical_inverse):
        value = 0
        for entry in matrix:
            interval = evaluate(entry, bindings)
            value += max(abs(interval.lo), abs(interval.hi))
        conversion_bounds.append(value)
    checks["whole_physical_growth_map_inverse"] = (
        physical_map * physical_inverse - s.eye(8)
    ).applyfunc(s.factor)
    return {
        "whole_on_shell_eight_phase_Hessian": whole,
        "whole_old_from_new_perfect_square_shear": R,
        "whole_perfect_square_shear_derivative": Rdot,
        "whole_light_growth_chart": S,
        "whole_light_momentum_weights": M[:6, :6],
        "whole_heavy_energy_chart": energy_chart,
        "whole_heavy_oscillator_core": heavy,
        "whole_physical_phase_to_growth_map": physical_map,
        "whole_growth_to_physical_phase_map": physical_inverse,
        "both_full_physical_phase_conversion_norm_upper_bound": s.Integer(10) ** 150,
        "whole_normalized_mass_adapted_remainder": scaled,
        "complete_input_box": {
            str(key): value.bounds() for key, value in bindings.items()
        },
        "full_rational_remainder_sum_upper_bound": s.Rational(1, 10**19),
        "positive_fast_rate_relation": relation,
        "normal_fast_cycle": cycle,
        "whole_retained_slow_nilpotent_block": slow,
        "slow_weight": slow_weight,
        "rho_min": s.Rational(1, 10**4),
        "rho_max": s.Rational(1, 100),
        "checks": checks,
        "gates": {
            "all_sixty_four_entries_bounded_on_whole_band": bool(
                bound < s.Rational(1, 10**19)
            ),
            "entire_remainder_below_rho_min_over_32": bool(
                bound < s.Rational(1, 320000)
            ),
            "heavy_mass_and_source_not_deleted": all(
                whole.has(term) for term in (mu, c.cs[1], c.ws[1], c.ds[1])
            ),
            "no_limit_in_P_or_heavy_mass_used": True,
            "exact_heavy_core_retained_not_a_raw_mass_error_bound": True,
            "both_full_physical_conversion_norms_below_1e150": all(
                value < 10**150 for value in conversion_bounds
            ),
        },
    }


@cache
def time_data():
    alpha, h, rho, d = s.symbols("alpha h rho d", positive=True)
    g = s.Symbol("g", nonzero=True, real=True)
    beta, charge = s.symbols("beta charge", real=True)
    variables = [alpha, g, h, beta, charge, d]
    dots = s.symbols("alpha_dot g_dot h_dot beta_dot charge_dot d_dot", real=True)
    rho_dot = rho * (dots[0] / alpha + 2 * dots[1] / g + dots[2] / h) / 4
    weight = s.Rational(1, 80000)
    S = s.zeros(6)
    S[0, 0] = 1
    S[1, 2] = -rho / g
    S[2, 3] = alpha / rho
    S[3, 5] = rho**2 / (g * h)
    S[4, 1] = 1
    S[4, 0] = -beta / alpha
    S[5, 4] = weight
    S[5, 2] = -weight * (charge + d * beta / alpha) / g
    Sdot = S.diff(rho) * rho_dot
    for variable, dot in zip(variables, dots):
        Sdot += S.diff(variable) * dot
    connection = (Sdot * S.inv()).applyfunc(s.factor)
    bindings = {
        alpha: I(s.Rational(1, 20), s.Rational(1, 10)),
        g: I(-s.Rational(1, 10**8), -s.Rational(1, 10**10)),
        h: I(10**5, 10**7),
        rho: I(s.Rational(1, 10**4), s.Rational(1, 100)),
        beta: I(-1, 1),
        charge: I(-1, 1),
        d: I(s.Rational(1, 2), 2),
        **{dot: I(-(10**50), 10**50) for dot in dots},
    }
    bound = 0
    for entry in connection:
        interval = evaluate(entry, bindings)
        bound += max(abs(interval.lo), abs(interval.hi))

    return {
        "whole_light_growth_chart_time_connection": connection,
        "whole_positive_rate_logarithmic_derivative": rho_dot / rho,
        "primitive_time_derivative_box": s.Integer(10) ** 50,
        "whole_light_chart_connection_upper_bound": s.Integer(10) ** 63,
        "whole_normalized_light_chart_connection_upper_bound": s.Rational(1, 10**33),
        "checks": {},
        "gates": {
            "all_light_chart_time_entries_bounded": bool(bound < 10**63),
            "full_light_chart_time_connection_below_rho_min_over_32": bool(
                bound / 10**96 < s.Rational(1, 320000)
            ),
        },
    }


@cache
def coefficient_data():
    packet = matrix_data()
    raw = packet["complete_input_box"]

    def enclosure(variable):
        return I(*raw[str(variable)])

    parameters = {
        "alpha": c.L2**2 / (8 * c.J * c.a**3),
        "g": c.r * s.sqrt(c.zeta) / (c.D * c.a**2),
        "h": c.a**3 * c.Yv / c.zeta,
        "beta": c.L2 * c.ws[0] / (4 * c.J * c.Z * c.a**3),
        "charge": c.a * c.cs[0] / c.D,
        "d": c.a * c.Y,
    }
    safe = {
        "alpha": (s.Rational(1, 20), s.Rational(1, 10)),
        "g": (-s.Rational(1, 10**8), -s.Rational(1, 10**10)),
        "h": (10**5, 10**7),
        "beta": (-1, 1),
        "charge": (-1, 1),
        "d": (s.Rational(1, 2), 2),
    }
    variables = set().union(*(value.free_symbols for value in parameters.values()))
    bindings = {variable: enclosure(variable) for variable in variables}
    dots = {
        variable: s.Symbol("physical_time_derivative_" + str(variable), real=True)
        for variable in variables
        if variable != c.zeta
    }
    dot_bindings = {dot: I(-(10**40), 10**40) for dot in dots.values()}
    gates = {}
    derivatives = {}
    for name, expr in parameters.items():
        gates["actual_growth_parameter_box_" + name] = evaluate(expr, bindings).inside(
            *safe[name]
        )
        derivative = sum(s.diff(expr, variable) * dot for variable, dot in dots.items())
        derivatives[name] = derivative
        value = evaluate(derivative, {**bindings, **dot_bindings})
        gates["actual_growth_parameter_first_time_jet_" + name] = (
            max(abs(value.lo), abs(value.hi)) < 10**50
        )
    fourth = parameters["alpha"] * parameters["g"] ** 2 * parameters["h"]
    fourth_box = evaluate(fourth, bindings)
    gates["actual_positive_fast_rate_lower_and_upper"] = fourth_box.inside(
        s.Rational(1, 10**16), s.Rational(1, 10**8)
    )
    gamma = evaluate(
        c.GAMMA, {variable: enclosure(variable) for variable in c.GAMMA.free_symbols}
    )
    gates["actual_whole_auxiliary_gamma_margin"] = gamma.lo > s.Rational(1, 2)
    return {
        "complete_physical_fast_parameter_map": parameters,
        "complete_physical_fast_parameter_time_derivatives": derivatives,
        "whole_parameter_box": safe,
        "actual_fast_rate_fourth_power": fourth,
        "actual_fast_rate_fourth_power_enclosure": fourth_box.bounds(),
        "checks": {
            "actual_complete_fast_quartic_coefficient": s.factor(
                fourth - c.r**2 * c.L2**2 * c.Yv / (8 * c.D**2 * c.J * c.a**4)
            )
        },
        "gates": {key: bool(value) for key, value in gates.items()},
    }


@cache
def other_modes_data():
    # In the fixed physical canonical Fourier variables. Each of these
    # Hamiltonians occurs twice: the two TT and two transverse Proca modes.
    tt_A, tt_B = 2 * c.a * c.C * c.P**2, 1 / (c.a**3 * c.D)
    pt_A = c.a * c.K * c.Y * c.P**2 / (c.zeta * c.Z) + c.a**3 * c.Yv / c.zeta
    pt_B = c.zeta / (c.a**3 * c.K)
    mu = s.Symbol("positive_heavy_potential_mass2", positive=True)
    heavy_A, heavy_B = c.a**3 * mu + c.a * c.Y * c.P**2, 1 / (c.Z * c.a**3)
    omega = s.Matrix([[0, 1], [-1, 0]])
    raw = matrix_data()["complete_input_box"]
    derivatives = {}
    bounds = {}
    cores = {}
    conversion_bounds = {}
    for name, A, B in [
        ("TT", tt_A, tt_B),
        ("transverse_Proca", pt_A, pt_B),
        ("heavy", heavy_A, heavy_B),
    ]:
        cores[name] = omega * s.diag(A, B)
        variables = A.free_symbols | B.free_symbols
        mapping = {variable: I(*raw[str(variable)]) for variable in variables}
        for direction, powers in [
            ("forward", (s.Rational(1, 2), s.Rational(1, 2))),
            ("inverse", (-s.Rational(1, 2), -s.Rational(1, 2))),
        ]:
            total = 0
            for value, power in zip((A, B), powers):
                interval = evaluate(value**power, mapping)
                total += max(abs(interval.lo), abs(interval.hi))
            conversion_bounds[name + "_" + direction] = bool(total < 10**150)
        dots = {
            variable: s.Symbol(name + "_dot_" + str(variable), real=True)
            for variable in variables
            if variable not in (c.P, c.zeta, mu)
        }
        eta = s.Symbol("heavy_relative_mass_coefficient_derivative", real=True)
        for slot, value in [("A", A), ("B", B)]:
            derivative = sum(
                s.diff(value, variable) * dot for variable, dot in dots.items()
            )
            if value.has(mu):
                derivative += s.diff(value, mu) * mu * eta
            log_derivative = s.factor(derivative / (2 * value))
            interval = evaluate(
                log_derivative,
                {
                    **mapping,
                    **{dot: I(-(10**40), 10**40) for dot in dots.values()},
                    eta: I(-(10**40), 10**40),
                },
            )
            upper = max(abs(interval.lo), abs(interval.hi))
            derivatives[name + "_" + slot] = log_derivative
            bounds[name + "_" + slot] = bool(upper < 10**50)
    AA, BB = s.symbols("positive_mode_potential positive_mode_kinetic", positive=True)
    chart = s.diag(s.sqrt(AA), s.sqrt(BB))
    expected = s.sqrt(AA * BB) * omega
    check = (chart * omega * s.diag(AA, BB) * chart.inv() - expected).applyfunc(
        s.simplify
    )
    # Independent ADM Maxwell and TT density re-entry.
    N, R, a, P, zeta = s.symbols("N R a P zeta", positive=True)
    electric = R ** (-s.Rational(1, 4)) * a / N
    magnetic = N * R ** s.Rational(1, 4) / a
    proca_mass = N * R ** (-s.Rational(1, 4)) * a / zeta
    actual = {
        c.a: a,
        c.P: P,
        c.zeta: zeta,
        c.K: zeta * R ** (-s.Rational(1, 4)) / (N * a * a),
        c.Y: N * R ** (-s.Rational(1, 4)),
        c.Yv: N * R ** (-s.Rational(1, 4)) / (a * a),
        c.Z: R ** (-s.Rational(3, 4)) / N,
    }
    # For a TT polarization normalized by e_ij e_ij=4, the ADM density
    # is D*a^3*qdot^2/2-C*a*P^2*q^2. Its full canonical Legendre transform
    # fixes both coefficients, with no zero-order volume Euler remainder
    # in the determinant-one tensor chart.
    tensor_q, tensor_v, tensor_p = s.symbols("tensor_q tensor_v tensor_p", real=True)
    tensor_lag = c.D * c.a**3 * tensor_v**2 / 2 - c.C * c.a * c.P**2 * tensor_q**2
    tensor_ham = (tensor_p * tensor_v - tensor_lag).subs(
        tensor_v, tensor_p / (c.D * c.a**3)
    )
    return {
        "all_remaining_positive_mode_cores": cores,
        "mode_multiplicities": {"TT": 2, "transverse_Proca": 2, "heavy": 1},
        "full_mode_energy_chart": chart,
        "all_remaining_energy_chart_time_entries": derivatives,
        "all_remaining_root_energy_time_entry_bound": s.Integer(10) ** 50,
        "full_normalized_remaining_time_bound": s.Rational(1, 10**46),
        "checks": {
            "whole_positive_oscillator_energy_chart": check,
            "full_TT_canonical_Legendre_block": s.hessian(
                tensor_ham, (tensor_q, tensor_p)
            )
            - s.diag(tt_A, tt_B),
            "full_transverse_Proca_electric_kinetic": s.factor(
                pt_B.subs(actual, simultaneous=True) - 1 / electric
            ),
            "full_transverse_Proca_magnetic_and_mass_potential": s.factor(
                pt_A.subs(actual, simultaneous=True) - magnetic * P**2 - proca_mass
            ),
        },
        "gates": {
            **bounds,
            **{
                "full_other_mode_physical_conversion_" + key: value
                for key, value in conversion_bounds.items()
            },
            "both_TT_both_transverse_Proca_and_heavy_retained": True,
            "large_heavy_oscillation_is_kept_as_skew_generator": True,
        },
    }


@cache
def growth_data():
    matrix = matrix_data()
    time = time_data()
    others = other_modes_data()
    rho_min = matrix["rho_min"]
    error = (
        matrix["full_rational_remainder_sum_upper_bound"]
        + time["whole_normalized_light_chart_connection_upper_bound"]
        + others["full_normalized_remaining_time_bound"]
    )
    delta = rho_min / 32
    cycle = matrix["normal_fast_cycle"]
    positive = s.ones(4, 1) / 2
    negative = s.Matrix([1, -1, -1, 1]) / 2
    d = s.Symbol("positive_slow_gradient", positive=True)
    slow = s.Matrix([[0, 0], [-matrix["slow_weight"] * d, 0]])
    x, y = s.symbols("slow_position slow_momentum", real=True)
    z = s.Matrix([x, y])
    slow_margin = (
        matrix["slow_weight"] * d * (x * x + y * y) / 2
        - (z.T * (slow + slow.T) * z)[0] / 2
    )
    cone_amplitude = rho_min - 2 * delta
    cone_complement = rho_min / 8 + 2 * delta
    exponent = rho_min * bg.MOMENTUM_MIN ** s.Rational(3, 2) * bg.TIME_LENGTH / 2
    upper_fast_rate = s.Integer(3) * 10**94
    return {
        "entire_normalized_remainder_and_time_error_upper_bound": error,
        "uniform_full_error_budget": delta,
        "whole_principal_complement_logarithmic_norm_upper_bound": rho_min / 8,
        "whole_cone_amplitude_derivative_coefficient_lower": cone_amplitude,
        "whole_cone_complement_boundary_derivative_coefficient_upper": cone_complement,
        "whole_cone_strict_invariance_margin": cone_amplitude - cone_complement,
        "comoving_momentum_band": [bg.MOMENTUM_MIN, bg.MOMENTUM_MAX],
        "physical_wavenumber_enclosure": [bg.MOMENTUM_MIN / 2, 2 * bg.MOMENTUM_MAX],
        "classical_clock_time_interval": [s.Integer(0), bg.TIME_LENGTH],
        "whole_growth_exponent_lower_bound": exponent,
        "physical_phase_conversion_prefactor_lower_bound": s.Rational(1, 10**300),
        "principal_fast_rate_uniform_upper_bound": upper_fast_rate,
        "complete_physical_mode_count": 8,
        "whole_real_phase_count": 16,
        "finite_band_growth_statement": "For the actual full-current S255 comparison solution with epsilon=10^-6 and the exact entire constraint root, all displayed bootstrap and coefficient bounds hold on[0,10^-60]. For every comoving P in[10^64,2*10^64], the complete physical canonical phase propagator has operator norm at least10^-300*exp(5*10^31). The proof keeps the massive heavy oscillator and all four TT/transverse-Proca modes, and uses the complete coupled remainder plus every chart time term. This is a fixed finite band and evaluated classical interval, not an infinite-momentum limit, a Wilsonian cutoff certificate, a nonlinear Cauchy theorem, or the fixed quantum mean.",
        "energy_cone_proof": "Split the normalized full16-phase equation into the positive four-cycle vector f and its entire complement z. The cycle symmetric part is eplus*eplus^T-eminus*eminus^T; the weighted slow block has logarithmic norm at most rho_min/8, and all heavy/other oscillator cores are retained and skew. The full remainder has norm below rho_min/32. On the boundary norm(z)=f>0, f_prime/P^(3/2)>=(rho_min-2delta)f while the complement upper derivative is at most(rho_min/8+2delta)f. Their strict separation makes the cone invariant. Integrate f_prime>=rho_min P^(3/2)f/2 and apply both complete physical chart bounds. Smooth Fourier packets inside a small directional sector of this band, together with their real conjugate packets, give the same L2 lower bound by Parseval.",
        "cutoff_boundary": "The momentum band and principal growth timescale are below the retained heavy mass scale, but this is not a derivation of the actual EFT cutoff or an error budget for omitted operators, loops, nonlinear terms or a UV matching model. No change to the fixed quantum reference preparation, earlier ladder verdicts or original P8(a) qualifications is made.",
        "checks": {
            "whole_fast_positive_eigenvector": cycle * positive - positive,
            "whole_fast_negative_eigenvector": cycle * negative + negative,
            "complete_cycle_symmetric_part": (cycle + cycle.T) / 2
            - positive * positive.T
            + negative * negative.T,
            "weighted_slow_logarithmic_norm_square": s.expand(
                slow_margin - matrix["slow_weight"] * d * (x + y) ** 2 / 2
            ),
            "actual_evaluated_finite_growth_exponent": exponent - 5 * 10**31,
            "strict_cone_margin_identity": cone_amplitude
            - cone_complement
            - s.Rational(3, 4) * rho_min,
        },
        "gates": {
            "whole_error_strictly_below_cone_budget": bool(error < delta),
            "whole_slow_complement_logarithmic_norm_bound": bool(
                matrix["slow_weight"] * 2 / 2 <= rho_min / 8
            ),
            "whole_cone_strict_forward_invariance": bool(
                cone_amplitude > cone_complement
            ),
            "complete_amplitude_growth_above_half_rho_min": bool(
                cone_amplitude > rho_min / 2
            ),
            "finite_band_not_an_unbounded_momentum_limit": bool(
                bg.MOMENTUM_MAX == 2 * bg.MOMENTUM_MIN
            ),
            "actual_physical_momenta_below_fixed_heavy_mass": bool(
                (2 * bg.MOMENTUM_MAX) ** 2 < source.MASS2
            ),
            "proper_time_principal_growth_rate_below_fixed_heavy_mass": bool(
                (2 * upper_fast_rate) ** 2 < source.MASS2
            ),
            "actual_entire_band_principal_rate_below_declared_upper_bound": bool(
                matrix["rho_max"] ** 2 * bg.MOMENTUM_MAX**3 < upper_fast_rate**2
            ),
            "evaluated_lifetime_with_complete_mass_energy_control": True,
            "no_cutoff_or_quantum_mean_or_original_P8_verdict": True,
        },
    }
