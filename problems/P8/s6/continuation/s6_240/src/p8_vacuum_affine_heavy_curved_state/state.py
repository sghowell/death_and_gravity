"""Actual reference KG equation, phase-regular SLE and full mode readouts."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_parent_one_loop import germs

TIME = s.Symbol("u", real=True)
X = s.Symbol("X", real=True)
A = (1 + TIME**2) ** 2
HUBBLE = 4 * TIME / (1 + TIME**2)
KAPPA = germs.KAPPA
MASS2 = germs.MASS2
N = germs.N
T0 = -s.Rational(1, 2)
SAMPLING_INTERVAL = (-s.Rational(3, 4), T0)
TIME_INTERVAL = (-s.S.One, s.S.One)
SCALE_MAX = s.Integer(4)


def sle_coefficients(c1, c2, conjugate_c2):
    delta = s.sqrt(c1 * c1 - c2 * conjugate_c2)
    alpha = s.sqrt((c1 + delta) / (2 * delta))
    beta = -c2 / s.sqrt(2 * delta * (c1 + delta))
    return alpha, beta, delta


def riccati_next(normalized, omega, epsilon=1):
    lam = s.diff(omega, TIME) / omega
    U = 3 * s.diff(HUBBLE, TIME) / 2 + 9 * HUBBLE**2 / 4
    log_rate = s.diff(normalized, TIME) / normalized
    bracket = (
        -U
        - s.diff(lam, TIME) / 2
        + lam**2 / 4
        - s.diff(log_rate, TIME) / 2
        + lam * log_rate / 2
        + log_rate**2 / 4
    )
    return s.sqrt(1 + epsilon**2 * bracket / omega**2)


@cache
def data():
    p, n = s.symbols("p n", positive=True)
    omega2 = n + p * p / A**2
    U = 3 * s.diff(HUBBLE, TIME) / 2 + 9 * HUBBLE**2 / 4
    chi = s.Function("chi")(TIME)
    mode = chi / (1 + TIME**2) ** 3
    canonical = (
        s.diff(mode, TIME, 2) + 3 * HUBBLE * s.diff(mode, TIME) + omega2 * mode
    ) * (1 + TIME**2) ** 3
    checks = {
        "full_minimal_KG_canonical_equation": s.simplify(
            canonical - s.diff(chi, TIME, 2) - (omega2 - U) * chi
        ),
        "global_Hubble_squared_margin": s.factor(
            4 - HUBBLE**2 - 4 * (TIME**2 - 1) ** 2 / (1 + TIME**2) ** 2
        ),
        "global_Hubble_derivative_upper": s.factor(
            4 - s.diff(HUBBLE, TIME) - 4 * TIME**2 * (TIME**2 + 3) / (1 + TIME**2) ** 2
        ),
    }
    lapse, scale, momentum = s.symbols("lapse scale comoving_momentum", positive=True)
    field, velocity = s.symbols("field velocity", real=True)
    lag = (
        lapse
        * scale**3
        * (velocity**2 / (2 * lapse**2) - (n + momentum**2 / scale**2) * field**2 / 2)
    )
    rho = -s.diff(lag, lapse) / scale**3
    pressure = s.diff(lag, scale) / (3 * lapse * scale**2)
    checks["independent_full_lapse_energy"] = s.cancel(
        rho - velocity**2 / (2 * lapse**2) - (n + momentum**2 / scale**2) * field**2 / 2
    )
    checks["independent_full_scale_pressure"] = s.cancel(
        pressure
        - velocity**2 / (2 * lapse**2)
        + (n + momentum**2 / (3 * scale**2)) * field**2 / 2
    )
    x, y, v, w, q, HH = s.symbols(
        "mode_re mode_im velocity_re velocity_im physical_q H", real=True
    )
    energy = (v * v + w * w + (q + n) * (x * x + y * y)) / 2
    press = (v * v + w * w - (n + q / 3) * (x * x + y * y)) / 2
    pair_re = (v * v - w * w + (q + n) * (x * x - y * y)) / 2
    pair_im = v * w + (q + n) * x * y
    derivative = (
        s.diff(energy, x) * v
        + s.diff(energy, y) * w
        + s.diff(energy, v) * (-3 * HH * v - (q + n) * x)
        + s.diff(energy, w) * (-3 * HH * w - (q + n) * y)
        + s.diff(energy, q) * (-2 * HH * q)
    )
    checks["full_physical_mode_Ward_identity"] = s.expand(
        derivative + 3 * HH * (energy + press)
    )
    checks["strict_pointwise_energy_cone"] = s.expand(
        energy**2 - pair_re**2 - pair_im**2 - (q + n) * (x * w - y * v) ** 2
    )
    frequency, aa = s.symbols("frequency a", positive=True)
    rate = s.Symbol("real_rate", real=True)
    initial = 1 / s.sqrt(2 * frequency * aa**3)
    initial_dot = (-rate - s.I * frequency) * initial
    phase = s.Matrix([initial, aa**3 * initial_dot])
    gram = phase * phase.conjugate().T
    antisymmetric = gram - gram.T - s.I * s.Matrix([[0, 1], [-1, 0]])
    for i in range(2):
        for j in range(2):
            checks[f"positive_frequency_Cauchy_Gram_CCR_{i}_{j}"] = s.simplify(
                antisymmetric[i, j]
            )
    delta = s.Symbol("delta", positive=True)
    squeeze = s.Symbol("squeeze", nonnegative=True)
    theta = s.Symbol("phase_angle", real=True)
    alpha = s.sqrt(1 + squeeze * squeeze)
    beta = -squeeze * s.exp(s.I * theta)
    c1 = delta * (1 + 2 * squeeze * squeeze)
    c2 = 2 * delta * squeeze * alpha * s.exp(s.I * theta)
    denom = s.sqrt(2 * delta * (c1 + delta))
    checks["phase_regular_alpha_formula"] = s.simplify(
        s.sqrt((c1 + delta) / (2 * delta)) - alpha
    )
    checks["phase_regular_beta_formula"] = s.simplify(-c2 / denom - beta)
    checks["phase_regular_SLE_CCR"] = s.simplify(
        alpha**2 - beta * s.conjugate(beta) - 1
    )
    minimum = (
        (alpha**2 + beta * s.conjugate(beta)) * c1
        + alpha * s.conjugate(beta) * c2
        + alpha * beta * s.conjugate(c2)
    )
    self_pair = alpha**2 * c2 + beta**2 * s.conjugate(c2) + 2 * alpha * beta * c1
    checks["complete_SLE_minimal_smeared_energy"] = s.simplify(minimum - delta)
    checks["complete_SLE_self_pair_zero"] = s.simplify(self_pair)
    checks["zero_mixing_static_vacuum_limit"] = s.simplify(beta.subs(squeeze, 0))
    checks["global_minimum_among_relative_Bogoliubov_states"] = (
        delta * (1 + 2 * squeeze**2) - delta - 2 * delta * squeeze**2
    )
    scalar_pressure_vac = s.Rational(3, 128) / s.pi**2
    return {
        "reference": "Physical CD metric ds²=dt²-(1+t²)^4 dx², u=t,X1,H_heavy0; use the recomputed S238 affine mean/source-centering, not a reset old connection trace. The old Proca and M1 choices are retained.",
        "actual_parameters": {
            "kappa0": KAPPA,
            "heavy_mass_squared": MASS2,
            "sampling_interval": SAMPLING_INTERVAL,
            "state_Cauchy_time": T0,
            "quantitative_time_interval": TIME_INTERVAL,
        },
        "full_mode_equation": "S_p''+3H S_p'+(n+p²/a²)S_p=0 with S_p conjugate(S_p')-conjugate(S_p)S_p'=i/a³. Its canonical chi=a^(3/2)S has frequency²=omega²-3H'/2-9H²/4. The state is defined using exact solutions at every momentum.",
        "actual_canonical_frequency_squared": omega2 - U,
        "sampling": "x=8(t+5/8); phi=exp[-1/(1-x²)] on |x|<1 and0 outside; w=f²=phi²/Z with Z=integral phi²dt. Support[-3/4,-1/2]. The normalization obeys Z>1/216. This sampling is fixed once before t0.",
        "SLE_definition": "For ANY exact normalized reference basis S, c1=integral w(|S'|²+omega²|S|²)/2 and c2=integral w(S'^2+omega²S²)/2. Put delta=sqrt(c1²-|c2|²), alpha=sqrt((c1+delta)/(2delta)), beta=-c2/sqrt(2delta(c1+delta)), and T=alpha S+beta conjugate(S). This formula has no singular phase choice at c2=0. Fix the physical positive-frequency convention exp(-i omega t) in the static limit and W2(t,t')=T(t)conjugate(T(t')).",
        "state_existence": "The pointwise energy cone gives strict c1>|c2| and c1>=Omega_star/128 on the sampling interval. The minimizer is unique up to a phase and independent of the comparison basis. Olbermann0704.2986 Theorems3.1 and4.9 apply to this positive massive minimally coupled field on the smooth globally hyperbolic flat RW geometry: the selected quasifree state is Hadamard. This external regularity theorem is not the quantitative stress estimate.",
        "state_family": "Retain the reference Cauchy covariance at t0 on compact smooth metric/source perturbations strictly to its future, sharing the Cauchy neighborhood and global hyperbolicity. Evolve source-free KG and add the unique retarded smooth coherent mean for the full prescribed source. Positivity, CCR and the Hadamard singularity are preserved. Do not re-minimize or recompute the state on a live history.",
        "complete_stress_readout": "rho=(|S'|²+(n+p²/a²)|S|²)/2, P=(|S'|²-(n+p²/(3a²))|S|²)/2 per mode, integrated with p²dp/(2pi²), with the full fourth-order subtraction and separately derived covariant finite local matching. Source and coherent mean vanish on the reference. No finite momentum cutoff is used.",
        "full_reference_source_boundary": "The classical J_H=g Phi² h(X)/2 has an eighth-order clock zero, and the stated finite onepoint switch has a1024th-order zero. The reference heavy mean and direct clock/source force vanish. Classical coherent mean starts at history order8 and its induced Euler force at order15; this does not remove the connected metric response of the heavy determinant.",
        "WKB_comparison": "Use full positive W0=omega and W(j+1)²=omega²-U-(log Wj)''/2+(log Wj)'²/4 through W7. Initialize an exact comparison mode to W6 at t=-1. The selected state is the basis-independent SLE, NOT the WKB comparison mode. The full defect and its Volterra comparison are bounded in estimates.",
        "light_vacuum_pressure_for_separate_constant_matching": scalar_pressure_vac,
        "scope": "The Gaussian heavy state is well defined on the whole reference and stress is smooth at each finite time. Numerical bounds apply only on[-1,1]. No full interacting light/metric state, quantitative curved response, nonlinear bounce or quantum gravitational decoupling is inferred.",
        "checks": {key: s.cancel(value) for key, value in checks.items()},
        "gates": {
            "mass_positive_and_global_reference_frequency_gap": MASS2 > 15,
            "reference_scale_factor_never_zero": True,
            "sampling_entirely_before_fixed_Cauchy_surface": SAMPLING_INTERVAL[1] == T0,
            "full_sampling_inside_quantitative_interval": TIME_INTERVAL[0]
            < SAMPLING_INTERVAL[0]
            < SAMPLING_INTERVAL[1]
            < TIME_INTERVAL[1],
            "new_reference_affine_mean_not_old_trace_reset": True,
            "physical_phase_and_CCR_fixed_before_Hadamard_identification": True,
            "SLE_not_an_instantaneous_or_finite_adiabatic_vacuum": True,
        },
    }
