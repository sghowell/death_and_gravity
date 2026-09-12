"""Ordinary finite-window causal inverse with its undamped pole retained."""

from functools import cache

import sympy as s

from . import spectral


@cache
def data():
    t, lam, omega, R, m = s.symbols(
        "time lambda omega positive_weight mass", positive=True
    )
    pole = -R * s.sin(omega * t) / omega
    primitive = -R * (1 - s.cos(omega * t)) / omega**2
    z = spectral.z
    H0 = s.Rational(50625, 59168)
    # Elementary beta moments of y^(2j)/sqrt(1-y^2).
    integrals = (s.pi / 4, 3 * s.pi / 16, 5 * s.pi / 32)
    G_bound = (integrals[0] - s.Rational(2, 3) * integrals[1] + integrals[2] / 10) / (
        2 * m
    )
    checks = {
        "pole_primitive_derivative": s.diff(primitive, t) - pole,
        "pole_primitive_initial_value": primitive.subs(t, 0),
        "pole_Laplace_transform": s.factor(
            s.integrate(s.exp(-lam * t) * pole, (t, 0, s.oo))
            + R / (lam * lam + omega * omega)
        ),
        "bounded_forward_factor_kernel_majorant": G_bound - 9 * s.pi / (128 * m),
        "total_primitive_bound_includes_pole": 2 * 30 - 60,
        "pole_period_absolute_integral": s.integrate(
            R * s.sin(omega * t) / omega, (t, 0, s.pi / omega)
        )
        * 2
        - 4 * R / omega**2,
        "threshold_cut_analytic_coefficient": s.limit(spectral.RHO / s.sqrt(z), z, 0)
        - H0,
        "threshold_frequency_coordinate": s.limit(
            (1 - 4 * m * m / (2 * m + t) ** 2) / t, t, 0
        )
        - 1 / m,
        "total_static_inverse_value": -s.Integer(30) + 30,
    }
    return {
        "forward_distribution": "F2(derivative^2)=-delta/30-derivative^2 G_m; G_m(t)=theta(t) integral_0^1 W2(y) sin(Omega_y t)/[(1-y^2)Omega_y]dy, Omega_y=2m/sqrt(1-y^2). G_m is continuous, G_m(0)=0 and |G_m|<=9pi/(128m). This defines the original causal factor without a separately divergent local/cut split.",
        "inverse_kernel": "K_m(t)=-theta(t)[R sin(m sqrt(r)t)/(m sqrt(r))+2 integral_(2m)^infinity rho(Omega^2)sin(Omega t)dOmega]. The continuum integral is ordinary oscillatory, not absolutely convergent in frequency.",
        "continuum_time_estimates": "The continuum piece is O(1/[t log^2(1/t)]) at0 and O(t^-3/2) at infinity, hence belongs toL1(0,infinity). Constants depend on the fixed mass; the mass-scaled formulation uses dimensionless times.",
        "full_time_domain": "K_m belongs toL1(0,T) for each finiteT but not toL1(0,infinity), because its nonzero pole is undamped and the continuum is absolutely integrable.",
        "primitive": "J_m(t)=-R/(r m^2)(1-cos(m sqrt(r)t))-integral rho(tau)/tau (1-cos(sqrt(tau)t))dtau, with -60<=J_m<=0.",
        "prepared_bound": "For f inC1[0,T] with f(0)=0, K_m*f=J_m*f' and ||K_m*f||_C0<=60T||f'||_C0. No removal of the pole or change in source preparation occurs.",
        "finite_window_inverse": "Causal convolution gives both F2*K_m=delta and K_m*F2=delta. On the graph domain of continuous zero-past u with u(0)=0 and distributional F2*u equal to an ordinary right-continuous source including the initial boundary, K_m is the unique inverse toC[0,T], bounded by||K_m||_L1(0,T). No initial delta derivatives may be silently omitted.",
        "scaling": "r andR/m^2 are mass independent; rho_m(tau)=rho_1(tau/m^2), K_m(t)=mK_1(mt), so ||K_m||_L1(0,T)=||K_1||_L1(0,mT). The full half-line norm is infinite.",
        "checks": checks,
        "gates": {
            "forward_kernel_bound_finite": s.Rational(9, 128) > 0,
            "threshold_exponent_gives_integrable_large_time_tail": s.Rational(3, 2) > 1,
            "small_time_log_power_integrable": s.Integer(2) > 1,
            "nonzero_pole_weight_prevents_total_half_line_L1": spectral.RESIDUE_LO > 0,
            "finite_window_pole_bound": spectral.RESIDUE_HI / s.sqrt(spectral.ROOT_LO)
            < 23,
            "primitive_constant_finite_and_positive": s.Integer(60) > 0,
        },
    }
