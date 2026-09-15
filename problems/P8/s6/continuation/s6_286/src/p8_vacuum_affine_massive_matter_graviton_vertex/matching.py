"""Positive massive spectral measure, exact complex bounds and scoped pole matching."""

from functools import cache

import sympy as s

from . import source, vertex

MU, N, T, Z, V, G, K = (
    source.MU,
    source.N,
    source.T,
    source.Z,
    source.V,
    source.G,
    source.K,
)


def metric_contraction(energy, transfer, F1, F2, mass=MU):
    u = 4 * mass - energy - transfer
    h = mass - transfer / 4
    d = (energy - u) / 4
    return (
        4 * (d * d - h * h / 2) * F1 * F1
        + 2 * transfer * h * F1 * F2
        - s.Rational(3, 2) * transfer**2 * F2**2
    )


def bound_Pi_second():
    return source.CUBIC**2 / (432 * source.HEAVY_MASS2**2)


def bound_forward():
    return bound_Pi_second() / (3 * source.KAPPA)


@cache
def data():
    ss, aa, bb, f1, f2 = s.symbols("s channel_a channel_b F1 F2")
    nu = s.Symbol("crossing_v")
    Ntree = 2 * MU * MU - 2 * MU * T - ss * (4 * MU - ss - T)
    full = metric_contraction(ss, T, f1, f2)
    linear1 = s.diff(full, f1).subs({f1: 1, f2: -s.Rational(1, 2)})
    linear2 = s.diff(full, f2).subs({f1: 1, f2: -s.Rational(1, 2)})
    xi = s.Symbol("finite_Ricci_derivative_xi")
    y = s.Symbol("y", nonnegative=True)
    sigma = s.Symbol("spectral_s", positive=True)
    w = (1 - Z) ** 2 * (1 - V * V) / 4
    Al = vertex.denominator(MU, N, 0)
    Ah = vertex.denominator(N, MU, 0)
    x = s.Symbol("x", real=True)
    checks = {
        "whole_tree_metric_contraction": s.factor(
            full.subs({f1: 1, f2: -s.Rational(1, 2)}) - Ntree
        ),
        "whole_linear_F1_coefficient": s.factor(
            linear1 - (2 * Ntree + MU * T + T * T / 2)
        ),
        "whole_linear_F2_coefficient": s.factor(linear2 - T * (T + 2 * MU)),
        "entire_crossing_v_squared_projection_of_F1": s.expand(
            linear1.subs(ss, 2 * MU - T / 2 + nu)
        ).coeff(nu, 2)
        - 2,
        "entire_crossing_v_squared_projection_of_F2": s.expand(
            linear2.subs(ss, 2 * MU - T / 2 + nu)
        ).coeff(nu, 2),
        "literal_forward_Pi_second_normalization": -2 * s.Rational(1, 6) / K
        + 1 / (3 * K),
        "independent_Ricci_derivative_forward_shift": -2 * (xi / 2) / K + xi / K,
        "constant_improvement_full_crossing_sum": s.expand(
            (aa + 2 * MU) + (bb + 2 * MU) + (4 * MU - aa - bb + 2 * MU) - 10 * MU
        ),
        "light_triangle_whole_complex_gap": s.factor(
            Al - 4 * MU * w - (N * Z + MU * (1 - Z) ** 2 * V * V)
        ),
        "heavy_triangle_whole_complex_gap": s.factor(
            Ah - 4 * N * w - (N * Z * (1 - Z) + MU * Z * Z + N * (1 - Z) ** 2 * V * V)
        ),
        "first_sheet_light_spectral_quadratic": s.factor(
            (Al / (1 - Z) ** 2).subs(Z, y / (1 + y)) - (MU + N * y + N * y * y)
        ),
        "first_sheet_heavy_spectral_quadratic": s.factor(
            (Ah / (1 - Z) ** 2).subs(Z, y / (1 + y)) - (N + N * y + MU * y * y)
        ),
        "whole_spectral_cut_measure_Jacobian": s.factor(
            (Z * Z / (1 - Z)).subs(Z, y / (1 + y)) / (1 + y) ** 2 - y * y / (1 + y) ** 3
        ),
        "whole_light_geometric_remainder": s.factor(
            1 / (Al - w * T)
            - 1 / Al
            - T * w / Al**2
            - T * T * w * w / (Al * Al * (Al - w * T))
        ),
        "whole_heavy_geometric_remainder": s.factor(
            1 / (Ah - w * T)
            - 1 / Ah
            - T * w / Ah**2
            - T * T * w * w / (Ah * Ah * (Ah - w * T))
        ),
        "whole_light_ratio_control": s.factor(
            1 / (4 * MU) - w / Al - (N * Z + MU * (1 - Z) ** 2 * V * V) / (4 * MU * Al)
        ),
        "whole_heavy_ratio_control": s.factor(
            1 / (4 * N)
            - w / Ah
            - (N * Z * (1 - Z) + MU * Z * Z + N * (1 - Z) ** 2 * V * V) / (4 * N * Ah)
        ),
        "whole_scalar_Pi_second_pointwise_majorant": s.factor(
            (1 - Z) ** 2 / N**2
            - vertex.self_second_integrand()
            - MU * (1 - Z) ** 4 * (2 * N * Z + MU * (1 - Z) ** 2) / (N * N * Al * Al)
        ),
        "whole_scalar_Pi_second_integrated_majorant": s.integrate(
            (1 - Z) ** 2, (Z, 0, 1)
        )
        - s.Rational(1, 3),
        "equal_mass_independent_primary_benchmark": s.simplify(
            s.integrate(Z * Z * (1 - Z) ** 2 / (Z * Z - Z + 1) ** 2, (Z, 0, 1))
            - (s.Rational(5, 3) - 8 * s.sqrt(3) * s.pi / 27)
        ),
        "elementary_pi_lower_polynomial_remainder": s.factor(
            1 / (1 + x * x) - sum((-x * x) ** j for j in range(8)) - x**16 / (1 + x * x)
        ),
    }
    q = s.Symbol("geometric_ratio", positive=True)
    checks["all_positive_Taylor_geometric_remainder"] = s.factor(
        1 / (1 - q) - 1 - q - q * q / (1 - q)
    )
    checks["whole_spacelike_ratio_lower_margin"] = s.factor(
        1 / (1 + q) - 1 + q / (1 + q)
    )
    c = G * G / (16 * s.pi**2)
    A = s.Symbol("positive_Feynman_anchor", positive=True)
    W = s.Symbol("positive_Feynman_weight", positive=True)
    rho_location = A / W
    checks["once_subtracted_positive_spectral_pushforward"] = s.factor(
        (1 / W) * T / (rho_location * (rho_location - T)) - (1 / (A - W * T) - 1 / A)
    )
    density_light = (
        2
        * s.pi
        * c
        / sigma
        * y
        * y
        / (1 + y) ** 3
        / s.sqrt(1 - 4 * (MU + N * y + N * y * y) / sigma)
    )
    density_heavy = (
        2
        * s.pi
        * c
        / sigma
        * y
        * y
        / (1 + y) ** 3
        / s.sqrt(1 - 4 * (N + N * y + MU * y * y) / sigma)
    )
    pi_lower = 4 * sum(s.Rational((-1) ** j, 2 * j + 1) for j in range(8))
    bounds = {
        "Pi_second_upper": bound_Pi_second(),
        "F1prime_upper": bound_Pi_second() / 6,
        "forward_b2_magnitude_upper": bound_forward(),
        "whole_abs_t_le_2_b2_magnitude_upper": 2 * bound_forward(),
        "whole_abs_t_le_2_F1_correction_upper": 2 * bound_Pi_second() / 3,
    }
    return {
        "whole_metric_contraction": full,
        "whole_linear_F1_metric_coefficient": linear1,
        "whole_linear_F2_metric_coefficient": linear2,
        "whole_massive_F1_positive_spectral_density_integrands": (
            density_light,
            density_heavy,
        ),
        "whole_generated_forward_b2_piece": -G
        * G
        / (48 * s.pi**2 * K)
        * s.Integral(vertex.self_second_integrand(), (Z, 0, 1)),
        "actual_exact_arithmetic_upper_bounds": bounds,
        "positive_spectral_proof": "For n>=mu>0 the light and heavy first-sheet supports start at4mu and4n respectively. Set y=z/(1-z); the on-shell Feynman anchors divided by(1-z)^2 are mu+n*y+n*y^2 and n+n*y+mu*y^2. Their monotone positive quadratics exclude any earlier first-sheet singularity. For each support integrate the displayed positive density over0<=y<=ymax, with4(a+n*ymax+b*ymax^2)=sigma. The integrable endpoint square root is retained. This is the generated matter F1 measure, not the whole four-point absorptive density.",
        "complete_analytic_bound": "The entire OS-subtracted F1 is a once-subtracted positive Stieltjes integral. With r=abs(t)/(4mu)<1, abs(f1)<=abs(t)*f1prime(0)/(1-r), and abs(f1-t*f1prime(0))<=abs(t)^2*f1prime(0)/[4mu(1-r)]. Every higher coefficient is positive. The t0 continuation and limits therefore follow for these gapped graphs, without importing the massless-graviton form-factor limit.",
        "whole_spacelike_bound": "For t=-tau<=0, f1(-tau)/(-tau) is between f1prime(0)/(1+tau/(4mu)) and f1prime(0). It is positive, so the generated t-channel finite v^2 piece is negative and its magnitude is at most its forward magnitude. No absolute value sign is inferred for the full physical b20.",
        "t_channel_matching": "At linear matter-loop order the full t-channel metric contraction gives delta_b2_t(t)=-2*f1(t)/(kappa*t), analytically continued at0. Thus delta_b2_t(0)=-Pi_second(mu)/(3kappa). F2, including metric-H mixing and constant improvements, has no v^2 coefficient in this t-channel expression. The result is not the whole s/u crossed amplitude, improved low-energy subtraction or Regge remainder.",
        "actual_bound": "At mu1,n=10^200/512+2,g=1/8192,kappa=10^800,0<Pi_second<10^-405. The magnitude of the generated forward t-channel b2 piece and of its entire continuation on abs(t)<=2 are both below10^-1205. These are strict exact rational upper bounds using pi>3 and the full Feynman denominator, not a physical cutoff or an omitted-loop bound.",
        "primary_comparison": "Noumi and Tokuda, Phys.Rev.D104066022(2021), Eq12 and AppendixC, use c2=2b2 and R(-t)=F1(t)/2. Their cubic one-loop result agrees with the equal-mass benchmark after the identical-field symmetry factor is applied. Their symmetry expansion alone does not fix the sign; their Regge order estimate is not adopted as a P8 bound.",
        "checks": checks,
        "gates": {
            "strict_actual_positive_mass_hierarchy": source.HEAVY_MASS2 > 2,
            "elementary_pi_lower_bound_proved": pi_lower > 3,
            "strict_actual_Pi_second_upper": 0
            < bound_Pi_second()
            < s.Rational(1, 10**405),
            "strict_actual_whole_complex_t_channel_upper": 0
            < 2 * bound_forward()
            < s.Rational(1, 10**1205),
            "strict_actual_F1_whole_complex_disk_upper": 0
            < 2 * bound_Pi_second() / 3
            < s.Rational(1, 10**405),
            "whole_positive_measure_not_four_point_unitarity_claim": True,
            "full_crossed_amplitude_local_IR_and_Regge_obligations_retained": True,
        },
    }
