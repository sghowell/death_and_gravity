"""Full flat potential and exact prescribed-source stress, with explicit scope."""

from functools import cache

import sympy as s
from p8_exceptional_vacuum import analytic

from . import family


@cache
def data():
    u = family.u
    original = family.original
    k, n, g, C, Gam = family.K0, family.MASS2, family.G, family.CONTACT, family.GAMMA
    D = family.model.GAP
    q = family.model.VALLEY_QUARTIC.subs({family.model.D: D, family.model.g: g})
    tree = original.original.data()["original_retuned_tree_scalar"].subs(family.X, 0)
    poly = (
        8800 * u**14
        + 53600 * u**12
        + 136800 * u**10
        + 180050 * u**8
        + 123200 * u**6
        + 38099 * u**4
        + 4573 * u**2
        + 224
    )
    tv = -poly / (200 * (1 + u * u) ** 8)
    fv = -u * u / 2 + (-family.N / s.Integer(3) + s.Rational(28, 25)) * u**4
    actual_fv = analytic.vacuum_lower_function().subs(family.X, 0)
    blend = original.data()["F"].subs(family.X, 0)
    phi, H, V0 = s.symbols("Phi H V_old_relative", real=True)
    ng, gg, cg, ga = s.symbols("n g C gamma", real=True)
    full = (
        V0
        - ga * phi**4 / 3
        - cg * phi**4 / 24
        + ng * H * H / 2
        - gg * H * phi * phi / 2
    )
    beta = -Gam / 3 - C / 24 - g * g / (8 * (n - s.Rational(1, 3)))
    generic_beta = -ga / 3 - cg / 24 - gg * gg / (8 * (ng - s.Rational(1, 3)))
    square = (
        V0
        + H * H / 6
        + (ng - s.Rational(1, 3))
        * (H - gg * phi**2 / (2 * (ng - s.Rational(1, 3)))) ** 2
        / 2
        + generic_beta * phi**4
    )
    # Metric variation of the complete source term is retained in the stress.
    Y, h, hy, t = s.symbols("positive_Y h h_Y t", real=True)
    J = gg * Y * t * t * h / 2
    lag = (
        s.symbols("Hdot", real=True) ** 2 / 2 - ng * H * H / 2 + H * gg * phi**2 * h / 2
    )
    d = gg * Y * h / ng**2
    hpart = J / ng - d
    hvalley = J / ng
    velocity = s.diff(hpart, t)
    JY = gg * phi**2 * hy / 2
    energy = lambda field: (
        velocity**2 / 2 + ng * field**2 / 2 + field * (2 * Y * JY - J)
    )
    pressure = lambda field: velocity**2 / 2 - ng * field**2 / 2 + field * J
    de = s.cancel(energy(hpart) - energy(hvalley))
    dp = s.cancel(pressure(hpart) - pressure(hvalley))
    eta2 = g * g * k / (s.sqrt(family.LOCALIZER) * n**3)
    amplitude_bound = s.Integer(10) ** 190
    low_metric_error = 2 * g * g * amplitude_bound**2 / n**2
    exp_bound = s.factorial(8) * 2**8 / family.LOCALIZER**8
    high_metric_error = 4 * g * g * k * family.LOCALIZER * exp_bound / n**2
    metric_error = max(low_metric_error, high_metric_error)
    unrestricted_counterexample_lower = g * g * k / (9 * n**2)
    checks = {
        "literal_full_original_tree_at_zero_gradient": s.cancel(tree - tv),
        "literal_vacuum_lower_at_zero_gradient": s.cancel(actual_fv - fv),
        "literal_full_original_potential_blend": s.cancel(
            blend - (s.exp(-(u**4)) * fv + (1 - s.exp(-(u**4))) * tv)
        ),
        "complete_new_coercive_square": s.cancel(full - square),
        "heavy_square_valley_quartic": s.cancel(-C / 24 - g * g / (8 * n) - q),
        "new_classical_quadratic_light_mass": s.diff(-k * blend, u, 2).subs(u, 0) / k
        - 1,
        "full_forced_heavy_equation": s.diff(hpart, t, 2) + ng * hpart - J,
        "source_second_box_zero": s.diff(J, t, 4),
        "particular_and_valley_same_velocity": s.diff(hpart - hvalley, t),
        "full_source_metric_energy_difference": s.cancel(
            de
            - (
                gg * gg * Y**2 * h * h / (2 * ng**3)
                - gg * gg * Y**2 * phi**2 * h * hy / ng**2
            )
        ),
        "full_pressure_difference": s.cancel(dp + gg * gg * Y**2 * h * h / (2 * ng**3)),
        "complete_source_lagrangian_kept": s.diff(lag, H, 2) + ng,
    }
    return {
        "literal_original_tree_at_Y_zero": tv,
        "literal_original_vacuum_lower_at_Y_zero": fv,
        "original_potential_relative_to_its_retained_origin": "-kappa0 F_base(Phi/sqrt(kappa0),0); it is nonnegative and has its only zero at Phi=0 by the two complete negative terms and the positive blend.",
        "new_coercive_square": square,
        "actual_quartic_q": q,
        "actual_coercive_quartic_beta": beta,
        "global_classical_potential": "Vnew-Vnew(0)>=H^2/6+q Phi^4/2 on every real constant Phi,H with Y=0. The literal original potential is retained. This proves a coercive unique classical origin, not a quantum vacuum or healthy arbitrary-gradient theory.",
        "prescribed_background": "Flat metric, constant 0<=Y<=6kappa0/5 and Phi=sqrt(Y)t, restricted pointwise to either |Phi|<=10^190 at any X=Y/kappa0, or 1/2<=X<=6/5 and |Phi|<=sqrt(kappa0). This connected coefficient-domain corridor contains the vacuum anchor and clock tube in the nonnegative-X sector; no open spacelike-gradient neighborhood is inferred from this homogeneous stress bound. This is not an interpolating solution or a coupled scalar trajectory. The exact particular heavy solution is chosen with no additional homogeneous heavy wave.",
        "exact_source": J,
        "exact_particular_heavy": hpart,
        "retained_valley_heavy": hvalley,
        "full_source_metric_energy_difference": de,
        "full_pressure_difference": dp,
        "source_response_proxy_squared_upper": eta2,
        "low_gradient_chart_field_amplitude_bound": amplitude_bound,
        "source_metric_stress_low_chart_upper": low_metric_error,
        "source_metric_stress_high_chart_upper": high_metric_error,
        "source_metric_stress_corridor_upper": metric_error,
        "full_error_bound": "On the stated two-chart corridor only, |rho_part-rho_valley|/Y<=eta^2/2+max(2g^2*10^380/n^2,4g^2*kappa0*A*8!*2^8/(n^2*A^8))<10^-8 and |p_part-p_valley|/Y<=eta^2/2<10^-8, with continuous zero-gradient interpretation. The source metric derivative is included. The full, possibly large valley kinetic energy is retained.",
        "unrestricted_stress_counterexample_lower": unrestricted_counterexample_lower,
        "unrestricted_stress_counterexample": "At X=A^-1/2 and Phi=sqrt(kappa0), the negative h_Y term makes the relative energy error positive and greater than g^2*kappa0/(9n^2)>10^395, even though the heavy offset energy proxy is tiny. The chosen linear path reaches this point at a correspondingly large time. There is no unrestricted-strip small heavy-stress theorem.",
        "not_inferred": "No coupled bounce, global mass gap, heavy-state transport, loop stress, cutoff, nonlinear stability, exact UV or original V/G/B closure follows from this prescribed-source identity.",
        "checks": {key: s.cancel(value) for key, value in checks.items()},
        "gates": {
            "all_original_negative_polynomial_coefficients_positive": all(
                c > 0 for c in s.Poly(poly, u).all_coeffs() if c != 0
            ),
            "vacuum_lower_quartic_negative": -family.N / s.Integer(3)
            + s.Rational(28, 25)
            < 0,
            "actual_quartic_after_original_potential_correction_positive": q > Gam / 3,
            "actual_coercive_heavy_square_positive": n > s.Rational(1, 3),
            "actual_coercive_quartic_exceeds_half_q": beta > q / 2,
            "actual_source_response_proxy_squared_below_1e_minus_8": eta2
            < s.Rational(1, 10**8),
            "actual_low_chart_source_metric_error_below_1e_minus_21": low_metric_error
            < s.Rational(1, 10**21),
            "actual_high_chart_source_metric_error_below_1e_minus_2500": high_metric_error
            < s.Rational(1, 10**2500),
            "global_switch_derivative_majorant": 8 / s.sqrt(family.LOCALIZER) + 1 < 2,
            "high_chart_switch_derivative_majorant": s.Rational(48, 5)
            + s.Rational(72, 25) * family.LOCALIZER
            < 4 * family.LOCALIZER,
            "counterexample_switch_factor_exceeds_one_half": (
                1 - 1 / s.sqrt(family.LOCALIZER)
            )
            ** 16
            > s.Rational(1, 2),
            "unrestricted_source_metric_error_exceeds_1e395": unrestricted_counterexample_lower
            > 10**395,
            "actual_complete_prescribed_stress_error_below_1e_minus_8": eta2 / 2
            + metric_error
            < s.Rational(1, 10**8),
        },
    }
