"""Full dimensional energy/pressure finite parts and a common covariant action."""

from functools import cache

import sympy as s


@cache
def data():
    M, m, Q = s.symbols("M m Q", positive=True)
    Md, Mdd, H, Hd = s.symbols("Mdot Mddot H Hdot", real=True)
    eps = s.Symbol("epsilon")
    d = 3 - 2 * eps
    ell = s.Symbol("ell", real=True)
    K = Md + H * M
    A = K * K / 4 - s.Rational(3, 2) * H * M * K
    B = -M * (Mdd + (Hd - H * H) * M) / 2
    C = s.Rational(3, 2) * M * M * K * Md
    J0 = 4 / Q * (1 / eps - ell - s.Rational(2, 3))
    J1 = 4 / Q * (1 / eps - ell - s.Rational(16, 15))
    J2 = 8 / (5 * Q * M * M)
    rho_bare = K * K / Q * (1 / eps - ell - s.Rational(2, 3))
    P_bare = (A * J1 + B * J0 + C * J2) / d
    rho_ct = (
        -Md * Md - d * (d - 1) * M * M * H * H / 6 - s.Rational(2, 3) * d * H * M * Md
    ) / (Q * eps)
    P_ct = (
        -Md * Md / 3
        + s.Rational(2, 3) * M * Mdd
        + M * M * ((d - 1) * Hd / 3 + d * (d - 1) * H * H / 6)
        + s.Rational(2, 3) * (d - 1) * H * M * Md
    ) / (Q * eps)
    rho_MS = s.simplify(s.limit(rho_bare + rho_ct, eps, 0))
    P_MS = s.simplify(s.limit(P_bare + P_ct, eps, 0))
    rho_expected = (-K * K * ell - s.Rational(2, 3) * Md * Md + M * M * H * H) / Q
    Dfun = M * M * (ell - 1) + m * m
    Ddot = 2 * M * Md * ell
    Dddot = 2 * (Md * Md + M * Mdd) * ell + 4 * Md * Md
    rho_ref = rho_MS - m * m * H * H / Q
    P_ref = P_MS + m * m * (H * H + s.Rational(2, 3) * Hd) / Q
    P_cov = (
        -Md * Md * (ell + s.Rational(2, 3))
        + Dfun * (2 * Hd + 3 * H * H) / 3
        + Dddot / 3
        + 2 * H * Ddot / 3
    ) / Q
    rho_cov = (-Md * Md * (ell + s.Rational(2, 3)) - Dfun * H * H - H * Ddot) / Q
    checks = {
        "rho_dimensional_pole_cancellation": s.simplify(
            s.limit(eps * (rho_bare + rho_ct), eps, 0)
        ),
        "pressure_dimensional_pole_cancellation": s.simplify(
            s.limit(eps * (P_bare + P_ct), eps, 0)
        ),
        "full_rho_MS_finite_part": s.simplify(rho_MS - rho_expected),
        "same_covariant_rho_with_Newton_reference": s.simplify(rho_ref - rho_cov),
        "same_covariant_pressure_with_Newton_reference": s.simplify(P_ref - P_cov),
        "flat_energy_recovered": s.simplify(
            rho_MS.subs({H: 0, Hd: 0}) + Md * Md * (ell + s.Rational(2, 3)) / Q
        ),
        "flat_pressure_improvement_recovered": s.simplify(
            (P_MS - rho_MS).subs({H: 0, Hd: 0}) - Dddot / (3 * Q)
        ),
        "zero_field_constant_mass_rho_reference": s.simplify(
            rho_ref.subs({M: m, Md: 0, Mdd: 0, ell: 0})
        ),
        "zero_field_constant_mass_pressure_reference": s.simplify(
            P_ref.subs({M: m, Md: 0, Mdd: 0, ell: 0})
        ),
        "fixed_Newton_reference_contribution_to_rho": s.simplify(
            rho_ref - rho_MS + m * m * H * H / Q
        ),
        "fixed_Newton_reference_contribution_to_pressure": s.simplify(
            P_ref - P_MS - m * m * (H * H + s.Rational(2, 3) * Hd) / Q
        ),
    }
    return {
        "symbols": dict(
            zip(
                ("M", "Mdot", "Mddot", "H", "Hdot", "m", "Q", "ell"),
                (M, Md, Mdd, H, Hd, m, Q, ell),
            )
        ),
        "dimension": d,
        "second_order_mass_plus_geometry": K,
        "regulated_rho_second_order": rho_bare,
        "regulated_pressure_second_order": P_bare,
        "regulated_covariant_rho_counterterm": rho_ct,
        "regulated_covariant_pressure_counterterm": P_ct,
        "finite_MS_rho_second_order": rho_MS,
        "finite_MS_pressure_second_order": P_MS,
        "Newton_referenced_rho_second_order": rho_ref,
        "Newton_referenced_pressure_second_order": P_ref,
        "finite_curvature_function_D": Dfun,
        "finite_Ddot": Ddot,
        "finite_Dddot": Dddot,
        "common_finite_action": "L2=-(partial M)^2[log(M^2/mu^2)+2/3]/Q-M^2[log(M^2/mu^2)-1]R/(6Q). Add the explicitly fixed free-fermion Newton reference -m^2R/(6Q) per copy at mu=m.",
        "dimensional_rule": "Use the full d=3-2epsilon momentum measure and isotropic factor1/d, fixed four-component trace, and the d-dependent metric variation of the covariant counterterms before taking finite parts. The flat result alone did not fix these gravitational finite products.",
        "checks": checks,
    }
