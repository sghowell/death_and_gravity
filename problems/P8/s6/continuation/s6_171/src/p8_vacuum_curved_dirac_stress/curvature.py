"""Explicit evanescent curvature prescription and finite Euler stress."""

from functools import cache

import sympy as s


@cache
def data():
    d = s.Symbol("spatial_dimension", positive=True)
    eps = s.Symbol("epsilon")
    H, Hd = s.symbols("H Hdot", real=True)
    Q = s.Symbol("Q", positive=True)
    R = -2 * d * Hd - d * (d + 1) * H * H
    Ric = d * d * (Hd + H * H) ** 2 + d * (Hd + d * H * H) ** 2
    Riem = 4 * d * (Hd + H * H) ** 2 + 2 * d * (d - 1) * H**4
    Euler = s.expand(Riem - 4 * Ric + R * R)
    Weyl = s.factor(Riem - 4 * Ric / (d - 1) + 2 * R * R / (d * (d - 1)))
    C = d * (d - 1) * (d - 2)
    expected = C * (4 * H * H * Hd + (d + 1) * H**4)
    lapse = s.Symbol("lapse", positive=True)
    A = s.Symbol("a", positive=True)
    v = s.Symbol("coordinate_H", real=True)
    # Local lapse action after one integration by parts; H=v/lapse.
    reduced = -C * (d - 3) * A**d * v**4 / (3 * lapse**3)
    rho_E = s.simplify(-s.diff(reduced, lapse).subs(lapse, 1) / A**d).subs(v, H)
    P_E = C * (d - 3) * (H**4 + 4 * H * H * Hd / d)
    anomaly_a = s.Rational(11, 360)
    anomaly_c = s.Rational(1, 20)
    rho_CT = s.simplify(
        s.limit(anomaly_a * rho_E.subs(d, 3 - 2 * eps) / (2 * Q * eps), eps, 0)
    )
    P_CT = s.simplify(
        s.limit(anomaly_a * P_E.subs(d, 3 - 2 * eps) / (2 * Q * eps), eps, 0)
    )
    x, y, z = s.symbols("R_squared Ricci_squared Riemann_squared")
    dim = d + 1
    Epoly = z - 4 * y + x
    Cpoly = z - 4 * y / (dim - 2) + 2 * x / ((dim - 1) * (dim - 2))
    literal = (5 * x - 8 * y - 7 * z) / 360
    difference = s.factor(literal - (anomaly_a * Epoly - anomaly_c * Cpoly))
    expected_difference = (dim - 4) * (
        y / (10 * (dim - 2)) - (dim + 1) * x / (60 * (dim - 1) * (dim - 2))
    )
    finite_basis = s.simplify(
        s.limit(difference.subs(d, 3 - 2 * eps) / (2 * eps * Q), eps, 0)
    )
    checks = {
        "complete_dimensional_Euler_invariant": s.factor(Euler - expected),
        "Weyl_tensor_square_vanishes_in_every_dimension_on_FLRW": Weyl,
        "Euler_integration_by_parts": s.factor(
            expected
            - s.Rational(4, 3) * C * (d * H**4 + 3 * H * H * Hd)
            + C * (d - 3) * H**4 / 3
        ),
        "lapse_variation_before_dimension_limit": s.factor(rho_E + C * (d - 3) * H**4),
        "dimensional_Euler_conservation": s.factor(
            s.diff(rho_E, H) * Hd + d * H * (rho_E + P_E)
        ),
        "four_dimensional_Euler_variation_zero": rho_E.subs(d, 3),
        "finite_evanescent_Euler_energy": s.factor(rho_CT - 11 * H**4 / (60 * Q)),
        "finite_evanescent_Euler_pressure": s.factor(
            P_CT + 11 * (H**4 + s.Rational(4, 3) * H * H * Hd) / (60 * Q)
        ),
        "finite_Euler_conservation": s.factor(
            s.diff(rho_CT, H) * Hd + 3 * H * (rho_CT + P_CT)
        ),
        "finite_Euler_trace": s.factor(
            rho_CT - 3 * P_CT - anomaly_a * Euler.subs(d, 3) / Q
        ),
        "literal_and_EC_pole_basis_difference": s.factor(
            difference - expected_difference
        ),
        "finite_pole_basis_difference": s.factor(
            finite_basis - (-y / s.Integer(20) + x / s.Integer(72)) / Q
        ),
        "finite_basis_in_four_dimensional_C_E_R_basis": s.factor(
            finite_basis
            - (-(z - 2 * y + x / 3) / 40 + (z - 4 * y + x) / 40 - x / 360) / Q
        ),
        "same_rank_four_heat_kernel_mass_quartic_normalization": s.Rational(1, 2)
        * 4
        * s.Rational(1, 2)
        - 1,
        "spin_curvature_gamma_trace_coefficient": s.Rational(4, 180)
        - s.Rational(1, 24)
        + s.Rational(7, 360),
        "heat_kernel_R_squared_coefficient": 4
        * (s.Rational(1, 72) + s.Rational(1, 32) - s.Rational(1, 24))
        - s.Rational(5, 360),
    }
    return {
        "physical_scalar_curvature_in_D": R,
        "Ricci_squared_in_D": Ric,
        "Riemann_squared_in_D": Riem,
        "Euler_invariant_in_D": Euler,
        "reduced_lapse_action_Euler": reduced,
        "Euler_energy_before_D_limit": rho_E,
        "Euler_pressure_before_D_limit": P_E,
        "finite_Euler_energy_per_copy": rho_CT,
        "finite_Euler_pressure_per_copy": P_CT,
        "literal_heat_kernel_curvature_polynomial": literal,
        "literal_minus_EC_pole_basis": difference,
        "finite_literal_minus_EC_counterterm": finite_basis,
        "selected_prescription": "GY14-SAT8-MR/EC-N0 uses fixed four-component trace and curvature poles [11E_D/360-C_D^2/20]/(2Qepsilon), with no additional finite R^2 term in THIS E_D,C_D basis. Total box R terms are boundary terms under compact metric variations. The evanescent continuation is explicit; flat data did not fix it.",
        "finite_basis_warning": "The literal rank-four heat-kernel polynomial and this E_D,C_D basis have the same four-dimensional pole residue but differ by a finite local term. On FLRW their stress difference is the variation of -R^2/(360Q); they are not silently identified.",
        "source_context": "The rank-four heat-kernel calculation uses the local coefficient and spin-curvature traces in arxiv1703.00908v2 equations106,109,110. Euclidean continuation is calibrated to the inherited +M^4/(Qepsilon) counterterm; R_E continues to -R_P8. The finite Euler stress sign is derived by lapse variation here, not imported from an unconverted trace-anomaly convention.",
        "checks": checks,
    }
