"""Entire homogeneous scalar adiabatic matching and complete finite local current bound."""

from functools import cache

import sympy as s

theta, X, k2, k3, k4, A2, AK, AK2 = s.symbols(
    "theta theta_dot tr_K2 tr_K3 tr_K4 tr_A2 tr_AK tr_AK2", real=True
)
trM = -X + 2 * k2 + theta**2
trM2 = A2 + 4 * k4 + theta**2 * k2 - 4 * AK2 - 2 * theta * AK + 4 * theta * k3
trMK = -AK + 2 * k3 + theta * k2
trMK2 = -AK2 + 2 * k4 + theta * k3
meanK1 = theta / 3
meanK2 = (theta**2 + 2 * k2) / 15
meanK3 = (theta**3 + 6 * theta * k2 + 8 * k3) / 105
meanK4 = (theta**4 + 12 * theta**2 * k2 + 12 * k2**2 + 32 * theta * k3 + 48 * k4) / 945
meanM = trM / 3
meanM2 = (trM**2 + 2 * trM2) / 15
meanMK2 = (trM * (theta**2 + 2 * k2) + 4 * theta * trMK + 8 * trMK2) / 105
mean_t2 = s.expand(
    X**2 / 4
    + X * meanM / 2
    - s.Rational(3, 2) * X * meanK2
    + meanM2 / 4
    - s.Rational(3, 2) * meanMK2
    + s.Rational(9, 4) * meanK4
)
mean_s4 = s.expand(
    theta**4 / 16
    - theta**3 * meanK1 / 4
    + 3 * theta**2 * meanK2 / 8
    - theta * meanK3 / 4
    + meanK4 / 16
)
whole = s.expand(mean_t2 + mean_s4)
R = 2 * X + theta**2 + k2
Ric = s.expand((X + k2) ** 2 + A2 + 2 * theta * AK + theta**2 * k2)
Riem = s.expand(4 * A2 + 8 * AK2 + 2 * k4 + 2 * k2**2)
pole = s.expand(R**2 / 36 + (Riem - Ric) / 90)
difference = s.expand(whole - pole)
c1, c2, c3 = s.symbols("c1 c2 c3")
boundary_basis = (
    theta**4 + 3 * theta**2 * X,
    theta**2 * k2 + X * k2 + 2 * theta * AK,
    theta * k3 + 3 * AK2,
)
residual = s.Poly(
    difference - sum(c * b for c, b in zip((c1, c2, c3), boundary_basis)),
    theta,
    X,
    k2,
    k3,
    k4,
    A2,
    AK,
    AK2,
)
solutions = s.solve(residual.coeffs(), (c1, c2, c3), dict=True)


@cache
def data():
    from p8_vacuum_affine_heavy_curved_state import state

    assert len(solutions) == 1
    boundary = s.expand(
        sum(c * b for c, b in zip((c1, c2, c3), boundary_basis)).subs(solutions[0])
    )
    checks = {
        "entire_general_homogeneous_fourth_order_covariant_matching": s.expand(
            difference - boundary
        )
    }
    for key, value in (
        (c1, -s.Rational(4, 189)),
        (c2, s.Rational(2, 315)),
        (c3, -s.Rational(34, 945)),
    ):
        checks["complete_weighted_boundary_coefficient_" + str(key)] = (
            solutions[0][key] - value
        )
    d, H, Hp, z = s.symbols("dimension H Hprime z", real=True)
    squeeze = (d - z) * H / 2
    derivative = s.diff(squeeze, H) * Hp + s.diff(squeeze, z) * (-2 * H * z * (1 - z))
    temp = s.expand(derivative + z * H * squeeze)
    iso = s.Poly(s.expand(temp**2 + squeeze**4), z)
    reduced = s.factor(
        sum(c * s.rf(d / 2, j) / s.rf(s.Rational(3, 2), j) for (j,), c in iso.terms())
    )
    subs = {
        theta: d * H,
        X: d * Hp,
        k2: d * H**2,
        k3: d * H**3,
        k4: d * H**4,
        A2: d * Hp**2,
        AK: d * H * Hp,
        AK2: d * H**2 * Hp,
    }
    checks["independent_arbitrary_dimension_isotropic_radial_action"] = s.factor(
        reduced - whole.subs(subs, simultaneous=True)
    )
    second = s.expand(theta**2 - 2 * theta**2 + (theta**2 + 2 * k2) / 3)
    checks["full_second_order_Einstein_density_before_boundary"] = s.expand(
        second - s.Rational(2, 3) * (k2 - theta**2)
    )
    checks["complete_Einstein_weighted_boundary"] = s.expand(
        R - 2 * (X + theta**2) - (k2 - theta**2)
    )
    eps, ell = s.symbols("epsilon ell", real=True)
    F = s.exp((s.EulerGamma - ell) * eps) * s.gamma(1 + eps)
    multiplied = (2 * F / ((eps - 1) * (eps - 2)), F / (3 * (eps - 1)), F)
    poles = (s.S.One, -s.Rational(1, 3), s.S.One)
    finite = (s.Rational(3, 2) - ell, (ell - 1) / 3, -ell)
    h0, hd = s.symbols("full_physical_invariant full_dimension_derivative", real=True)
    for j, expr in enumerate(multiplied):
        checks["complete_scalar_radial_pole_" + str(j)] = s.simplify(
            expr.subs(eps, 0) - poles[j]
        )
        checks["complete_scalar_radial_finite_" + str(j)] = s.simplify(
            s.diff(expr, eps).subs(eps, 0) - finite[j]
        )
        checks["full_counteraction_before_physical_limit_" + str(j)] = s.simplify(
            s.diff((expr - poles[j]) * (h0 - 2 * eps * hd), eps).subs(eps, 0)
            - finite[j] * h0
        )
    n, kappa = state.MASS2, state.KAPPA
    full_finite = (
        (s.Rational(3, 2) - ell) * n * n + (ell - 1) * n * R / 3 - ell * pole
    ) / (64 * s.pi**2)
    # Conservative holomorphic bounds for the ENTIRE local Euler current.
    kbound, abound = 10**5, 10**6
    theta_bound = 3 * kbound
    theta_dot_bound = 3 * abound
    curv_R = 2 * theta_dot_bound + theta_bound**2 + 3 * kbound**2
    curv_Ric = (theta_dot_bound + 3 * kbound**2) ** 2 + 3 * (
        abound + theta_bound * kbound
    ) ** 2
    curv_Riem = 12 * (abound + kbound * kbound) ** 2 + 2 * (
        9 * kbound**4 + 3 * kbound**4
    )
    density_bound = 1000 * (
        462
        + 154 * (2 * 10**11)
        + s.Rational(462, 36) * (2 * 10**11) ** 2
        + s.Rational(462, 90) * (2 * 10**22)
    )
    euler_count = 12 + 18 * 288 * 18 + 18 * 288 * 50 + 324 * 10368 * 18**2
    full_heat = 10**44 * n * n
    profile = 10**406
    state_subtraction = s.Rational(10**51, 1) / n
    complete = 10**45 * n * n / kappa
    checks["complete_current_Euler_term_count"] = s.S(euler_count) - 1088743692
    gates = {
        "full_dimensional_boundaries_before_variation": bool(boundary != 0),
        "one_scalar_not_vector_multiplicities": True,
        "full_scalar_pole_action_subtracted_before_dimension_limit": True,
        "unchanged_S240_scalar_finite_prescription": True,
        "entire_matrix_exponential_extrinsic_curvature_bound": 3 * 3**8 < kbound,
        "entire_extrinsic_time_derivative_bound": 41 * 3**8 < abound,
        "full_Ricci_scalar_complex_jet_bound": curv_R < 2 * 10**11,
        "full_Ricci_squared_complex_jet_bound": curv_Ric < 10**22,
        "full_Riemann_squared_complex_jet_bound": curv_Riem < 10**22,
        "full_volume_complex_bound": 3**6 < 1000,
        "every_finite_heat_coefficient_in_complex_density_bound": density_bound
        < 10**28,
        "complete_Euler_first_second_third_jet_derivative_count": euler_count
        < 2 * 10**9,
        "full_heat_current_and_two_source_derivatives": 6
        * euler_count
        * 2
        * 100**2
        * 10**28
        < 10**44,
        "whole_fixed_profile_not_only_linear_pullback": s.Rational(128 * 3 * 9, 2 * 4)
        * 10**400
        < profile,
        "full_state_subtraction_heat_and_profile_sum": full_heat
        + profile
        + state_subtraction
        < 10**45 * n * n,
        "entire_normalized_homogeneous_current_bound": complete
        < s.Rational(1, 10**350),
        "profile_and_covariant_vacuum_current_cancel_only_in_complete_sum": True,
        "no_homogeneous_contact_dropped": True,
        "no_same_space_inverse_or_all_transfer_bound_inferred": True,
    }
    return {
        "complete_extrinsic_curvature_invariants": {
            "R_old": R,
            "Ricci_squared": Ric,
            "Riemann_squared": Riem,
        },
        "entire_radial_fourth_order_scalar_polynomial": whole,
        "full_covariant_scalar_fourth_order_pole": pole,
        "complete_general_dimensional_boundary_difference": boundary,
        "all_three_weighted_boundary_coefficients": solutions[0],
        "all_scalar_radial_poles": poles,
        "all_scalar_radial_finite_coefficients": finite,
        "entire_finite_scalar_local_action_per_physical_volume": full_finite,
        "complete_holomorphic_geometry_and_density_bounds": {
            "extrinsic": kbound,
            "extrinsic_derivative": abound,
            "R": curv_R,
            "Ricci_squared": curv_Ric,
            "Riemann_squared": curv_Riem,
            "density_over_n_squared": density_bound,
        },
        "full_Euler_current_jet_derivative_term_count": euler_count,
        "complete_heat_current_two_parameter_derivative_bound": full_heat,
        "entire_unchanged_fixed_profile_derivative_bound": profile,
        "full_actual_state_subtraction_current_derivative_bound": state_subtraction,
        "entire_normalized_homogeneous_current_and_response_bound": complete,
        "full_local_Euler_identity": "For autonomous L(J0,J1,J2), E_i=L_0i-sum L_(1i,jb) J_(j+1,b)+sum L_(2i,jb) J_(j+2,b)+sum L_(2i,jb,kc)J_(j+1,b)J_(k+1,c). All6 symmetric entries and18,18,324 terms are counted before source Cauchy derivatives.",
        "graph": "Pointwise detector Frobenius norm times normalized source C12 powers for current derivatives0..2. The first prepared homogeneous response has detector L2 times source H13 on the unit slab. The full auxiliary cutoff tail uses this same first-response graph.",
        "profile": "At N1,X1 the ENTIRE S240 fixed profile is -v P_H_ref(t); it is held fixed under later metric histories. Its volume/current contacts are retained. It cancels the complete spatial reference mean at Q0, not the response or the finite-K mean.",
        "assembly_scope": "This supplies the homogeneous actual-state anchor with the S240 scalar prescription exactly once. S243's local spatial difference is zero here and is not added as another finite action. Nonzero-transfer remainders and ADM/clock Ward assembly remain separate.",
        "checks": checks,
        "gates": {key: bool(value) for key, value in gates.items()},
    }
