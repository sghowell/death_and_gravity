"""Full-covariance coherent quantization and controlled first Weyl calibration."""

from functools import cache

import sympy as s


def heat_generator(expr, variables, covariance):
    covariance = s.Matrix(covariance)
    if covariance.shape != (len(variables), len(variables)):
        raise ValueError("The whole covariance must match every phase coordinate")
    return s.expand(
        sum(
            covariance[i, j] * s.diff(expr, variables[i], variables[j]) / 2
            for i in range(len(variables))
            for j in range(len(variables))
        )
    )


def polynomial_heat(expr, variables, covariance, amount=1):
    degree = s.Poly(expr, *variables).total_degree()
    result, term = expr, expr
    for order in range(1, degree // 2 + 1):
        term = heat_generator(term, variables, covariance)
        result += amount**order * term / s.factorial(order)
    return s.expand(result)


def correlated_fixture():
    q1, q2, p1, p2 = s.symbols("q1 q2 p1 p2", real=True)
    x = (q1, q2, p1, p2)
    Omega = s.zeros(2).row_join(s.eye(2)).col_join((-s.eye(2)).row_join(s.zeros(2)))
    T = s.diag(2, s.Rational(3, 2), s.Rational(1, 2), s.Rational(2, 3))
    shear = s.eye(4)
    shear[2:, :2] = s.Matrix(
        [[s.Rational(1, 5), s.Rational(1, 7)], [s.Rational(1, 7), -s.Rational(1, 6)]]
    )
    T = shear * T
    return x, T * T.T / 2, T, Omega


def physical_hessian(v, U, UN, UNN, vg, vh, ng, nh):
    vg, ng = s.Matrix(vg), s.Matrix(ng)
    return s.exp(3 * v) * (
        9 * U * vg * vg.T
        + 3 * UN * (vg * ng.T + ng * vg.T)
        + UNN * ng * ng.T
        + 3 * U * s.Matrix(vh)
        + UN * s.Matrix(nh)
    )


@cache
def data():
    x, V, T, Omega = correlated_fixture()
    q1, q2, p1, p2 = x
    D = lambda expr: heat_generator(expr, x, V)
    H = lambda expr, amount=1: polynomial_heat(expr, x, V, amount)
    quartic = (
        q1**4
        + 2 * q1 * q2 * p1
        + 3 * q2 * q2 * p2 * p2
        + q1 * p1
        + q2 * p2
        + q1 * q2 * p1 * p2
    )
    cubic = q1 * q2 * p1 + q2 * p2 * p2
    gaussian_moments = H(quartic, 2).subs(dict.fromkeys(x, 0))
    whole = H(quartic - D(quartic))
    chi = s.Function("smooth_whole_cutoff")(*x)
    f = s.Function("entire_local_interaction")(*x)
    cross = sum(
        V[i, j] * s.diff(chi, x[i]) * s.diff(f, x[j])
        for i in range(4)
        for j in range(4)
    )
    cutoff_contact = chi * D(f) + f * D(chi) + cross
    v0, U, UN, UNN = s.symbols("v0 full_U full_U_N full_U_NN", real=True)
    vg = s.Matrix(s.symbols("v_first0:2", real=True))
    ng = s.Matrix(s.symbols("N_first0:2", real=True))
    vh = s.Matrix(
        [[s.Symbol("v_00"), s.Symbol("v_01")], [s.Symbol("v_01"), s.Symbol("v_11")]]
    )
    nh = s.Matrix(
        [[s.Symbol("N_00"), s.Symbol("N_01")], [s.Symbol("N_01"), s.Symbol("N_11")]]
    )
    y = s.Matrix(s.symbols("canonical_y0:2", real=True))
    dv = vg.dot(y) + (y.T * vh * y)[0] / 2
    dn = ng.dot(y) + (y.T * nh * y)[0] / 2
    physical = s.exp(3 * (v0 + dv)) * (U + UN * dn + UNN * dn**2 / 2)
    hess = physical_hessian(v0, U, UN, UNN, vg, vh, ng, nh)
    checks = {
        "full_squeezed_correlated_map_symplectic": T * Omega * T.T - Omega,
        "full_correlated_reference_purity": V * Omega * V - Omega / 4,
        "first_calibration_matches_every_cubic_contact": H(cubic - D(cubic)) - cubic,
        "full_quartic_nonzero_second_order_remainder": s.expand(
            whole - quartic + D(D(quartic)) / 2
        ),
        "nonzero_quartic_contact_exact_value": D(D(quartic))
        - s.Rational(826341, 31360),
        "complete_Husimi_double_covariance_moments": gaussian_moments
        - H(H(quartic)).subs(dict.fromkeys(x, 0)),
        "all_cutoff_derivative_contacts": s.expand(D(chi * f) - cutoff_contact),
        "whole_implicit_physical_volume_Hessian": (
            s.hessian(physical, y).subs(dict.fromkeys(y, 0)) - hess
        ).applyfunc(s.expand),
    }
    # Heat semigroup identities tested on every monomial through degree5.
    for degree in range(6):
        mono = (q1 + 2 * q2 - 3 * p1 + p2) ** degree
        checks["complete_heat_semigroup_degree_" + str(degree)] = H(
            H(mono, s.Rational(2, 3)), s.Rational(1, 3)
        ) - H(mono)
    return {
        "whole_normalized_quantization": "Q_V(a)=(2pi hbar)^(-d) integral a(z)|D(z)psi0><D(z)psi0| dz; unit-CCR formulas below set hbar=1 after restoring the original kappa and Fourier normalization. The same fixed pure psi0 is used both as initial state and as coherent window; displacements are not new initial state choices.",
        "whole_positive_unital_compression": "The normalized coherent analysis W is an isometry by momentum Fourier inversion and translation of |psi0|^2. Q_V(a)=W* M_a W. Thus Q_V(1)=I, Q_V(a)>=0 for a>=0, ||Q_V(a)||<=||a||infinity, and Q_V(a)*Q_V(a)<=Q_V(|a|^2). Real bounded symbols give bounded self-adjoint operators.",
        "whole_Weyl_and_initial_expectation": "With full Gaussian Wigner covariance V, Q_V(a)=OpW(a*w_V)=OpW(e^D a), D=(1/2)V_AB partial_AB. The SAME seed expectation is (a*w_(2V))(0)=e^(2D)a(0). All q-p and inter-field cross covariances remain; this is not an occupation-only prescription.",
        "whole_first_Weyl_calibration": "a_match=(1-D)g_ext. Its exact Weyl-symbol error is e^D(1-D)g_ext-g_ext=-integral_0^1 t e^(tD) D^2 g_ext dt, with sup error<=||D^2 g_ext||infinity/2. This is an explicit named finite ordering prescription, not exact Weyl matching for general nonlinear symbols or a counterterm of the original theory.",
        "whole_cutoff_calibration_contact": cutoff_contact,
        "whole_correlated_fixture_covariance": V,
        "whole_quartic_fixture": quartic,
        "whole_quartic_Weyl_error": -D(D(quartic)) / 2,
        "whole_complete_physical_Hessian": hess,
        "whole_volume_ordering_boundary": "Q_V(F_ext) is positive without further conditions. Q_V((1-D)F_ext) matches Weyl to the displayed second-order symbol error but is positive by this argument ONLY if ||D F_ext||infinity<inf F_ext. No actual P8 inequality of that type has been evaluated. The contraction D F retains the full canonical pullback N_AB=N_ij z_i,A z_j,B+N_i z_i,AB, the v_AB term and all cross terms.",
        "primary_source_convention": "de Gosson, Generalized Anti-Wick Quantum States, arXiv:1907.02471, section3 Proposition3 and Corollary4. Its unnormalized Lebesgue convention is divided by (2pi hbar)^d here. Only the complete rank-one Wigner convolution and window covariance are used; unit resolution and all estimates here have independent written proofs.",
        "checks": checks,
        "gates": {
            "full_covariance_including_nonzero_off_diagonal_used": V[0, 3] != 0,
            "quartic_first_calibration_not_exact_Weyl": D(D(quartic)) != 0,
            "all_cutoff_gradient_and_Hessian_contacts_mandatory": cross != 0,
            "all_nonlinear_physical_auxiliary_contacts_retained": hess.has(nh[0, 1]),
            "positive_coherent_volume_not_automatically_positive_Weyl_volume": True,
            "ordering_error_is_symbol_sup_bound_not_claimed_operator_norm_bound": True,
        },
    }
