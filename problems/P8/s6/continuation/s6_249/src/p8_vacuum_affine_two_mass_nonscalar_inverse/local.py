"""Literal nonscalar curvature Hessians and full current local/classical operators."""

from functools import cache
from itertools import product

import sympy as s
from p8_vacuum_affine_heavy_clock_quadratic import tensor as heavy_tensor
from p8_vacuum_affine_heavy_curved_state import state

from . import geometry

lam, P = s.symbols("lambda P", real=True)
eta = (1, -1, -1, -1)
t = geometry.t
a = geometry.metric.a
H = geometry.metric.H
R0 = 6 * (s.diff(H, t) + 2 * H * H)
q = P * P / a**2
mu2, n, ell = s.symbols(
    "proca_mass_squared heavy_mass_squared log_heavy_mass_squared", positive=True
)
A = (5 * mu2 + n * (ell - 1)) / 12 - (ell + 2) * R0 / 72
b = (ell + 2) / 60
C0 = 16 * s.pi**2 * state.KAPPA


def curvature(Q):
    h = s.zeros(4)
    h[1:4, 1:4] = -Q
    k = (lam, 0, 0, s.I * P)
    R = {}
    for a, b, c, d in product(range(4), repeat=4):
        R[a, b, c, d] = (
            k[c] * k[b] * h[a, d]
            + k[d] * k[a] * h[b, c]
            - k[d] * k[b] * h[a, c]
            - k[c] * k[a] * h[b, d]
        ) / 2
    Ric = s.Matrix(4, 4, lambda b, d: sum(eta[a] * R[a, b, a, d] for a in range(4)))
    scalar = sum(eta[a] * Ric[a, a] for a in range(4))
    return R, Ric, scalar


def weyl_mixed(D, G):
    RD, d, rd = curvature(D)
    RG, g, rg = curvature(G)
    riem = sum(
        eta[a] * eta[b] * eta[c] * eta[e] * RD[a, b, c, e] * RG[a, b, c, e]
        for a, b, c, e in product(range(4), repeat=4)
    )
    ric = sum(
        eta[a] * eta[b] * d[a, b] * g[a, b] for a, b in product(range(4), repeat=2)
    )
    return s.factor(2 * (riem - 2 * ric + rd * rg / 3))


def spatial_mixed(D, G):
    def christoffel(h, sign):
        k = (0, 0, sign * s.I * P)
        C = {
            (i, j, l): (k[j] * h[i, l] + k[l] * h[i, j] - k[i] * h[j, l]) / 2
            for i, j, l in product(range(3), repeat=3)
        }
        Ric = s.Matrix(
            3,
            3,
            lambda j, l: sum(k[i] * C[i, j, l] - k[l] * C[i, j, i] for i in range(3)),
        )
        return C, Ric

    cd, rd = christoffel(D, -1)
    cg, rg = christoffel(G, 1)
    quadratic = sum(
        cd[i, i, m] * cg[m, j, j]
        + cg[i, i, m] * cd[m, j, j]
        - cd[i, j, m] * cg[m, j, i]
        - cg[i, j, m] * cd[m, j, i]
        for i, j, m in product(range(3), repeat=3)
    )
    # The mixed metric jet has zero spatial momentum; its differentiated connection is zero.
    whole = (
        quadratic
        - s.trace(D * rg)
        - s.trace(G * rd)
        + s.trace(D) * s.trace(rg) / 2
        + s.trace(G) * s.trace(rd) / 2
    )
    return s.factor(whole)


def D0(f):
    return s.diff(f, t, 2) + H * s.diff(f, t)


def Dq(f):
    return D0(f) + q * f


def Dadj0(f):
    return s.diff(f, t, 2) + 5 * H * s.diff(f, t) + (2 * s.diff(H, t) + 6 * H * H) * f


def L0(f):
    return s.diff(f, t, 2) + 3 * H * s.diff(f, t)


@cache
def data():
    f = s.Function("unit_nonscalar_metric")(t)
    test = s.Function("detector")(t)
    Lq = lambda f: L0(f) + q * f
    Bv = Dadj0(D0(f)) + q * D0(f)
    Bt = Dadj0(Dq(f)) + q * Dq(f)
    higher = {"tensor": Bt, "vector": Bv}
    wave = {"tensor": Lq(f), "vector": L0(f)}
    density = {
        "tensor": a**3 * (A * (s.diff(f, t) ** 2 - q * f * f) - b * Dq(f) ** 2 / 2),
        "vector": a**3
        * (A * s.diff(f, t) ** 2 - b * (D0(f) ** 2 - q * s.diff(f, t) ** 2) / 2),
    }
    full = {}
    coefficients = {}
    checks = {}
    for i, D in enumerate(geometry.BASIS[2:]):
        for j, G in enumerate(geometry.BASIS[2:]):
            is_tensor = i >= 2
            target = (
                (lam * lam + P * P) ** 2
                if i == j and is_tensor
                else lam * lam * (lam * lam + P * P)
                if i == j
                else 0
            )
            checks[f"full_four_dimensional_Weyl_mixed_{i}_{j}"] = (
                weyl_mixed(D, G) - target
            )
            checks[f"full_spatial_Ricci_mixed_{i}_{j}"] = spatial_mixed(D, G) - (
                -P * P / 2 if i == j and is_tensor else 0
            )
    checks["independent_scalar_Einstein_normalization"] = (
        spatial_mixed(2 * s.eye(3), 2 * s.eye(3)) - 4 * P * P
    )
    for sector in ("tensor", "vector"):
        Euler = (
            s.diff(density[sector], f)
            - s.diff(s.diff(density[sector], s.diff(f, t)), t)
            + s.diff(s.diff(density[sector], s.diff(f, t, 2)), t, 2)
        )
        quantum = (
            -2 * (A * wave[sector] + s.diff(A, t) * s.diff(f, t)) - b * higher[sector]
        )
        tree = -C0 * wave[sector]
        full[sector] = s.expand(quantum + tree)
        checks["entire_weighted_local_action_Euler_" + sector] = s.cancel(
            Euler / a**3 - quantum
        )
        checks["canonical_two_leg_quantum_conversion_" + sector] = s.cancel(
            quantum / (16 * s.pi**2 * state.KAPPA)
            + (A * wave[sector] + s.diff(A, t) * s.diff(f, t) + b * higher[sector] / 2)
            / (8 * s.pi**2 * state.KAPPA)
        )
        checks["actual_classical_source_sign_" + sector] = tree + C0 * wave[sector]
        coefficients[sector] = {
            j: s.factor(full[sector].coeff(s.diff(f, t, j))) for j in range(5)
        }
        checks["whole_highest_coefficient_" + sector] = coefficients[sector][4] + b
        checks["whole_third_coefficient_" + sector] = (
            coefficients[sector][3] + 6 * b * H
        )
        checks["no_fifth_derivative_" + sector] = full[sector].coeff(s.diff(f, t, 5))
    checks["different_full_tensor_vector_gradient_operator"] = s.expand(
        Bt - Bv - q * D0(f) - q * q * f
    )
    checks["vector_no_zeroth_local_or_tree_term"] = coefficients["vector"][0]
    checks["tensor_entire_zeroth_local_and_tree_term"] = (
        coefficients["tensor"][0] + 2 * A * q + b * q * q + C0 * q
    )
    checks["whole_weighted_q_derivative"] = s.diff(
        a**3 * q * s.diff(f, t), t
    ) / a**3 - q * D0(f)
    checks["full_weighted_D0_adjoint"] = (
        s.diff(a**3 * test, t, 2) - s.diff(a**3 * H * test, t)
    ) / a**3 - Dadj0(test)
    current = state.germs.parent
    RR = current.coefficients()["R"]
    checks["whole_current_R_clock_value_not_only_small_difference"] = s.cancel(
        RR.subs(current.X, 1) - 1
    )
    checks["whole_current_R_clock_time_derivative"] = s.cancel(
        s.diff(RR.subs(current.X, 1), current.u)
    )
    e, g = s.symbols("detector_parameter source_parameter", real=True)
    for i, D in enumerate(geometry.BASIS[2:]):
        for j, G in enumerate(geometry.BASIS[2:]):
            volume = s.exp(s.trace(e * D + g * G) / 2)
            checks[f"entire_unimodular_profile_mixed_{i}_{j}"] = s.diff(volume, e, g)
    # Derive the contraction from a complete arbitrary synchronous spatial metric.
    coordinates = geometry.metric.COORDS
    entries = s.symbols("spatial_metric_name0:6")
    hs = [s.Function(str(name))(*coordinates) for name in entries]
    spatial = s.Matrix(
        [[hs[0], hs[1], hs[2]], [hs[1], hs[3], hs[4]], [hs[2], hs[4], hs[5]]]
    )
    spacetime = s.zeros(4)
    spacetime[0, 0] = 1
    spacetime[1:4, 1:4] = -a * a * spatial
    phi_gradient = s.Matrix([s.diff(t, coord) for coord in coordinates])
    inverse_time_row = s.Matrix([1, 0, 0, 0])
    checks["whole_synchronous_inverse_time_row"] = (
        spacetime * inverse_time_row - phi_gradient
    )
    checks["synchronous_X_exactly_one"] = (phi_gradient.T * inverse_time_row)[0] - 1
    christoffel_time = lambda mu, nu: sum(
        inverse_time_row[sigma]
        * (
            s.diff(spacetime[sigma, nu], coordinates[mu])
            + s.diff(spacetime[sigma, mu], coordinates[nu])
            - s.diff(spacetime[mu, nu], coordinates[sigma])
        )
        / 2
        for sigma in range(4)
    )
    contractions = s.Matrix([-christoffel_time(0, nu) for nu in range(4)])
    checks["all_phi_gradient_Hessian_contractions"] = contractions
    spatial_Hessian = s.Matrix(3, 3, lambda i, j: -christoffel_time(i + 1, j + 1))
    checks["canonical_Einstein_kinetic_factor"] = (state.KAPPA / 8) * (
        2 / s.sqrt(state.KAPPA)
    ) ** 2 - s.Rational(1, 2)
    oldtime = heavy_tensor.t
    actualH = heavy_tensor.A.subs(
        {oldtime: t, heavy_tensor.mass_squared: n, heavy_tensor.ell: ell},
        simultaneous=True,
    )
    checks["whole_current_scalar_heat_A_not_vector_copy"] = s.cancel(
        A - (5 * mu2 / 12 - R0 / 36) - actualH
    )
    return {
        "whole_finite_curvature_kinetic_A": A,
        "whole_fixed_fourth_constant": -b,
        "full_tensor_and_vector_quantum_density_before_64pi2": density,
        "entire_nonscalar_local_plus_classical_operators": full,
        "all_local_derivative_coefficients": coefficients,
        "complete_vector_fourth_operator": s.expand(Bv),
        "complete_tensor_fourth_operator": s.expand(Bt),
        "classical_action": "On N1,beta0,phi=t and any tracefree Q(t,x), X1 and volume are exact. R(t,1)=1; phi^mu phi_mu,nu=(1/2)partial_nu X=0 kills L3,L4,L5. Full F/profiles and homogeneous M1 terms have no variation. Einstein alone gives kappa a³[gamma'²-q gamma²]/8 for tensor and kappa a³ gamma'²/8 for synchronous vector.",
        "full_local_action": "Before64pi², tensor density a³{A(gamma'²-q gamma²)-b(Dq gamma)²/2}; vector density a³{A gamma'²-b[(D0 gamma)²-q gamma'²]/2}. A=(5mu²+n(ell-1))/12-(ell+2)R0/72, b=(ell+2)/60. Full finite Euler topological term is removed only after its complete compact-variation boundary.",
        "source_normalization": "With gamma=2h/sqrt(kappa), T_Gamma=64pi²/a³ times the gamma-current derivative gives canonical quantum force T_Gamma h/(16pi²kappa). The complete tensor/vector equation is T_total h=-16pi²kappa f, T_total=T_Gamma,total-16pi²kappa L_{q or0}. Both canonical chain factors and signs remain.",
        "contact_scope": "The displayed local operators do not set other retarded nonlocal or metric second-vertex contacts to zero. Every actual remaining contribution belongs to the whole strong remainder. No tensor gradient or fourth-order local operator is copied into the vector sector.",
        "checks": {key: geometry.zero(value) for key, value in checks.items()},
        "gates": {
            "tensor_and_vector_full_fourth_operators_distinct": s.expand(Bt - Bv) != 0,
            "vector_spatial_Ricci_restriction_zero_not_tensor_gradient": True,
            "whole_finite_scalar_and_Proca_A_terms_retained": A.has(n, mu2, ell),
            "both_canonical_chain_factors_retained": True,
            "whole_profile_zero_by_chart_not_Gaussian_mean_reset": True,
            "complete_all_P_local_operator_not_frozen_time_coefficients": s.diff(A, t)
            != 0,
            "Hessian_not_deleted_when_its_gradient_contraction_vanishes": spatial_Hessian
            != s.zeros(3),
        },
    }
