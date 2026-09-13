"""Exact finite-regulator heavy integration and full one-light-loop generator."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_scalar_tree_matching import model

NAME = model.NAME
m2, g, C = s.symbols("heavy_mass_squared trilinear contact", real=True)
phi0, phi1, H0, H1 = s.symbols("phi0 phi1 H0 H1", real=True)


@cache
def data():
    phi, H = s.Matrix([phi0, phi1]), s.Matrix([H0, H1])
    KH = s.Matrix([[m2 + 1, -1], [-1, m2 + 1]])
    GH = KH.inv()
    KL = s.Matrix([[2, -1], [-1, 2]])
    f = phi.applyfunc(lambda v: v * v)
    action = (
        (phi.T * KL * phi)[0] / 2
        + (H.T * KH * H)[0] / 2
        - g * (H.T * f)[0] / 2
        - C * sum(v**4 for v in phi) / 24
    )
    center = g * GH * f / 2
    effective = (
        (phi.T * KL * phi)[0] / 2
        - C * sum(v**4 for v in phi) / 24
        - g * g * (f.T * GH * f)[0] / 8
    )
    DP = s.diag(*phi)
    W = -C * s.diag(*f) / 2 - g * g * s.diag(*(GH * f)) / 2 - g * g * DP * GH * DP
    A = s.hessian(action, (phi0, phi1)).subs({H0: center[0], H1: center[1]})
    B = s.Matrix([[s.diff(action, p, h) for h in (H0, H1)] for p in (phi0, phi1)])
    I0, j, z = s.symbols("coincident_light_integral source_counterterm marker")
    checks = {
        "complete_heavy_square": s.cancel(
            action - effective - ((H - center).T * KH * (H - center))[0] / 2
        ),
        "full_stationary_heavy_action": s.cancel(
            action.subs({H0: center[0], H1: center[1]}) - effective
        ),
        "unchanged_positive_quartic_coefficient": s.cancel(
            -model.C / 24 - model.g**2 / (8 * (model.D + 2)) - model.VALLEY_QUARTIC
        ),
        "heavy_constant_mode_inverse": s.cancel(sum(GH.row(0)) - 1 / m2),
        "one_loop_heavy_tadpole_counterterm_cancellation": -g * g * I0 / (2 * m2)
        + g * (g * I0 / 2) / m2,
        "source_square_quadratic_cross": s.expand(
            -((g * z * z / 2 - j) ** 2) / (2 * m2)
            + g * g * z**4 / (8 * m2)
            + j * j / (2 * m2)
            - g * j * z * z / (2 * m2)
        ),
        "full_log_quadratic_marker_coefficient": s.diff(s.log(1 + z) / 2, z).subs(z, 0)
        - s.Rational(1, 2),
        "full_log_quartic_marker_coefficient": s.diff(s.log(1 + z) / 2, z, 2).subs(z, 0)
        / 2
        + s.Rational(1, 4),
        "mixed_two_point_wick_factor": s.Rational(2 * 2 * 2, 2**2 * 2) - 1,
        "heavy_two_point_identical_wick_factor": s.Rational(2 * 2, 2**2 * 2)
        - s.Rational(1, 2),
        "light_contact_tadpole_wick_factor": s.Rational(4 * 3, s.factorial(4))
        - s.Rational(1, 2),
        "heavy_one_point_wick_factor": s.Rational(1, 2) - s.Rational(1, 2),
    }
    for label, matrix in {
        "full_light_hessian": s.hessian(effective, (phi0, phi1)) - KL - W,
        "ordered_full_two_field_Schur_complement": A - B * GH * B.T - KL - W,
        "heavy_gradient_equation": KH * center - g * f / 2,
    }.items():
        for i, value in enumerate(matrix):
            checks[label + "_" + str(i)] = s.cancel(value)
    return {
        "model": NAME,
        "exact_Euclidean_light_action": "S_eff=1/2<phi,K_phi phi>-C integral phi^4/24-g^2<phi^2,G_H phi^2>/8; K_H=-Delta+M_H^2, G_H=K_H^-1. The Gaussian normalization is (det K_H)^(-1/2).",
        "regulated_scope": "An exact identity at a specified finite Euclidean regulator and finite volume with positive K_H>=M_H^2 and canonical light operator. The same quadratic-form inequality holds on its continuum test-function domain. No continuum interacting quantum measure, reflection-positive continuum limit, S matrix or high-energy arc is constructed.",
        "positive_nonlocal_action_bound": "Since0<G_H<=1/M_H^2, S_eff>=1/2<phi,K_phi phi>+q integral phi^4, q=g^2(D-1)/[6D^2(D+2)]>0 for the S233 bare classical parameters. This finite-regulator statement is not a claim that all renormalized bare counterterms preserve this bound uniformly as a regulator is removed.",
        "full_light_hessian_operator": "K_phi+W[phi], W=-C M_(phi^2)/2-g^2 M_(G_H phi^2)/2-g^2 M_phi G_H M_phi. M denotes multiplication. Both local and bilocal terms are retained.",
        "complete_formal_one_light_loop": "Gamma_1LPI=S_eff+hbar/2 Tr log(K_phi+W)+counterterms+O(hbar^2); its field-dependent one-loop terms through phi^4 are hbar/2 Tr(G_phi W)-hbar/4 Tr(G_phi W G_phi W). This is the complete one-light-particle-irreducible loop generator, not an evaluated renormalized four-point matching amplitude.",
        "mixed_loop_warning": "The flat heavy determinant has no light dependence, but the remaining light trace contains all mixed heavy-light loops. In a curved metric the heavy determinant is metric-dependent and cannot be dropped from stress or gravity.",
        "one_point_condition": "At first loop order set the heavy one-point function to zero with Euclidean linear counterterm j1 H, j1=g I_phi/2, where I_phi is the regulated coincident light propagator. Completing the full source square gives a quadratic term g j1 phi^2/(2M_H^2), cancelling the local heavy-source tadpole term in Tr(G_phi W)/2. The quartic tadpole and mixed bubble remain.",
        "boundary_prescription": "The Euclidean identity analytically continues to the Feynman/in-out observable with the full heavy propagator. Replacing it by a retarded kernel inside a one-copy action is not this identity; causal expectation values require their own in-in construction.",
        "finite_matrix_fixture": {
            "K_H": KH,
            "K_phi": KL,
            "W": W,
            "effective_action": effective,
        },
        "checks": checks,
        "gates": {
            "source_pinned_separate_model": NAME == "V2S-T1",
            "actual_nonlocal_quartic_lower_bound_positive": model.VALLEY_QUARTIC.subs(
                {model.D: model.GAP, model.g: model.G}
            )
            > 0,
            "both_hessian_terms_and_heavy_source_tadpole_retained": True,
            "complete_one_light_loop_not_heavy_determinant_only": True,
            "finite_regulator_identity_not_continuum_UV_construction": True,
        },
    }
