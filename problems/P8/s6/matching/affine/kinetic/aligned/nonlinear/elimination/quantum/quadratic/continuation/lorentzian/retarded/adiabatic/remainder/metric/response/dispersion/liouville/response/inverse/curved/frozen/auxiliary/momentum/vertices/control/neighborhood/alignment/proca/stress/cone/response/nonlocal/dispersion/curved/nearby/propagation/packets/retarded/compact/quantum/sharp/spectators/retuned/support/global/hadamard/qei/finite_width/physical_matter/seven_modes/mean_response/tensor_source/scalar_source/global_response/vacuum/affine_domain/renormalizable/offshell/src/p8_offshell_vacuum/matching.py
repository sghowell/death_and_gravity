"""Actual polynomial parameters and exact untruncated heavy resolvent remainder."""

from functools import cache

import sympy as sp
from p8_affine_vacuum_domain import family, local
from p8_exceptional_vacuum import analytic, heavy
from p8_polynomial_vacuum import model

from . import jets


@cache
def data():
    d = model.data()
    lam, gamma = sp.symbols("positive_fixed_lambda positive_fixed_gamma", positive=True)
    D = d["D"]
    G2 = d["cubic_coupling_squared"]
    quartic = d["polynomial_quartic"]
    z, J = sp.symbols("centered_box_eigenvalue light_field_squared", real=True)
    c = G2 / D**2
    full = -quartic * J**2 / 24 + G2 * J**2 / (8 * (D + z))
    polynomial = (
        -quartic * J**2 / 24
        + G2 * J**2 * sum((-z) ** j / D ** (j + 1) for j in range(4)) / 8
    )
    separated = (
        c * J**2 / 12
        - c * J**2 * z / 8
        + lam * J**2 * z * z / 4
        - gamma * J**2 * z**3 / 8
    )
    rest = G2 * J**2 * z**4 / (8 * D**5 * (1 + z / D))
    parameters = {lam: analytic.VACUUM_LAMBDA_BAR, gamma: analytic.FIXED_GAMMA}
    radius = heavy.CHANNEL_RADIUS
    upper = (G2 * radius**4 / (8 * D**5 * (1 - radius / D))).subs(parameters)
    return {
        "D": D,
        "G_squared": G2,
        "polynomial_quartic": quartic,
        "redundant_c": sp.factor(c),
        "lambda": lam,
        "gamma": gamma,
        "scalar_resolvent_variable": "z is the eigenvalue of K=Box+2, not a scattering Mandelstam invariant",
        "literal_quartic_resolvent": full,
        "four_term_centered_truncation": polynomial,
        "separated_operator_coefficients": separated,
        "exact_operator_remainder": rest,
        "selected_spectral_radius": radius,
        "selected_quartic_action_remainder_upper": upper,
        "spectral_scope": "For J=Phi^2 in L2 spacetime with Fourier support |2-p^2|<=r<D, the real self-adjoint multiplier remainder obeys |Delta S4|<=upper ||J||_2^2. This is not a loop-momentum restriction, Euclidean derivative expansion or arbitrary cosmological history.",
        "checks": {
            "actual_constant_and_linear_centered_coefficients": sp.factor(
                polynomial - separated
            ),
            "exact_untruncated_resolvent_remainder": sp.factor(
                full - polynomial - rest
            ),
            "actual_quadratic_centered_coefficient": sp.factor(
                G2 / (8 * D**3) - lam / 4
            ),
            "actual_cubic_centered_coefficient": sp.factor(
                -G2 / (8 * D**4) + gamma / 8
            ),
            "actual_redundant_coefficient": sp.factor(c - 4 * lam**2 / gamma),
        },
        "bounds": {
            "selected_radius_inside_self_adjoint_resolvent_gap": radius
            < D.subs(parameters),
            "selected_radius_over_gap_below_one_e_minus_53": radius / D.subs(parameters)
            < sp.Rational(1, 10**53),
            "exact_action_remainder_positive_below_one_e_minus_400": 0
            < upper
            < sp.Rational(1, 10**400),
        },
    }


@cache
def actual_target():
    d = local.data()
    j = jets.data()
    m = data()
    lam, gamma = m["lambda"], m["gamma"]
    kap = sp.Symbol("positive_canonical_kappa", positive=True)
    phi, Y = sp.symbols("canonical_field canonical_gradient_square", real=True)
    old_lam = sp.Symbol("fixed_derivative_quartic", positive=True)
    substitute = {kap: family.N / gamma, old_lam: lam, phi: j["phi"], Y: j["X"]}
    lower = d["actual_quartic"].subs(substitute, simultaneous=True)
    # S6.109 proves RX(0)=0 and the removable Ia coefficients:
    # A3(0)=RXX(0), A4(0)=-RXX(0), A5(0)=0.
    rxx = (4 * d["source_c_over_X_at_vacuum"]).subs(family.u, 0)
    a3 = rxx * gamma / family.N
    full = sp.expand(lower + a3 * (j["L3"] - j["L4"]))
    jl, jg = sp.symbols("quartic_lambda quartic_gamma", real=True)
    target = j["target_quartic"].subs({jl: lam, jg: gamma}, simultaneous=True)
    return {
        "actual_new_affine_vacuum_lower_quartic": lower,
        "actual_canonical_A3_quartic": a3,
        "actual_full_quartic_target": full,
        "scope": "Only the actual canonical vacuum quartic jet in the fixed-interaction co-scaled limit. Not finite-M graviton matching or the entire analytic scalar action.",
        "checks": {
            "actual_canonical_quadratic_is_the_same_mass_one_free_action": sp.expand(
                d["actual_quadratic"].subs(substitute, simultaneous=True)
                - (j["X"] - j["phi"] ** 2) / 2
            ),
            "actual_scalar_lower_quartic": sp.expand(
                lower - lam * j["X"] ** 2 + gamma * j["phi"] ** 4 / 3
            ),
            "actual_new_removable_DHOST_coefficient": sp.factor(a3 + 2 * gamma),
            "full_quartic_target_comes_from_actual_affine_vacuum_jet": sp.expand(
                full - target
            ),
        },
    }
