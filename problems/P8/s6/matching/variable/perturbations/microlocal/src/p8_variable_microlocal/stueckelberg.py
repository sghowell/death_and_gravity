"""Literal unreduced action and exact time boundaries; no H division.

All fields below are perturbative scalar coefficients. The transformation
is an invertible differential change of the already gauge-fixed eight
original fields; it is not an extra source, a new action or an assertion
that the five stationary leading fields are exact auxiliary solutions.
"""

from functools import cache

import sympy as sp
from p8_variable_constraints import action

PSIG, PSIF, PI, CHI, XG, XF, PHIG, PHIF = sp.symbols(
    "psi_g psi_f pi chi_B xi_g xi_f Phi_g Phi_f", real=True)
VPSIG, VPSIF, VPI, VCHI, VXG, VXF = sp.symbols(
    "psi_g_prime psi_f_prime pi_prime chi_B_prime xi_g_prime xi_f_prime", real=True)
FIELDS = (PSIG, PSIF, PI, CHI, XG, XF, PHIG, PHIF)
VELOCITIES = (VPSIG, VPSIF, VPI, VCHI, VXG, VXF, sp.S.Zero, sp.S.Zero)
PHYSICAL = (PI, XG, CHI)
PHYSICAL_VELOCITIES = (VPI, VXG, VCHI)
LEADING_AUXILIARIES = (PSIG, PSIF, PHIG, PHIF, XF)


def total_time(expression):
    """Time derivative of a field polynomial without acceleration slots."""
    u = action.derive()["u"]
    return sp.diff(expression, u)+sum(
        sp.diff(expression, q)*v for q, v in zip(FIELDS, VELOCITIES, strict=True))


def linear_velocity_boundary(lagrangian, indices):
    """Integrate a compatible linear-velocity row by polynomial identities.

No quadrature in background time is performed. This primitive is in
perturbation fields, whose coefficients remain time dependent.
"""
    boundary = sp.S.Zero
    for i in indices:
        target = sp.diff(lagrangian, VELOCITIES[i])-sp.diff(boundary, FIELDS[i])
        boundary += target*FIELDS[i]-sp.diff(target, FIELDS[i])*FIELDS[i]**2/2
    return boundary


@cache
def derive():
    d = action.derive()
    a, b, c, h, w = (d[key] for key in ("a", "b", "c", "h", "chi_speed"))
    K = action.K
    old_q = (PSIG-h*XG, PSIF+h*XF, -K*PI, CHI-w*XG)
    old_v = tuple(total_time(value) for value in old_q)
    substitution = dict(zip(d["coordinates"], old_q, strict=True))
    substitution.update(zip(action.lagrangian()["velocities"], old_v, strict=True))
    substitution.update({
        action.NG: PHIG-VXG, action.NF: PHIF-VXF,
        action.BG: XG/a**2, action.BF: VPI+c**2*XF/b**2,
    })
    literal = action.lagrangian()["L"].subs(substitution, simultaneous=True)
    L2 = sp.diff(literal, K, 2)/2
    L1 = sp.diff(literal, K).subs(K, 0)
    L0 = literal.subs(K, 0)
    # Remove all L1 field velocities except pi'. Compatibility is replayed
    # by checks, rather than assumed by the polynomial primitive.
    F1 = linear_velocity_boundary(L1, (0, 1, 3, 4, 5))
    reduced_L1 = L1-total_time(F1)
    # xi_f' is linear in L0 too. Remove it before its stationary value,
    # which contains pi', can introduce an artificial pi'' expression.
    F0 = linear_velocity_boundary(L0, (5,))
    reduced_L0 = L0-total_time(F0)
    return {**d, "fields": FIELDS, "velocities": VELOCITIES,
            "physical_fields": PHYSICAL, "physical_velocities": PHYSICAL_VELOCITIES,
            "leading_auxiliaries": LEADING_AUXILIARIES,
            "old_coordinates": old_q, "old_velocities": old_v,
            "old_auxiliary_substitution": {q: substitution[q] for q in action.lagrangian()["auxiliaries"]},
            "literal_L": literal, "K2_coefficient": L2,
            "L1": L1, "L0": L0, "F1": F1, "F0": F0,
            "boundary_subtracted_L1": reduced_L1,
            "boundary_subtracted_L0": reduced_L0,
            "boundary_convention": "literal L = K*L1bar + L0bar + d_u(K*F1+F0)"}


@cache
def checks():
    d = derive()
    residuals = {"literal_K2_cancellation": sp.factor(d["K2_coefficient"])}
    for i in (0, 1, 3, 4, 5):
        residuals[f"L1_no_velocity_{i}"] = sp.factor(sp.diff(d["boundary_subtracted_L1"], VELOCITIES[i]))
    residuals["L0_no_xif_velocity"] = sp.factor(sp.diff(d["boundary_subtracted_L0"], VXF))
    for name in ("F1", "F0"):
        for i, v in enumerate((*VELOCITIES[:6], PHIG, PHIF)):
            residuals[f"{name}_actual_configuration_boundary_{i}"] = sp.factor(sp.diff(d[name], v))
    return residuals
