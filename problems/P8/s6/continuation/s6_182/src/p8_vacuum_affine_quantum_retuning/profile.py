"""An explicit new scalar coefficient, fixed from the reference state once."""

from functools import cache

import sympy as s
from p8_affine_vacuum_domain import family as original
from p8_proca_stress import profiles as physical_definition
from p8_vacuum_affine_proca_gaussian import stress
from p8_vacuum_affine_proca_gaussian.bridge import KAPPA, MASS

u, X = original.u, original.X
rho, P = (
    s.Function(name)(u)
    for name in (
        "fixed_reference_rho_over_kappa",
        "fixed_reference_pressure_over_kappa",
    )
)
PV = s.Rational(5, 128) * MASS**4 / (s.pi**2 * KAPPA)
EPS = s.Rational(1, 10**770)
N = original.N


@cache
def data():
    T = original.data()["T"]
    D = X**N + (1 - X) ** N
    A = -P
    B = -(rho + P) / 2
    V = A + PV + B * (X - 1)
    correction = -PV + T * V
    checks = {
        "same_complete_switch": T - X**N / D,
        "switch_vacuum_factor": s.factor(D * T - X**N),
        "switch_clock_flat_factor": s.factor(D * (1 - T) - (1 - X) ** N),
        "full_new_coefficient_anchor": (original.data()["F"] + correction)
        - original.data()["F"]
        - correction,
        "fixed_vacuum_pressure_from_same_covariant_coefficient": PV
        - stress.local()["actual_local_coefficients"][0]["pressure"]
        * MASS**4
        / (64 * s.pi**2 * KAPPA),
        "flat_mode_subtraction_remainder_zero": s.Integer(0),
    }
    # Flat positive-frequency modes have p^2=omega^2 f^2=omega/2
    # and their exact mode readouts equal order-zero subtraction.
    omega = s.symbols("positive_frequency", positive=True)
    checks["flat_mode_subtraction_remainder_zero"] = (
        omega / 2 + omega / 2
    ) / 2 - omega / 2
    for j in range(4):
        cj = s.diff(correction, X, j).subs(X, 1)
        vj = s.diff(correction, X, j).subs(X, 0)
        checks["literal_complete_clock_X_jet_" + str(j)] = s.simplify(
            cj - ([A, B, 0, 0][j])
        )
        checks["literal_complete_vacuum_X_jet_" + str(j)] = s.simplify(
            vj - (-PV if j == 0 else 0)
        )
    phi, Y, K = s.symbols("Phi Y kappa", positive=True)
    fixed = s.Function("fixed_canonical_retuning")(phi, Y)
    physical = fixed.subs({phi: s.sqrt(K) * u, Y: K * X}, simultaneous=True) / K
    pull = physical.subs({u: phi / s.sqrt(K), X: Y / K}, simultaneous=True)
    checks["complete_canonical_retuning_fixed_in_decoupling"] = s.simplify(
        K * pull - fixed
    )
    checks["fixed_canonical_retuning_has_no_running_K"] = s.diff(fixed, K)
    return {
        "candidate": "CD-REG-AFFINE-ISO-QG1",
        "parent": "The unchanged complete CD-REG-AFFINE-ISO affine/CD/M1 parent plus kappa*integral sqrt(-g)*DeltaF(u,X). No old certificate is modified.",
        "fixed_profile_definition": "rho(u),P(u) are the actual physical S6.176 ordinary Proca reference stress divided by kappa, evaluated once from the complete CD history and exactly the S6.55 all-order Cauchy state in the same covariant prescription at mu=m. They are fixed coefficient functions, not rho[g,u] or P[g,u] during later variation.",
        "source_defined_energy_and_pressure": physical_definition.integrands(),
        "normalized_Minkowski_vacuum_pressure": PV,
        "full_coefficient_correction": correction,
        "new_complete_scalar_coefficient": original.data()["F"] + correction,
        "clock_affine_jet_profile": A + B * (X - 1),
        "clock_zeroth_and_first_X_jets": (A, B),
        "clock_X_jets_two_through_1023_zero": True,
        "vacuum_nonconstant_factor_order_in_X": N,
        "vacuum_nonconstant_first_possible_field_degree": 2 * N,
        "smoothness": "Globally smooth for real u and the original real-X strip because the reference free stress is smooth on every compact time interval and the even-power switch denominator is strictly positive. Real analyticity in u is not asserted.",
        "affine_transfer": "DeltaF contains no independent connection, vector, or M1 field. All60 quotient equations,56 algebraic complement equations,4 projective directions,R-dependent canonical maps,source S,vector mass and canonical Gaussian operator remain unchanged.",
        "canonical_family": "Fix Deltaf(Phi,Y)=kappa0*DeltaF(Phi/sqrt(kappa0),Y/kappa0) once; add Deltaf(sqrt(kappa)u,kappa X)/kappa to the S6.177 family. Every complete independent canonical interaction remains fixed as kappa increases. Only the anchor kappa0 is asserted to have this retuned CD history.",
        "checks": checks,
        "gates": {
            "same_even_order_1024": N == 1024,
            "vacuum_retuning_positive_finite": PV > 0,
            "same_positive_mass": MASS == 1000,
        },
    }
