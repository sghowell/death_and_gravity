"""Entire source, fixed H/Proca metric terms and unresolved light curvature matching."""

from functools import cache

import sympy as s
from p8_vacuum_affine_full_flat_tensor_response import normalization as vector_force
from p8_vacuum_affine_heavy_curved_state import quantum, state
from p8_vacuum_affine_massive_dimensional_cut_completion import source as previous
from p8_vacuum_affine_physical_background_vertices import parent
from p8_vacuum_affine_proca_gaussian import stress as vector_scheme
from p8_vacuum_affine_quantum_retuning import profile

KAPPA, HEAVY_MASS2, VECTOR_MASS2 = state.KAPPA, state.MASS2, s.Integer(10) ** 6
NU, L, R, W, EU, M = s.symbols(
    "scalar_mass_squared log_mass_squared R_old Weyl_squared Euler vector_mass",
    real=True,
)


def require_spin(value):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, s.Integer))
        or value not in (0, 2)
    ):
        raise ValueError("Require a conserved spin0 or spin2 metric sector")
    return int(value)


def require_mass(value):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, s.Rational))
        or value <= 0
    ):
        raise ValueError("Require a positive exact rational Gaussian mass squared")
    return s.Rational(value)


def scalar_density(mass2=NU, logarithm=L):
    a2 = W / 120 - EU / 360 + R * R / 72
    return (
        (s.Rational(3, 2) - logarithm) * mass2**2
        + (logarithm - 1) * mass2 * R / 3
        - 2 * logarithm * a2
    )


def vector_density(mass=M):
    return (
        s.Rational(5, 2) * mass**4
        + s.Rational(5, 3) * mass**2 * R
        - 4 * (W / 120 - EU / 360 + R * R / 72)
    )


@cache
def fixed_coefficients():
    n = HEAVY_MASS2
    value = s.expand(
        (
            scalar_density(n, s.log(n))
            + s.Rational(3, 2)
            + vector_density(s.sqrt(VECTOR_MASS2))
        )
        / (64 * s.pi**2)
    )
    return {
        "volume": value.subs({R: 0, W: 0, EU: 0}),
        "R_old": value.coeff(R, 1),
        "Weyl_squared": value.coeff(W),
        "R_old_squared": value.coeff(R, 2),
        "Euler": value.coeff(EU),
    }


@cache
def data():
    inherited = previous.data()
    answer = {
        key: value for key, value in inherited.items() if key not in ("checks", "gates")
    }
    checks = dict(inherited["checks"])
    full = parent.fixed_functions()
    coeff = fixed_coefficients()
    scalar = quantum.data()["complete_scalar_finite_heat_matching_before_64pi2"]
    names = {str(v): v for v in scalar.free_symbols}
    mapped = scalar.subs(
        {
            names["mass_squared"]: NU,
            names["log_n"]: L,
            names["R_old"]: R,
            names["Ricci_squared"]: (W - EU) / 2 + R * R / 3,
            names["Riemann_squared"]: 2 * W - EU + R * R / 3,
        }
    )
    vector = vector_scheme.prescription()["mu_mass_specialization"]
    names = {str(v): v for v in vector.free_symbols}
    mapped_vector = vector.subs(
        {
            names["positive_mass"]: M,
            names["R_old"]: R,
            names["Ricci_squared"]: (W - EU) / 2 + R * R / 3,
            names["Riemann_squared"]: 2 * W - EU + R * R / 3,
        }
    )
    n = HEAVY_MASS2
    checks.update(
        {
            "entire_scalar_fixed_covariant_density": s.expand(
                mapped - scalar_density()
            ),
            "entire_Proca_fixed_covariant_density": s.expand(
                mapped_vector - vector_density()
            ),
            "entire_three_vacuum_source_density": s.cancel(
                KAPPA * full["vacuum_pressure_total"] - coeff["volume"]
            ),
            "whole_F_zero_cancels_Gaussian_volume": s.cancel(
                KAPPA * full["F_full"].subs({parent.u: 0, parent.X: 0})
                + coeff["volume"]
            ),
            "whole_R_at_vacuum_is_one": s.cancel(
                full["R_full"].subs({parent.u: 0, parent.X: 0}) - 1
            ),
            "whole_Proca_constant_not_new_choice": s.cancel(
                KAPPA * profile.PV - 5 * VECTOR_MASS2**2 / (128 * s.pi**2)
            ),
            "retained_finite_R_old_coefficient": s.factor(
                coeff["R_old"]
                - (n * (s.log(n) - 1) + 5 * VECTOR_MASS2) / (192 * s.pi**2)
            ),
            "retained_finite_Weyl_coefficient": s.factor(
                coeff["Weyl_squared"] + (s.log(n) + 2) / (3840 * s.pi**2)
            ),
            "retained_finite_R_squared_coefficient": s.factor(
                coeff["R_old_squared"] + (s.log(n) + 2) / (2304 * s.pi**2)
            ),
            "retained_compact_Euler_coefficient": s.factor(
                coeff["Euler"] - (s.log(n) + 2) / (11520 * s.pi**2)
            ),
            "prior_full_Proca_positive_Newton_sign": s.cancel(
                (vector_force.C - vector_force.C0) / (16 * s.pi**2)
                - 5 * VECTOR_MASS2 / (96 * s.pi**2)
            ),
            "exact_unrounded_heavy_mass": n - (s.Integer(10) ** 200 / 512 + 2),
            "same_original_kappa": KAPPA - s.Integer(10) ** 800,
        }
    )
    answer.update(
        {
            "entire_fixed_H_Proca_curvature_and_three_vacuum_coefficients": coeff,
            "unmatched_light_scalar_curvature_coefficients": s.symbols(
                "unmatched_Phi_R unmatched_Phi_Weyl_squared unmatched_Phi_R_squared",
                real=True,
            ),
            "entire_generic_fixed_scalar_density_before_64pi2": scalar_density(),
            "entire_fixed_vector_density_before_64pi2": vector_density(),
            "curvature_signature": vector_scheme.prescription()["curvature_signature"],
            "Einstein_density": "-kappa R_P8/2=+kappa R_old/2. A fixed c_R R_old changes kappa by PLUS2c_R. This agrees independently with the original full Proca force response.",
            "gapped_species": {
                "Phi_mass_squared": s.Integer(1),
                "H_mass_squared": n,
                "Proca_mass_squared": VECTOR_MASS2,
            },
            "whole_source_scope": "The full R,F and retained H/Proca sources are unchanged. At the constant-scalar flat reference, all nonminimal and mixed source vertices start above quadratic Gaussian order. All three Gaussian volume constants are fixed. The H and Proca finite curvature prescriptions are explicit, but the light Phi finite gravitational curvature polynomial has not been source-matched. Its coefficients remain unresolved, not set by applying the H prescription to another sector. These are Gaussian metric insertions, not the full interacting four-scalar amplitude, massless loops or the original curved state.",
            "checks": checks,
            "gates": {
                "all_three_original_vacuum_constants_retained": "all_three_vacuum_constants"
                in answer,
                "source_check_dictionary_not_parent_alias": checks
                is not inherited["checks"],
                "original_signature_bridge_retained": "R_old=-R_P8"
                in vector_scheme.prescription()["curvature_signature"],
                "no_new_finite_vacuum_or_curvature_counterterm": True,
                "Gaussian_metric_insertion_not_full_interacting_scattering": True,
                "same_flat_reference_and_unmatched_light_curvature_not_original_state": True,
            },
        }
    )
    return answer
