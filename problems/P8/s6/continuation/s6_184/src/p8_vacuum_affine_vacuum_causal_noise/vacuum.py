"""Full fixed vacuum coefficient and conditional finite-time decoupling."""

from functools import cache

import sympy as s
from p8_vacuum_affine_quantum_retuning import profile

from . import source


@cache
def data():
    phi, Y = s.symbols("Phi Y", real=True)
    k0 = s.Symbol("kappa0", positive=True)
    n = source.N
    eps = profile.EPS
    X, u = profile.X, profile.u
    prior = profile.data()["full_coefficient_correction"]
    T = profile.original.data()["T"]
    V = -profile.P + profile.PV - (profile.rho + profile.P) * (X - 1) / 2
    changed = k0 * (prior + profile.PV)
    # Keep the symbolic anchor: expanding (10^800)^1024 would obscure
    # rather than improve the exact exponent comparison below.
    action = 8 * eps * 4**n / k0 ** (n - 1)
    scalar_terms = (
        8 * eps * 4**n / k0 ** (n - s.Rational(1, 2)),
        64 * n * eps * 4 ** (n - 1) / k0 ** (n - 1),
        16 * eps * 4**n / k0**n,
    )
    action_exp = -770 + 617 + 1 - 800 * (n - 1)
    force_exps = (
        -770 + 617 + 1 - 800 * (n - s.Rational(1, 2)),
        -770 + 617 + 5 - 800 * (n - 1),
        -770 + 617 + 2 - 800 * n,
    )
    worst = max(force_exps) + 1
    checks = {
        "literal_complete_nonconstant_fixed_correction": s.factor(changed - k0 * T * V),
        "literal_scalar_profile_X_derivative": s.diff(V, X)
        + (profile.rho + profile.P) / 2,
        "fixed_canonical_coefficient_independent_of_family_kappa": s.diff(
            changed, source.K
        ),
        "action_decimal_exponent": action_exp + 818552,
        "sum_force_decimal_exponent": worst + 818547,
        "full_profile_low_field_degree": 2 * n - 2048,
        "canonical_vacuum_constant_cancels_same_fixed_determinant": -profile.KAPPA
        * profile.PV
        + 5 * profile.MASS**4 / (128 * s.pi**2),
    }
    return {
        "new_fixed_canonical_nonconstant_coefficient": changed.subs(
            {u: phi / s.sqrt(k0), X: Y / k0}, simultaneous=True
        ),
        "constant_boundary": "The scalar constant -kappa0*p_v cancels the same-prescription flat Gaussian vacuum constant. It is retained before taking the gravitational limit, not silently dropped.",
        "full_scalar_action_difference_coefficient": action,
        "full_scalar_force_difference_terms": scalar_terms,
        "anchor_action_decimal_exponent_upper": action_exp,
        "anchor_force_decimal_exponent_upper": worst,
        "common_loose_decimal_upper": "1e-818500 (coefficient, not an enormous printed rational denominator)",
        "bounds_result": "After the explicit constant cancellation, the complete fixed canonical profile action is below1e-818500 integral U_Phi^2, and its scalar force pairing is below1e-818500 ||U_Phi||L2 ||U_eta||L2 on the stated real jet class.",
        "full_gaussian_scalar_boundary": "At fixed Minkowski g the same canonical three-polarization determinant has no Phi dependence. Its continuum scalar variations vanish after the fixed covariant prescription; the source mean/contact and connected linear force covariance above are retained exactly.",
        "actual_family_limit": "On every fixed admitted finite-time canonical class, the full source mean-force bound and noise variance decay as1/kappa, with fixed mass and full scalar functions. The fixed Deltaf itself persists in f+Deltaf; it is not declared to vanish with kappa.",
        "quantum_scope": "This is conditional Gaussian sector decoupling, not a regulator-removed interacting light/graviton/auxiliary measure limit or full vacuum amplitude theorem.",
        "checks": checks,
        "gates": {
            "reference_argument_in_numeric_stress_slab": s.Rational(1, 1)
            / s.sqrt(source.K0)
            < s.Rational(1, 2),
            "full_vacuum_real_X_interval_admitted": 4 / source.K0 < s.Rational(1, 4096),
            "Bernoulli_denominator_margin": 4 * (n + 1) / source.K0 < s.Rational(1, 2),
            "power_four_exponent_enclosure": 4**n < 10**617,
            "gradient_force_prefactor_enclosure": 64 * n < 10**5,
            "complete_action_below_loose_display": action_exp < -818500,
            "complete_force_below_loose_display": worst < -818500,
            "constant_matches_anchor_not_running_kappa": source.K0 == profile.KAPPA,
        },
    }
