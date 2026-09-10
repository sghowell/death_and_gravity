"""Leading heavy-fermion gauge threshold with its field-domain remainder."""

from functools import cache

import sympy as sp


@cache
def data():
    phi = sp.symbols("scalar_Phi", real=True)
    m, y, a = sp.symbols(
        "positive_fermion_mass positive_Yukawa positive_gauge_squared", positive=True
    )
    z = sp.symbols("nonnegative_field_ratio", nonnegative=True)
    mass = sp.diag(m + y * phi, m - y * phi)
    determinant = sp.factor((mass.T * mass).det())
    derivative = sp.diff(sp.log(determinant), phi, 2).subs(phi, 0)
    C = sp.factor(a * derivative / (192 * sp.pi**2))
    reference = 2 * sp.log(1 - y * y * phi * phi / (m * m))
    odd = sp.diff(reference, phi).subs(phi, 0)
    beta_threshold = sp.Rational(2, 3)
    running_reference = (
        a * beta_threshold / (32 * sp.pi**2) * sp.log(1 - y * y * phi * phi / (m * m))
    )
    coefficients = {str(n): -sp.Rational(2, n) for n in range(1, 9)}
    checks = {
        "direct_determinant": determinant - (m * m - y * y * phi * phi) ** 2,
        "differentiated_log_determinant": derivative + 4 * y * y / (m * m),
        "opposite_pair_odd_coupling_cancels": odd,
        "leading_even_threshold_coefficient": C + a * y * y / (48 * sp.pi**2 * m * m),
        "independent_running_threshold_dictionary": sp.simplify(
            running_reference - a * reference / (96 * sp.pi**2)
        ),
        "reference_log_derivative_without_branch_identity": sp.factor(
            sp.diff(sp.log(determinant) - reference, phi)
        ),
        "field_remainder_anchor": (2 * sp.log(1 - z) + 2 * z).subs(z, 0),
        "field_remainder_monotone_bound": sp.factor(
            sp.diff(z * z / (1 - z) + 2 * sp.log(1 - z) + 2 * z, z)
            - z * z / (1 - z) ** 2
        ),
    }
    for n in range(1, 9):
        checks[f"direct_log_field_coefficient_{n}"] = sp.diff(
            2 * sp.log(1 - z), z, n
        ).subs(z, 0) / sp.factorial(n) + sp.Rational(2, n)
    return {
        "flavor_mass_matrix": mass,
        "mass_squared_determinant": determinant,
        "vacuum_subtracted_log_determinant": reference,
        "leading_Phi_squared_F_squared_coefficient": C,
        "log_field_coefficients": coefficients,
        "field_series_absolute_remainder_upper": z * z / (1 - z),
        "field_domain": "z=y^2 Phi^2/m_F^2 in [0,1); for z<=1/4 the logarithm remainder after -2z is at most 4z^2/3",
        "momentum_scope": "Leading one-fermion-loop Wilson coefficient at zero external momenta. The field-series bound is NOT a momentum/derivative remainder estimate and does not bound the full finite-mass on-shell amplitude.",
        "checks": checks,
    }
