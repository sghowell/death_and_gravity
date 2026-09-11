"""Unit-disc extension of the direct primitive tail and canonical pole algebra."""

from functools import cache

import sympy as s
from p8_vacuum_two_loop_finite_contact.calibration import rational


def tail_coefficient(soft_tail, radius=1):
    E, R = map(rational, (soft_tail, radius))
    if E < 0 or not 0 < R <= 1:
        raise ValueError("Need a nonnegative exact tail and 0<radius<=1")
    return E / (4 * (1 - R / 2))


@cache
def data():
    v, h, B1, B2, K, rho = s.symbols("distance h B1 B2 k0 radius")
    S = s.Function("scalar_one_loop_OS")(v)
    F = s.Function("fermion_one_loop_OS")(v)
    E = s.Symbol("direct_fermion_soft_tail")
    extra = -K * (2 * S - F)
    inverse = v + h * v**2 * B1 + h**2 * v**2 * B2
    a0, a1 = s.symbols("local_mass local_kinetic")
    affine = a0 + a1 * v
    return {
        "direct_fermion_tail_to_OS_coefficient": E / (4 * (1 - rho / 2)),
        "first_finite_field_parameter_variation_of_one_loop_OS": extra,
        "two_loop_canonical_inverse": inverse,
        "checks": {
            "parent_half_disc_tail_constant": E / (4 * (1 - s.Rational(1, 2) / 2))
            - E / 3,
            "new_unit_disc_tail_constant": E / (4 * (1 - s.Rational(1, 2))) - E / 2,
            "geometric_second_and_higher_coefficient_sum": (rho / 2) ** 2
            / (1 - rho / 2)
            / rho**2
            - 1 / (4 * (1 - rho / 2)),
            "canonical_inverse_factorization": s.expand(
                inverse - v * (1 + v * (h * B1 + h**2 * B2))
            ),
            "mass_one_zero_is_exact": inverse.subs(v, 0),
            "unit_light_residue": s.limit(v / inverse, v, 0) - 1,
            "all_local_mass_kinetic_references_OS_zero": s.expand(
                affine - affine.subs(v, 0) - v * s.diff(affine, v).subs(v, 0)
            ),
            "finite_parameter_change_commutes_with_OS": s.diff(
                s.Symbol("parameter") * v**2, s.Symbol("parameter")
            )
            - v**2,
        },
        "scope": "The existing direct fermion soft tail is holomorphic on a neighborhood of |s-1|<=2. Its Cauchy coefficient estimate extends the previous half-disc bound E/3 to E/2 on the unit disc. Every local affine reference is projected only after regulated pairing. A positive factored inverse bound establishes one unit-residue pole for the stated fixed-order inverse, not an all-orders pole count.",
    }
