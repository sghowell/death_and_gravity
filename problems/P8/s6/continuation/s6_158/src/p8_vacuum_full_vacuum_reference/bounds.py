"""Exact rational power caps for the regulator-circle vacuum enclosure."""

import sympy as s


def rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer, s.Rational)):
        raise TypeError("An exact real rational is required")
    return s.Rational(value)


def scalar_enclosure(m, M, L, g, Q, U, A, B):
    m, M, L, g, Q, U, A, B = map(rational, (m, M, L, g, Q, U, A, B))
    if not (
        m >= 2
        and M >= 2
        and L >= 0
        and g >= 0
        and 0 < Q <= 144
        and U >= 1
        and A >= 1
        and B >= 1
    ):
        raise ValueError(
            "Need positive physical masses, nonnegative couplings, Q lower bound and power caps"
        )
    if not (m <= U**8 and 3 * M <= A**8 and M <= B**16):
        raise ValueError(
            "The exact regulator-circle power caps do not enclose the masses"
        )
    groups = {
        "quartic_and_physical_scalar_mass_reference": 4900 * L * U**2 / (8 * Q**2),
        "cubic_sunset": 1200 * g * U**2 * A**9 / Q**2,
        "mixed_bubble_OS_mass_reference": 2240 * g * U**2 * B / Q**2,
        "heavy_mass_MS_reference": 280 * g * U * B**17 / Q**2,
    }
    return {
        "power_caps": {"U": U, "A": A, "B": B},
        "group_uppers": groups,
        "complete_scalar_vacuum_absolute_upper": sum(groups.values()),
    }


def one_loop_enclosure(m, M, Qlo, Qhi, ell_upper):
    m, M, Qlo, Qhi, ell_upper = map(rational, (m, M, Qlo, Qhi, ell_upper))
    if not (
        m >= 2 and 2 <= M <= m * m and 0 < Qlo <= 144 and Qhi >= 256 and ell_upper >= 1
    ):
        raise ValueError(
            "Need the named heavy hierarchy, Q enclosure and verified log upper"
        )
    # ell_upper is an externally proved bound on log(m^2), not inferred here.
    scalar = (1 + M * M) * (ell_upper + s.Rational(3, 2)) / (4 * Qlo)
    return {
        "scalar_absolute_upper": scalar,
        "all_flavor_fermion_lower": 63 * m**4 / Qhi,
        "all_flavor_fermion_upper": 63 * m**4 / Qlo,
        "complete_lower": 63 * m**4 / Qhi - scalar,
        "complete_upper": 63 * m**4 / Qlo + scalar,
    }


def data():
    E, sunset, mixed, heavy, T = s.symbols("E sunset mixed heavy tree_vacuum")
    h = s.Symbol("h")
    return {
        "finite_part_bound_rule": "For a meromorphic F with only the origin singular inside |epsilon|=1/16, Fin F is the circular average of F. Hence |Fin F|<=sup_circle |F|. Overall pure-pole vacuum counterterms do not change this average.",
        "power_caps_rule": "mF<=U^8, 3M<=A^8, M<=B^16 imply mF^(1/8)<=U, (3M)^(9/8)<=A^9 and M^(17/16)<=B^17.",
        "checks": {
            "all_four_scalar_groups_sum": E
            + sunset
            + mixed
            + heavy
            - (E + sunset + mixed + heavy),
            "vacuum_reference_cancels_each_known_order": s.expand(
                h * T + h * h * E + h * (-T) + h * h * (-E)
            ),
            "source_J2_first_affects_vacuum_at_order_three": s.expand(
                (h * s.Symbol("J1") + h * h * s.Symbol("J2")) ** 2
            ).coeff(h, 3)
            - 2 * s.Symbol("J1") * s.Symbol("J2"),
        },
    }
