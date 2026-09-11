"""The source square cancels a reducible graph only with its mass insertion."""

from functools import cache

import sympy as s


@cache
def data():
    h, e, G, M, Phi, H, T0, T1, T2, ell, Q = s.symbols(
        "h epsilon G M Phi H T0 T1 T2 ell Q"
    )
    Tad = T0 / e + T1 + e * T2
    J = -G * Tad / 2
    potential = M * H**2 / 2 + H * (h * J + G * Phi**2 / 2)
    eliminated = s.expand(potential.subs(H, -(h * J + G * Phi**2 / 2) / M))
    dumbbell = -(G**2) * Tad**2 / (8 * M)
    source_mass = -G * J * Tad / (2 * M)
    source_square = -(J**2) / (2 * M)
    alpha1 = ell + 1
    alpha2 = 1 + ell + ell**2 / 2 + s.pi**2 / 12
    tad = -1 / (Q * e) - alpha1 / Q - e * alpha2 / Q
    finite_tad_square = s.expand(tad**2).coeff(e, 0)
    finite_source_square = -(G**2) * finite_tad_square / (8 * M)
    wrong_source_square = -(G**2) * alpha1**2 / (8 * M * Q**2)
    return {
        "full_regulated_first_source": J,
        "exact_H_eliminated_action": eliminated,
        "three_regulated_vacuum_occurrences": {
            "heavy_reducible_dumbbell": dumbbell,
            "source_induced_mass_insertion": source_mass,
            "source_square": source_square,
        },
        "sum_of_three_occurrences": s.factor(dumbbell + source_mass + source_square),
        "finite_first_source": G * (ell + 1) / (2 * Q),
        "finite_first_source_square_contribution": finite_source_square,
        "early_finite_part_source_square_defect": s.factor(
            finite_source_square - wrong_source_square
        ),
        "source_second_order_scope": "J2 is still an open full source-reference coefficient. Its square and cross term first enter the vacuum at h4 and h3, respectively; it does not affect this h2 vacuum result.",
        "checks": {
            "regulated_stationary_H_one_point": J + G * Tad / 2,
            "literal_source_induced_mass": s.diff(eliminated, Phi, 2).subs(Phi, 0)
            + h * G * J / M,
            "literal_source_square": eliminated.subs(Phi, 0) + h * h * J * J / (2 * M),
            "three_occurrences_cancel_before_epsilon_limit": s.factor(
                dumbbell + source_mass + source_square
            ),
            "finite_source_square_keeps_positive_epsilon": finite_tad_square
            - (2 * ell**2 + 4 * ell + 3 + s.pi**2 / 6) / Q**2,
            "early_finite_part_defect_is_pole_times_epsilon": s.factor(
                finite_source_square
                - wrong_source_square
                + G**2 * alpha2 / (4 * M * Q**2)
            ),
            "first_and_second_source_square_orders": s.expand(
                (h * s.Symbol("J1") + h * h * s.Symbol("J2")) ** 2
            ).coeff(h, 2)
            - s.Symbol("J1") ** 2,
            "no_second_source_in_vacuum_at_order_two": s.diff(
                s.expand((h * s.Symbol("J1") + h * h * s.Symbol("J2")) ** 2).coeff(
                    h, 2
                ),
                s.Symbol("J2"),
            ),
        },
    }
