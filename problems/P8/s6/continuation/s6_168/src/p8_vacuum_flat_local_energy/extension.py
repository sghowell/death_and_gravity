"""A separately named high-field finite mass-reference extension."""

from functools import cache

import sympy as s
from p8_vacuum_protected_yukawa_profile import order


def loop_floor(external):
    if type(external) is not int or not 0 <= external <= 10:
        raise TypeError("Require a native source-endpoint count from zero through ten")
    return int(s.ceiling(6 - s.Rational(external, 2)))


@cache
def data():
    phi, R, eps, k = s.symbols("Phi R epsilon k", positive=True)
    scale = s.Symbol("field_count")
    chi_squared = phi**2 * (1 + (phi / R) ** 8) ** (-s.Rational(1, 4))
    expansion = s.series(chi_squared.subs(phi, scale * phi), scale, 0, 19).removeO()
    I, V, E, J = s.symbols("I V E J")
    return {
        "branch_name": "GY14-SAT8-MR",
        "finite_functional_definition": "At formal loop order one replace only the finite fermion pole-mass reference potential -fF(1)Phi^2/2 by -fF(1)chi(Phi)^2/2, chi=f_R(Phi). Its coefficient is the previously fixed mass anchor, not fitted to a clock energy. The fermion determinant and its minimal pole counterterms already use the same chi. Keep the zero-field vacuum reference and minimal field/quartic finite prescription.",
        "not_a_frozen_rewrite": "This selects a separately named higher-field EFT prescription. Original GY14, frozen SAT8, and all earlier source/report bytes are unchanged. It is not inferred that SAT8 had already fixed this finite high-field extension.",
        "first_changed_finite_vertex": {
            "scalar_valence": 10,
            "formal_loop_grade": 1,
            "weighted_excess": 5,
        },
        "source_count_loop_floors": {e: loop_floor(e) for e in (0, 1, 2, 4, 8, 10)},
        "same_regulated_lift": "Use chi_D=Phi_D/[1+(Phi_D/R_D)^8]^(1/8), R_D=mu^-epsilon R. Apply the argument replacement to the inherited dimensional mass-reference coefficient before Laurent finite products; do not discard its positive-epsilon terms. The k-th mass-reference difference has dimension2-8k+8k epsilon, exactly D-(8k+2)(1-epsilon).",
        "finite_order_invariance": "The new loop-one degree-ten vertex has weighted excess5. L_eff=1-E/2+sum[(n-2)/2+j] plus nonnegative physical-map degree. Forest contraction preserves that excess. It cannot contribute to named E<=4 data through loop two, including physical Phi pole/residue/four-point, vacuum, H source and the low cut. No complete higher-field quantum functional or higher-loop invariance follows.",
        "classical_match": "The change carries one formal loop power, so the complete classical zero-fermion S6.164 full-target action match is unchanged. This is not a quantum target match or a parent clock solution.",
        "unsaturated_counterterm_control": "Along Phi=sqrt(kappa)t, the original -fF(1)Phi^2/2 grows without bound while the protected-mass determinant stays bounded. A bounded mass alone therefore cannot bound that unsaturated one-loop reference contribution. This is not an exclusion of the entire parent theory.",
        "checks": {
            "first_changed_mass_reference_field_degree": s.expand(expansion).coeff(
                scale, 10
            )
            + phi**10 / (4 * R**8),
            "next_changed_mass_reference_field_degree": s.expand(expansion).coeff(
                scale, 18
            )
            - 5 * phi**18 / (32 * R**16),
            "all_lower_mass_reference_field_jets_unchanged": sum(
                s.expand(expansion - scale**2 * phi**2).coeff(scale, n)
                for n in range(10)
            ),
            "first_counterterm_weight": s.Rational(10 - 2, 2) + 1 - 5,
            "four_source_change_loop_floor": loop_floor(4) - 4,
            "two_source_change_loop_floor": loop_floor(2) - 5,
            "vacuum_change_loop_floor": loop_floor(0) - 6,
            "ten_source_change_allowed_at_one_loop": loop_floor(10) - 1,
            "forest_contraction_weight_identity": s.expand(
                (2 * I + E - 2 * V) / 2 + J - ((E - 2) / 2 + I - V + 1 + J)
            ),
            "dimensional_mass_reference_coefficient": s.expand(
                4 - 2 * eps - (8 * k + 2) * (1 - eps) - (2 - 8 * k + 8 * k * eps)
            ),
            "same_inherited_connected_loop_identity": order.grade(4, (10,), (1,)) - 4,
        },
        "gates": {
            "all_named_up_to_four_source_two_loop_data_protected": all(
                loop_floor(e) > 2 for e in (0, 1, 2, 4)
            ),
            "new_higher_field_one_loop_data_not_identified": loop_floor(10) == 1,
            "finite_mass_reference_not_retuned_to_cancel_energy": True,
        },
    }
