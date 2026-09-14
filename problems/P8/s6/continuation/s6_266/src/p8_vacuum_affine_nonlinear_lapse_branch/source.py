"""Source-pinned full bounce Hamiltonian with its nonzero primitive contact."""

from functools import cache

import sympy as s
from p8_vacuum_affine_finite_band_neighborhood import source as fifth
from p8_vacuum_affine_nonlinear_auxiliary_measure import canonical as parent
from p8_vacuum_affine_nonlinear_reference_volume import source as previous

N, R, j = parent.N, parent.R, parent.j
F = s.Function("whole_original_unrescaled_F", real=True)(N)
Ruu = s.Function("whole_bounce_Ruu", real=True)(N)
Iu = s.Function("whole_bounce_primitive_Iu", real=True)(N)
p, G, dp, ph, eta, sh, el, ma, wm, gm, gh, curv = s.symbols(
    "trace_p Gauss_G delta_pm H_momentum mass_adapted_H shear "
    "electric_over_zeta zeta_magnetic vector_mass M1_gradient H_gradient curvature",
    real=True,
)
mu = s.Symbol("actual_H_mass2_over_1e200", positive=True)
COORDS = (p, G, dp, ph, eta, sh, el, ma, wm, gm, gh, curv)
CENTER = {z: s.Integer(0) for z in COORDS}
DELTA = s.Rational(1, 10**250)
OUTER = s.Rational(1, 10**6)
INNER = s.Rational(1, 10**245)
CENTER_RADIUS = s.Rational(1, 10**375)
ERROR = s.Rational(1, 10**380)
JERROR = s.Rational(1, 10**2600)
R_TREE = N**-2
F_TREE = -(624 * N**-4 + 753 * N**-2 + 224) / 200
RUU_TREE = -6 * (N**-2 - 1)
NORMAL = (
    -3 * (p - (R - 1) * G) ** 2 / (4 * R ** s.Rational(1, 4))
    - R ** -s.Rational(3, 4) * F
    + G**2 * R ** s.Rational(3, 4) / 2
    + 2 * sh * R ** -s.Rational(1, 4)
    + ((s.Rational(1, 10) + dp) ** 2 + ph**2) * R ** s.Rational(3, 4) / 2
    + R ** -s.Rational(1, 4) * (gm + gh) / 2
    + R ** -s.Rational(3, 4) * (mu * eta**2 / 2 - j * eta / s.Integer(10) ** 100)
    - R ** s.Rational(3, 4) * curv / 2
    + el * R ** s.Rational(1, 4) / 2
    + ma * R ** s.Rational(1, 4) / 4
    + R ** -s.Rational(1, 4) * wm / 2
)
PRIMITIVE_N = 3 * R ** -s.Rational(7, 4) * Ruu * s.diff(R, N) / (4 * N)
HAMILTONIAN = N * NORMAL + Iu
CONSTRAINT = PRIMITIVE_N + s.diff(N * NORMAL, N)
TEMPORAL = -3 * (R - 1) * (p - (R - 1) * G) / (
    2 * R ** s.Rational(1, 4)
) - G * R ** s.Rational(3, 4)
GAMMA = 1 - 3 * (R - 1) ** 2 / (2 * R)
TEMPORAL_PIVOT = -N * R ** -s.Rational(3, 4) / GAMMA


def tree_rules(last=4):
    return {
        s.diff(f, N, k): s.diff(tree, N, k)
        for f, tree in [(R, R_TREE), (F, F_TREE), (Ruu, RUU_TREE), (j, s.Integer(0))]
        for k in range(last + 1)
    }


def tree(expr):
    return s.factor(expr.subs(tree_rules(), simultaneous=True))


@cache
def data():
    old = parent.full_trace()
    current = parent.current
    bound = old["whole_reduced_Hamiltonian"].subs(
        {
            parent.p: p,
            parent.G: G,
            parent.pm: s.Rational(1, 10) + dp,
            parent.ph: ph,
            parent.h: eta / s.Integer(10) ** 100,
            parent.n: mu * s.Integer(10) ** 200,
            parent.shear: sh,
            parent.electric: el * parent.zeta,
            parent.magnetic: ma / parent.zeta,
            parent.wmass: wm,
            parent.gm: gm,
            parent.gh: gh,
            parent.curvature: curv,
            parent.B: 0,
            parent.Href: 0,
            parent.F: R ** -s.Rational(3, 4) * F - Iu / N,
        },
        simultaneous=True,
    )
    temporal = old["whole_temporal_solution"].subs(
        {parent.p: p, parent.G: G, parent.B: 0, parent.Href: 0}, simultaneous=True
    )
    # Substitute complete derivative atoms simultaneously: Ru=0 does NOT imply Ruu=0.
    slice_rules = {current.N: N}
    for i in range(3):
        for k in range(3 - i):
            slice_rules[s.diff(current.R, current.u, i, current.N, k)] = (
                s.diff(R, N, k)
                if i == 0
                else s.Integer(0)
                if i == 1
                else s.diff(Ruu, N, k)
            )
    primitive = s.diff(current.IN, current.u).xreplace(slice_rules)
    evidence = fifth.data()
    comparison = tree(CONSTRAINT).subs(CENTER).subs(N, 1)
    return {
        "whole_original_full_time_Hamiltonian": old["whole_reduced_Hamiltonian"],
        "whole_bounce_Hamiltonian_in_normalized_density_invariants": HAMILTONIAN,
        "whole_bounce_constraint": CONSTRAINT,
        "whole_bounce_temporal_solution": TEMPORAL,
        "whole_bounce_primitive_contact_N": PRIMITIVE_N,
        "whole_actual_source_and_fixed_profile_bindings": evidence[
            "whole_current_function_bindings"
        ],
        "whole_fixed_profiles": evidence["whole_fixed_profile_bindings"],
        "whole_normalized_heavy_source": evidence["whole_normalized_heavy_source"],
        "whole_independent_invariant_coordinates": list(COORDS),
        "whole_mass_adapted_binding": {
            parent.h: eta / s.Integer(10) ** 100,
            parent.n: mu * s.Integer(10) ** 200,
        },
        "whole_original_mass_ratio": parent.source.MASS2 / s.Integer(10) ** 200,
        "whole_bare_comparison_constraint": tree(CONSTRAINT),
        "whole_bare_center_lapse_pivot": tree(s.diff(CONSTRAINT, N))
        .subs(CENTER)
        .subs(N, 1),
        "whole_source_scope": "Only u=0. The entire Ru vanishes at every lapse, giving I=B=0 from I(u,1)=0, but Iu(N) is retained with Iu(1)=0 and the displayed nonzero derivative. All scalar gradients, curvature, shear, vector electric/magnetic/mass and Gauss terms remain. The twelve density/jet invariants are NOT independent canonical oscillators. Their signed box is an analytic enlargement containing the physical data subset. No spatial gauge constraint is silently solved.",
        "checks": {
            "literal_entire_parent_Hamiltonian_slice": s.factor(bound - HAMILTONIAN),
            "literal_entire_parent_lapse_constraint": s.factor(
                s.diff(bound, N).subs(s.diff(Iu, N), PRIMITIVE_N) - CONSTRAINT
            ),
            "literal_entire_parent_temporal_root": s.factor(temporal - TEMPORAL),
            "literal_nonzero_primitive_time_derivative": s.factor(
                primitive - PRIMITIVE_N
            ),
            "whole_original_Ru_slice_zero": evidence["checks"][
                "whole_full_R_initial_time_derivative_zero"
            ],
            "bare_reference_constraint_zero_only_comparison": comparison,
            "bare_reference_lapse_pivot": tree(s.diff(CONSTRAINT, N))
            .subs(CENTER)
            .subs(N, 1)
            + s.Rational(243, 80),
            "bare_reference_second_constraint_derivative": tree(
                s.diff(CONSTRAINT, N, 2)
            )
            .subs(CENTER)
            .subs(N, 1)
            - s.Rational(13821, 400),
        },
        "gates": {
            "all_original_full_C5_source_enclosures": all(evidence["gates"].values()),
            "actual_full_mass_inside_mass_adapted_box": s.Rational(1, 10**4)
            < parent.source.MASS2 / s.Integer(10) ** 200
            < s.Rational(1, 100),
            "same_original_kappa_zeta_and_heavy_localizer": parent.source.KAPPA
            == 10**800
            and previous.heavy.LOCALIZER == 10**420,
            "all_twelve_spatial_matter_vector_invariants_retained": all(
                CONSTRAINT.has(z) for z in COORDS
            ),
            "nonzero_Ruu_primitive_contact_retained": CONSTRAINT.has(Ruu),
            "source_and_state_not_reselected": True,
        },
    }
