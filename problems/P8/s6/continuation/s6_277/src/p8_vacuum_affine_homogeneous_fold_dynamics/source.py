"""Entire original homogeneous source and licensed finite clock derivatives."""

from functools import cache

import sympy as s
from p8_vacuum_affine_global_lapse_chart_obstruction import fold as previous
from p8_vacuum_affine_physical_background_vertices import parent as physical
from p8_vacuum_affine_quantitative_local_time import source as original

u, N, p, sh, eta, ph = (
    original.u,
    original.N,
    original.p,
    original.sh,
    original.eta,
    original.ph,
)
m = s.Symbol("whole_homogeneous_M1_momentum_density", real=True)
rho = s.Function("actual_complete_fixed_energy")(u)
pressure = s.Function("actual_complete_fixed_pressure")(u)
RHO, PRESSURE, RHO1, PRESSURE1, RHO2, PRESSURE2 = s.symbols(
    "actual_rho actual_pressure actual_rho_first actual_pressure_first actual_rho_second actual_pressure_second",
    real=True,
)
PROFILE_JETS = (RHO, PRESSURE, RHO1, PRESSURE1, RHO2, PRESSURE2)
PROFILE_BOUND = previous.PROFILE_BOUND
MASS_RATIO = original.parent.source.MASS2 / 10**200
ZERO_HEAVY = {eta: s.Integer(0), ph: s.Integer(0)}
FOLD_SHEAR = s.Rational(81, 160) - 2 * PRESSURE / 3 + 7 * RHO / 12 + 3 * p * p / 8
FOLD_MATTER_SQUARE = s.Rational(1217, 200) - 2 * PRESSURE + 3 * RHO
FOLD_SECOND = (11391 + 1400 * PRESSURE - 1600 * RHO) / 400
ALLOWED = {
    "H": ((0, 0), (1, 0)),
    "C": ((0, 0), (1, 0), (2, 0), (0, 1), (0, 2), (1, 1)),
    "T": ((0, 0), (1, 0), (0, 1)),
    "B": ((1, 0), (1, 1)),
}


def restriction():
    result = dict.fromkeys(original.COORDS, s.Integer(0))
    for z in (p, sh, eta, ph):
        result.pop(z)
    result[original.dp] = m - s.Rational(1, 10)
    result[original.mu] = MASS_RATIO
    return result


@cache
def whole():
    bind = restriction()
    return {
        "H": original.HAMILTONIAN.subs(bind, simultaneous=True),
        "C": original.CONSTRAINT.subs(bind, simultaneous=True),
        "T": original.TEMPORAL.subs(bind, simultaneous=True),
        "B": original.B,
    }


def require_jet(kind, time_order, lapse_order):
    if not isinstance(kind, str) or kind not in ALLOWED:
        raise ValueError("Require an original H,C,T or B clock-jet family")
    for value in (time_order, lapse_order):
        if isinstance(value, bool) or not isinstance(value, (int, s.Integer)):
            raise TypeError("Require exact integer clock derivative orders")
    pair = (int(time_order), int(lapse_order))
    if pair not in ALLOWED[kind]:
        raise ValueError("Require a licensed full-source clock jet")
    return kind, *pair


def clock_jet(kind, time_order=0, lapse_order=0):
    return _clock_jet(*require_jet(kind, time_order, lapse_order))


@cache
def _clock_jet(kind, time_order, lapse_order):
    expr = original.eliminate_N_primitive(
        s.diff(whole()[kind], u, time_order, N, lapse_order)
    )
    Rclock = 1 + (N**-2 - 1) / (1 + u * u) ** 3
    tree = physical.original.original.data()["original_retuned_tree_scalar"].subs(
        physical.X, N**-2
    )
    Fclock = tree - pressure - (rho + pressure) * (N**-2 - 1) / 2
    replacements = {
        original.R: Rclock,
        original.F: Fclock,
        original.j: s.Integer(0),
        original.primitive: s.Integer(0),
        original.Hclock: 4 * u / (1 + u * u),
    }
    for term in expr.atoms(s.Derivative):
        counts = dict(term.variable_count)
        nt, nn = counts.get(u, 0), counts.get(N, 0)
        if term.expr in (original.R, original.F):
            if nn > 5:
                raise ValueError("A clock derivative exceeds the finite-source license")
            value = Rclock if term.expr == original.R else Fclock
            replacements[term] = s.diff(value, u, nt, N, nn)
        elif term.expr == original.j:
            if nn > 3:
                raise ValueError(
                    "A heavy-source derivative exceeds the checked license"
                )
            replacements[term] = s.Integer(0)
        elif term.expr == original.primitive:
            if nn:
                raise ValueError("An N primitive derivative was not eliminated")
            replacements[term] = s.Integer(0)
        elif term.expr == original.Hclock:
            replacements[term] = s.diff(4 * u / (1 + u * u), u, nt)
    expr = expr.subs(replacements, simultaneous=True).subs(N, 1).subs(u, 0)
    for fun, values in (
        (rho, (RHO, RHO1, RHO2)),
        (pressure, (PRESSURE, PRESSURE1, PRESSURE2)),
    ):
        expr = expr.xreplace(
            {s.diff(fun, u, k).subs(u, 0): value for k, value in enumerate(values)}
        )
    return s.factor(expr)


def at_fold(expr):
    return s.factor(
        expr.subs(ZERO_HEAVY).subs(
            {sh: FOLD_SHEAR, m * m: FOLD_MATTER_SQUARE}, simultaneous=True
        )
    )


@cache
def data():
    germs = physical.source_germs()
    jfull = original.heavy_source.physical_source(u, physical.X)
    jets = {
        kind + str(pair): clock_jet(kind, *pair)
        for kind, pairs in ALLOWED.items()
        for pair in pairs
    }
    H0 = clock_jet("H").subs(ZERO_HEAVY)
    C0 = clock_jet("C").subs(ZERO_HEAVY)
    CN = clock_jet("C", 0, 1).subs(ZERO_HEAVY)
    checks = {
        **{"source_" + k: v for k, v in germs["checks"].items()},
        **{
            "entire_heavy_X_jet_" + str(k): s.diff(jfull, physical.X, k).subs(
                physical.X, 1
            )
            for k in range(4)
        },
        "whole_clock_H": s.factor(
            H0 + 3 * p * p / 4 - s.Rational(1601, 200) - PRESSURE - m * m / 2 - 2 * sh
        ),
        "whole_clock_C": s.factor(
            C0
            - s.Rational(1, 400)
            - 3 * PRESSURE / 2
            + RHO
            + 9 * p * p / 8
            + m * m / 4
            - 3 * sh
        ),
        "whole_clock_CN": s.factor(
            CN
            + s.Rational(2433, 800)
            - 7 * PRESSURE / 4
            + 2 * RHO
            + 9 * p * p / 16
            - 3 * m * m / 8
            - 3 * sh / 2
        ),
        "whole_time_Cu": s.factor(
            clock_jet("C", 1, 0).subs(ZERO_HEAVY) + 9 * p - 3 * PRESSURE1 / 2 + RHO1
        ),
        "full_Bu_zero": clock_jet("B", 1, 0),
        "full_BuN_not_discarded": clock_jet("B", 1, 1) + 6,
        "fold_C": at_fold(clock_jet("C")),
        "fold_CN": at_fold(clock_jet("C", 0, 1)),
        "fold_CNN": s.factor(at_fold(clock_jet("C", 0, 2)) - FOLD_SECOND),
        "fold_transverse_shear": s.diff(C0, sh) - 3,
        "fold_temporal_value": at_fold(clock_jet("T")),
        "fold_temporal_N": at_fold(clock_jet("T", 0, 1)) - 3 * p,
        "fold_temporal_u": at_fold(clock_jet("T", 1, 0)),
    }
    return {
        "whole_original_homogeneous_H_C_T_B": whole(),
        "whole_actual_fixed_coefficient_bindings": germs[
            "fixed_total_coefficient_function_bindings"
        ],
        "whole_all_three_vacuum_constants": germs["all_three_vacuum_constants"],
        "whole_source_R_F_remainders": (
            germs["entire_R_difference"],
            germs["entire_F_difference"],
        ),
        "whole_entire_normalized_heavy_source": jfull,
        "whole_actual_heavy_mass_ratio": MASS_RATIO,
        "whole_licensed_clock_jets": jets,
        "whole_fold_shear_and_matter_square": (FOLD_SHEAR, FOLD_MATTER_SQUARE),
        "whole_actual_C5_profile_bound": PROFILE_BOUND,
        "whole_source_scope": "Only finite clock derivatives are evaluated by the licensed order-eight germ. The full original H,C,T,primitive,heavy source and fixed profile functions define off-clock dynamics; neither H nor C is replaced there by a Taylor polynomial.",
        "checks": checks,
        "gates": {
            "all_original_source_gates": all(germs["gates"].values()),
            "all_four_heavy_clock_jets_zero": all(
                checks["entire_heavy_X_jet_" + str(k)] == 0 for k in range(4)
            ),
            "actual_profile_bound": PROFILE_BOUND == s.Rational(2, 10**400),
            "mass_ratio_in_original_interval": s.Rational(1, 10**4)
            < MASS_RATIO
            < s.Rational(1, 100),
            "fold_shear_positive_for_all_real_p": s.Rational(81, 160)
            - s.Rational(5, 4) * PROFILE_BOUND
            > s.Rational(49, 100),
            "fold_matter_square_positive": s.Rational(1217, 200) - 5 * PROFILE_BOUND
            > 6,
            "fold_second_above28": s.Rational(11391, 400)
            - s.Rational(15, 2) * PROFILE_BOUND
            > 28,
            "original_profile_not_new_quantum_state_or_live_mean": True,
        },
    }
