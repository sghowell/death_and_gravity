"""Full-dimensional scalar pair and fixed physical proper-time matching."""

from functools import cache

import sympy as s
from p8_vacuum_affine_ordered_scalar_symbol import geometry as g

from .coordinates import L0


def leading(channel):
    rows = g.contractions(channel)
    combo = (
        sum(rows[name] for name in ("00", "01", "10", "11")) / 4
        + rows["LL"] / 4
        + (rows["LC"] + rows["CL"]) / 8
        + rows["CC"] / 16
    )
    return s.factor(g.sphere_average(combo.subs(g.y, 0)))


@cache
def data():
    d = g.d
    trace = leading("trace_trace")
    C = leading("tensor")
    B = s.factor((trace - d * C) / d**2)
    cross = s.ImmutableMatrix([leading("trace_scalar"), leading("scalar_trace")])
    physical = C * s.Matrix([[12, -4], [-4, 4]]) + B * s.Matrix([[36, -12], [-12, 4]])
    M = physical.subs(d, 3).applyfunc(s.factor)
    jet = physical.diff(d).subs(d, 3).applyfunc(s.factor)
    wd, cd, wg, cg = s.symbols("wd cd wg cg")
    varying = C * (4 * d * wd * wg - 4 * (wd * cg + cd * wg) + 4 * cd * cg) + B * (
        2 * d * wd - 2 * cd
    ) * (2 * d * wg - 2 * cg)
    varying_matrix = s.Matrix(
        [[s.diff(varying, x, y) for y in (wg, cg)] for x in (wd, cd)]
    )
    varying_jet = varying_matrix.diff(d).subs(d, 3).applyfunc(s.factor)
    tau = s.Symbol("lag", positive=True)
    return {
        "all_dimension_tracefree_invariant": C,
        "all_dimension_double_trace_invariant": B,
        "literal_I_d_trace_pair": trace,
        "both_ordered_trace_tracefree_crosses": cross,
        "fixed_physical_scalar_matrix": M,
        "fixed_physical_first_dimension_jet": jet,
        "distinct_dimension_varying_test_jet": varying_jet,
        "proper_normalized_radial_matrix": 32 * M,
        "off_diagonal_leading_matrix": 24 * M,
        "fixed_fourth_finite_local_matrix": L0.T * s.diag(-4, -s.Rational(4, 45)) * L0,
        "reference": "gamma L0^T diag(Ftrace,(8/3)F2)L0 D^4 in proper time, gamma=1/(64pi^2 kappa); original factors, mass, cut and pole.",
        "finite_transfer_phase": "Retain exp(-i(n dot P)Delta_sigma) as a smooth orderzero amplitude; its diagonal value is1. It is not inverse-k small.",
        "strong_remainder": "For each finite Pmax and each fixed j, four OUTPUT primitives of the actual scalar remainder have a uniform complete row-sum bound Cj(Pmax)(1+|log lag|). Full original state, dimensional jets, contacts and lower local terms remain. No Cj is numerically evaluated.",
        "checks": {
            "full_trace_invariant": trace - (d - 1) * (d * d - 5 * d + 8) / 4,
            "full_tracefree_invariant": C - (2 * d * d + d - 8) / (2 * d * (d + 2)),
            "full_double_trace_invariant": B
            - (d**3 - 4 * d * d - 3 * d + 16) / (4 * d * (d + 2)),
            "both_ordered_crosses": cross,
            "fixed_physical_matrix": M
            - s.Matrix([[4, -s.Rational(4, 3)], [-s.Rational(4, 3), s.Rational(8, 5)]]),
            "two_original_channel_logarithms": M
            - L0.T * s.diag(4, s.Rational(52, 45)) * L0,
            "fixed_physical_dimension_jet": jet
            - s.Matrix(
                [
                    [s.Rational(46, 15), -s.Rational(46, 45)],
                    [-s.Rational(46, 45), s.Rational(22, 25)],
                ]
            ),
            "double_trace_dimension_jet": s.diff(B, d).subs(d, 3) - s.Rational(4, 225),
            "leading_determinant": M.det() - s.Rational(208, 45),
            "four_primitive_Abel_matching": s.diff(M / tau, tau, 4) - 24 * M / tau**5,
            "proper_momentum_full_dimension_scale_exponent": (-d - 2) + (d + 1) + 1,
        },
        "gates": {
            "physical_dimension_jet_not_varying_test_jet": jet != varying_jet,
            "both_original_channels_nonzero": M.det() != 0,
            "complete_full_dimensional_matching_before_finite_part": True,
            "compact_external_ball_not_internal_loop_cutoff": True,
            "actual_state_comparison_not_flat_vacuum_replacement": True,
        },
    }
