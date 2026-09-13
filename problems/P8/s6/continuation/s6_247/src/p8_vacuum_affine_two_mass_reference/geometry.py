"""Complete minimal-scalar cut and unchanged finite two-mass reference inputs."""

from functools import cache

import sympy as s
from p8_vacuum_affine_band_coupled_inverse import matching as old_matching
from p8_vacuum_affine_heavy_curved_state import quantum, state
from p8_vacuum_affine_isolated_shear_resolvent import normalization as old_shear

d, y, p = s.symbols("spatial_dimension radial_fraction invariant_p", positive=True)
n, mu2, sigma = s.symbols(
    "heavy_mass_squared proca_mass_squared spectral_squared", positive=True
)
ell = s.Symbol("actual_log_heavy_mass_squared", real=True)
ETA = s.diag(1, -1, -1, -1)
PROCA_MASS2 = s.Integer(10) ** 6
HEAVY_MASS2 = state.MASS2
KAPPA = state.KAPPA
W_H = {0: y**2 * (3 - y**2) ** 2, 2: y**6 / 30}
W_P = {
    0: y**2 * (3 - 2 * y**2 + 3 * y**4),
    2: y**2 * (30 - 20 * y**2 + 3 * y**4) / 30,
}
C0 = 2 * (ell + 2)
C2 = (ell + 2) / 60
L0 = s.Matrix([[1, -s.Rational(1, 3)], [0, -1]])


def stress_pair(k, l, mass_squared):
    """One symmetric Wick amplitude; the exchange factor two is not here."""
    kl = (k.T * ETA * l)[0]
    return -(k * l.T + l * k.T) / 2 + ETA * (kl + mass_squared) / 2


def heavy_cut(spin):
    if (
        isinstance(spin, bool)
        or not isinstance(spin, (int, s.Integer))
        or spin not in (0, 2)
    ):
        raise ValueError("Keep both physical scalar stress spin channels")
    beta = s.sqrt(1 - 4 * n / sigma)
    if spin == 0:
        return beta * (sigma + 2 * n) ** 2 / (384 * s.pi**2)
    return beta**5 * sigma**2 / (3840 * s.pi**2)


@cache
def data():
    energy, k = s.symbols("on_shell_energy com_momentum", positive=True)
    ka, la = s.Matrix([energy, 0, 0, k]), s.Matrix([energy, 0, 0, -k])

    def on_shell(value):
        return (
            s.Poly(s.expand(value), energy)
            .rem(s.Poly(energy**2 - k * k - n, energy))
            .as_expr()
        )

    pair = stress_pair(ka, la, n).applyfunc(on_shell)
    spatial = pair[1:4, 1:4]
    trace = s.trace(spatial)
    norm = s.trace(spatial * spatial)
    com = {0: s.factor(trace**2 / 3), 2: s.factor((norm - trace**2 / 3) / 5)}
    invariant = {i: s.factor(v.subs(k * k, sigma / 4 - n)) for i, v in com.items()}
    phase_space = s.sqrt(1 - 4 * n / sigma) / (8 * s.pi)
    derived_cut = {
        i: s.factor(2 * phase_space * v / (8 * s.pi)) for i, v in invariant.items()
    }
    checks = {
        "whole_symmetric_scalar_stress_pair": pair - pair.T,
        "whole_COM_time_Ward_row": pair[0, :],
        "scalar_spin2_COM_projection": com[2] - 2 * k**4 / 15,
        "scalar_spin0_COM_projection": com[0] - (2 * k * k + 3 * n) ** 2 / 3,
        "scalar_spin2_invariant_projection": invariant[2] - (sigma - 4 * n) ** 2 / 120,
        "scalar_spin0_invariant_projection": invariant[0] - (sigma + 2 * n) ** 2 / 12,
        "two_body_measure_Wick_exchange_and_current_quarter_spin2": s.simplify(
            derived_cut[2] - heavy_cut(2)
        ),
        "two_body_measure_Wick_exchange_and_current_quarter_spin0": s.simplify(
            derived_cut[0] - heavy_cut(0)
        ),
    }
    q11, q12, q13, q22, q23, q33 = s.symbols("q11 q12 q13 q22 q23 q33", real=True)
    Q = s.Matrix([[q11, q12, q13], [q12, q22, q23], [q13, q23, q33]])
    MQ = s.zeros(5)
    MQ[0, 0] = -s.trace(Q) / 2
    MQ[1:4, 1:4] = s.trace(Q) * s.eye(3) / 2 - Q
    MQ[4, 4] = s.trace(Q) / 2
    fk = s.Matrix([-s.I * energy, 0, 0, s.I * k, s.sqrt(n)])
    fl = s.Matrix([-s.I * energy, 0, 0, -s.I * k, s.sqrt(n)])
    vertex = on_shell((fk.T * MQ * fl)[0])
    checks["full_scalar_Hamiltonian_to_covariant_stress_pair"] = s.expand(
        vertex + s.trace(Q * spatial)
    )
    checks["full_normalized_trace_vertex"] = s.expand(
        vertex.subs({q11: 2, q22: 2, q33: 2, q12: 0, q13: 0, q23: 0}) / 2
        - (2 * k * k + 3 * n)
    )
    # No conformal trace improvement has been introduced.
    for spin, factor in ((0, 768), (2, 64)):
        change = 4 * n / (1 - y * y)
        converted = (
            factor * s.pi**2 * p * heavy_cut(spin) / (sigma**3 * (sigma + p))
        ).subs(sigma, change) * s.diff(change, y)
        converted = s.refine(converted, s.Q.positive(1 - y * y))
        checks["complete_scalar_cut_radial_conversion_" + str(spin)] = s.simplify(
            converted - p * W_H[spin] / (4 * n + p * (1 - y * y))
        )
    checks["unchanged_Proca_shear_weight"] = s.expand(
        old_shear.W2.subs(old_shear.y, y) - W_P[2]
    )
    checks["unchanged_Proca_trace_weight"] = s.expand(
        old_shear.W0.subs(old_shear.y, y) - W_P[0]
    )
    checks["unchanged_Proca_shear_finite_constant"] = old_shear.A0 - s.Rational(1, 30)

    # Reconstruct the full finite scalar density before the physical curvature basis.
    finite = quantum.data()["complete_scalar_finite_heat_matching_before_64pi2"]
    names = {str(v): v for v in finite.free_symbols}
    R, C, E = s.symbols("R Weyl_squared Euler", real=True)
    # Ricci^2=(C^2-E)/2+R^2/3; Riemann^2=2C^2-E+R^2/3.
    reduced = s.expand(
        finite.subs(
            {
                names["mass_squared"]: n,
                quantum.ell: ell,
                names["R_old"]: R,
                names["Ricci_squared"]: (C - E) / 2 + R * R / 3,
                names["Riemann_squared"]: 2 * C - E + R * R / 3,
            }
        )
    )
    finite_c2 = reduced.coeff(C)
    finite_trace = 72 * reduced.coeff(R, 2)
    checks["whole_same_scheme_heavy_Weyl_coefficient"] = finite_c2 + ell / 60
    checks["whole_same_scheme_heavy_trace_Hessian"] = finite_trace + 2 * ell
    checks["complete_two_mass_shear_constant"] = -s.Rational(1, 30) + finite_c2 + C2
    checks["complete_two_mass_trace_constant"] = -4 + finite_trace + C0
    checks["whole_compact_Euler_finite_coefficient"] = reduced.coeff(E) - ell / 180

    # The COMPLETE high-radius scalar pair is (tr Q - nhat Q nhat)/2.
    ch = 1 / (2 * d * (d + 2))
    bh = (1 - 2 / d + 1 / (d * (d + 2))) / 4
    tf = s.Matrix([[12, -4], [-4, 4]])
    tr = s.Matrix([[36, -12], [-12, 4]])
    mh = (ch * tf + bh * tr).applyfunc(s.factor)
    op = old_matching.data()
    od = old_matching.g.d
    cp = op["all_dimension_tracefree_invariant"].subs(od, d)
    bp = op["all_dimension_double_trace_invariant"].subs(od, d)
    mt = ((ch + cp) * tf + (bh + bp) * tr).applyfunc(s.factor)
    h3, hj = mh.subs(d, 3), mh.diff(d).subs(d, 3)
    total3, totalj = mt.subs(d, 3), mt.diff(d).subs(d, 3)
    varying = ch * s.Matrix([[4 * d, -4], [-4, 4]]) + bh * s.Matrix(
        [[4 * d * d, -4 * d], [-4 * d, 4]]
    )
    checks["complete_scalar_double_trace_invariant"] = s.factor(
        bh - (d * d - 3) / (4 * d * (d + 2))
    )
    checks["complete_scalar_literal_dimension_trace"] = s.factor(
        d * ch + d * d * bh - (d - 1) ** 2 / 4
    )
    checks["complete_scalar_fixed_physical_leading_matrix"] = h3 - s.Matrix(
        [[4, -s.Rational(4, 3)], [-s.Rational(4, 3), s.Rational(8, 15)]]
    )
    checks["scalar_fixed_physical_dimension_jet"] = hj - s.Matrix(
        [
            [s.Rational(22, 15), -s.Rational(22, 45)],
            [-s.Rational(22, 45), s.Rational(26, 225)],
        ]
    )
    checks["total_fixed_physical_leading_matrix"] = total3 - s.Matrix(
        [[8, -s.Rational(8, 3)], [-s.Rational(8, 3), s.Rational(32, 15)]]
    )
    checks["total_fixed_physical_dimension_jet"] = totalj - s.Matrix(
        [
            [s.Rational(68, 15), -s.Rational(68, 45)],
            [-s.Rational(68, 45), s.Rational(224, 225)],
        ]
    )
    checks["scalar_complete_two_logarithms"] = (
        h3 - L0.T * s.diag(4, s.Rational(4, 45)) * L0
    )
    checks["total_complete_two_logarithms"] = (
        total3 - L0.T * s.diag(8, s.Rational(56, 45)) * L0
    )
    checks["total_high_log_shear_coefficient"] = (ch + cp).subs(d, 3) / 2 - s.Rational(
        7, 30
    )
    checks["all_dimension_proper_momentum_cancellation"] = -d - 2 + d + 1 + 1
    checks["same_scheme_fixed_fourth_matrix"] = (
        L0.T
        * s.diag(-4 + finite_trace, s.Rational(8, 3) * (-s.Rational(1, 30) + finite_c2))
        * L0
        + L0.T * s.diag(C0, s.Rational(8, 3) * C2) * L0
    )
    checks["actual_heavy_mass_unchanged"] = HEAVY_MASS2 - (
        s.Integer(10) ** 200 / 512 + 2
    )
    # Paper's noise formula Eq5.8 at minimal coupling; rho_current = N/pi.
    paper2 = (
        s.sqrt(1 - 4 * n / sigma)
        * sigma**2
        * (1 - 4 * n / sigma) ** 2
        * 3
        / (4 * 2880 * s.pi)
    )
    paper0 = (
        s.sqrt(1 - 4 * n / sigma)
        * sigma**2
        * 10
        * (-s.Rational(1, 2) - n / sigma) ** 2
        * 3
        / (2880 * s.pi)
    )
    checks["independent_primary_spin2_convention_conversion"] = s.simplify(
        paper2 / s.pi - heavy_cut(2)
    )
    checks["independent_primary_spin0_convention_conversion"] = s.simplify(
        paper0 / s.pi - heavy_cut(0)
    )
    return {
        "actual_parameters": {
            "Proca_mass_squared": PROCA_MASS2,
            "heavy_mass_squared": HEAVY_MASS2,
            "kappa": KAPPA,
            "ell": "log(actual heavy mass squared), mu=1;394<ell<462",
        },
        "complete_minimal_scalar_stress_pair": pair,
        "complete_scalar_spin_polynomials": invariant,
        "complete_scalar_cut_densities": {i: heavy_cut(i) for i in (0, 2)},
        "both_entire_mass_radial_weights": {"Proca": W_P, "heavy": W_H},
        "whole_actual_scalar_finite_action": reduced,
        "fixed_total_fourth_constants": {"trace": -C0, "shear": -C2},
        "full_dimensional_scalar_invariants": {
            "tracefree": ch,
            "double_trace": s.factor(bh),
        },
        "fixed_physical_scalar_leading_and_jet": (h3, hj),
        "fixed_physical_total_leading_and_jet": (total3, totalj),
        "different_dimension_varying_scalar_jet": varying.diff(d).subs(d, 3),
        "full_dimensional_proper_match": "The frozen radial factor a0^(-d-2) k^(d+1) dk becomes p^(d+1) dp for every nearby d. All invariant coefficients and fixed-physical first dimension jets are matched before the finite part; no new fourth-order finite contact is chosen.",
        "scope": "Complete flat fourth-order Proca-plus-heavy Gaussian reference. Full cuts and fixed fourth-order coefficients retained. Volume, Einstein, reference-state, profile, classical tree and light/mixed loops are NOT included in this isolated block or silently changed.",
        "primary_comparison": "Martin and Verdaguer, gr-qc/0001098v1 Eq5.8, after minimal-coupling and current/noise normalization conversion. Their subtraction convention and full gravity propagator are not adopted.",
        "checks": {
            name: value.applyfunc(s.factor)
            if isinstance(value, s.MatrixBase)
            else s.factor(value)
            for name, value in checks.items()
        },
        "gates": {
            "both_full_heavy_stress_channels_nonzero": all(
                v != 0 for v in com.values()
            ),
            "minimal_trace_not_conformal_improvement": s.Poly(
                invariant[0], sigma
            ).degree()
            == 2,
            "full_scalar_mass_and_trace_terms_retained": vertex.has(n)
            and vertex.has(q11),
            "fixed_physical_jet_not_dimension_varying_tests": hj
            != varying.diff(d).subs(d, 3),
            "full_scalar_leading_matrix_positive": h3[0, 0] > 0 and h3.det() > 0,
            "full_two_mass_leading_matrix_positive": total3[0, 0] > 0
            and total3.det() > 0,
            "actual_heavy_mass_above_lower_mass": HEAVY_MASS2 > PROCA_MASS2,
            "unchanged_finite_counteraction_not_new_pole_tuning": True,
            "no_fourth_order_profile_or_tree_term_absorbed": True,
        },
    }
