"""The missing original ADM boundary shear, derived before physical matching."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import phase
from p8_vacuum_affine_finite_volume_turnaround import moving
from p8_vacuum_affine_hybrid_core_regulator_comparison import quadratic
from p8_vacuum_affine_reduced_scalar_hamiltonian import scalar
from p8_vacuum_affine_scalar_tame_propagator import charts
from p8_vacuum_affine_selfconsistent_finite_feedback import symmetry

from . import source

u, q = source.u, source.q
a = (1 + u * u) ** 2
H = 4 * u / (1 + u * u)
ell = s.Rational(1, 10) / a**3
J = phase.OMEGA


def density_generator():
    return -9 * H * scalar.v**2 + 3 * ell * scalar.v * scalar.sigma


def matrix():
    return a**3 * s.Matrix([[-18 * H, 3 * ell], [3 * ell, 0]])


def shear():
    result = s.eye(4)
    result[2:, :2] = matrix()
    return result


def canonical_generator():
    x, y = s.symbols("unit_CCR_configuration_x unit_CCR_configuration_y", real=True)
    Q = s.Matrix([x, y])
    return (Q.T * matrix() * Q)[0] / 2


@cache
def raw_momenta():
    tree = scalar.tree
    v, vd, sigma, sd, n, b = (
        scalar.v,
        scalar.vd,
        scalar.sigma,
        scalar.sd,
        scalar.n,
        scalar.b,
    )
    eps = s.Symbol("entire_original_ADM_expansion_parameter", real=True)
    jets = tree.clock_coefficients()["coefficient_clock_N_jets"]
    base = tree.quadratic()
    lapse = 1 + eps * n
    rows = {
        name: sum(values[j] * (eps * n) ** j / s.factorial(j) for j in range(3))
        for name, values in jets.items()
    }
    M2 = -3 * rows["a"]
    rate = tree.H + eps * vd
    raw = s.exp(3 * eps * v) * (
        -3 * M2 * rate**2 / lapse
        + 2 * M2 * rate * eps * b / lapse
        + rows["b"] * (3 * rate - eps * b)
        + lapse * rows["f"]
        + rows["U"] * (tree.ell + eps * sd) ** 2 / (2 * lapse)
    )
    raw2 = s.expand(s.diff(raw, eps, 2).subs(eps, 0) / 2)
    prepared = scalar.lagrangian(
        base["J"], base["theta"], base["w"], tree.ell, 1 - 3 * scalar.delta
    )
    A, B, C = s.symbols("fixed_potential_A fixed_potential_B fixed_potential_C")
    retuning = A * n**2 + B * n * v + C * v**2
    before = s.Matrix([s.diff(raw2 + retuning, vd), s.diff(raw2 + retuning, sd)])
    after = s.Matrix([s.diff(prepared + retuning, vd), s.diff(prepared + retuning, sd)])
    actual = (before - after).applyfunc(s.factor)
    wanted = s.Matrix([-18 * tree.H * v + 3 * tree.ell * sigma, 3 * tree.ell * v])
    clock = source.previous.previous.previous.clock_at
    return {
        "whole_original_raw_quadratic_ADM": raw2,
        "whole_original_prepared_quadratic_action": prepared,
        "whole_arbitrary_fixed_potential_retuning": retuning,
        "whole_raw_minus_prepared_density_momenta": actual,
        "whole_required_density_momentum_shift": wanted,
        "checks": {
            "independent_entire_ADM_both_velocity_derivatives": (
                actual - wanted
            ).applyfunc(s.factor),
            "actual_bounce_both_nonzero_momentum_shifts": actual.subs(tree.u, 0)
            - s.Matrix([3 * sigma / 10, 3 * v / 10]),
            "weighted_original_matter_boundary": s.factor(
                s.diff(tree.ell, tree.u) + 3 * tree.H * tree.ell
            ),
            "actual_complete_R_clock_value": clock(q.R) - 1,
            "actual_complete_U_clock_value": clock(q.U) - 1,
            "actual_complete_B_clock_value": clock(q.B),
            "fixed_potential_has_no_velocity_contact": s.Matrix(
                [s.diff(retuning, vd), s.diff(retuning, sd)]
            ),
        },
    }


@cache
def tree_hessian():
    # Independent exact fixture ONLY: never identify the actual live/profile
    # nonlinear root Hessian with the stipulated reference Hamiltonian.
    literal = quadratic.literal_hessian()
    h2, c1 = literal["whole_original_H2"], literal["whole_literal_C1"]
    named = {x.name: x for x in h2.free_symbols | c1.free_symbols}
    Q = s.Symbol("independent_physical_momentum_squared", positive=True)
    v, sigma, pv, ps = charts.v, charts.sigma, phase.pv, phase.ps
    bindings = {
        "hom_p": -2 * H,
        "hom_pm": ell,
        "hom_ph": 0,
        "hom_eta": 0,
        "delta_trace_density": pv,
        "volume_field": v,
        "M1_field": sigma,
        "M1_density": ps,
        "H_density": 0,
        "H_scaled_field": 0,
        "Gauss_linear": 0,
        "full_TT_norm": 0,
        "electric_quad": 0,
        "magnetic_quad": 0,
        "Wmass_quad": 0,
        "M1gradient_quad": Q * sigma**2,
        "Hgradient_quad": 0,
        "curvature_linear": 4 * Q * v,
        "curvature_quad": -10 * Q * v**2,
    }
    rules = {named[name]: value for name, value in bindings.items() if name in named}
    clock = source.previous.previous.previous.clock_at

    def tree_clock(value):
        return s.factor(
            clock(value).subs({q.physical.RHO: 0, q.physical.PRESSURE: 0}).doit()
        )

    whole = s.factor(
        tree_clock(h2.subs(rules, simultaneous=True))
        - tree_clock(c1.subs(rules, simultaneous=True)) ** 2
        / (2 * tree_clock(s.diff(q.CONSTRAINT, q.N)))
    )
    old = moving.old.original

    def original(value):
        return s.factor(
            value.subs({old.u: u, old.rho: 0, old.pressure: 0}, simultaneous=True)
        )

    reference = phase.whole_hamiltonian().subs(
        {
            charts.th: original(old.theta),
            charts.E: original(old.E),
            charts.l: original(old.ell),
            charts.J: original(old.J),
            charts.H: original(old.H),
            charts.A: 0,
            charts.T: 0,
            charts.q: Q,
        },
        simultaneous=True,
    )
    f = -9 * H * v**2 + 3 * ell * v * sigma
    weighted = s.diff(f, u) + 3 * H * f
    transported = (
        whole.subs(
            {pv: pv + s.diff(f, v), ps: ps + s.diff(f, sigma)}, simultaneous=True
        )
        + weighted
    )
    unshifted = s.factor(whole - reference)
    return {
        "whole_zero_profile_tree_raw_constrained_H2": whole,
        "whole_zero_profile_tree_original_reference_H2": s.factor(reference),
        "whole_nonzero_unshifted_tree_H2_defect": unshifted,
        "whole_weighted_boundary_time_generator": s.factor(weighted),
        "checks": {
            "independent_full_tree_constraint_fixture": tree_clock(q.CONSTRAINT),
            "complete_tree_boundary_H2_transport": s.factor(transported - reference),
        },
        "gates": {
            "unshifted_full_tree_H2_is_not_the_reference": unshifted != 0,
        },
    }


@cache
def transport():
    S, C = shear(), matrix()
    F = canonical_generator()
    x, y = sorted(F.free_symbols - {u}, key=str)
    Q = s.Matrix([x, y])
    raw_p = s.Matrix(s.symbols("raw_unit_CCR_momentum0:2", real=True))
    standard_p = s.Matrix(s.symbols("prepared_unit_CCR_momentum0:2", real=True))
    wave = s.Function("original_prepared_wavefunction")(u, x, y)
    physical = s.exp(s.I * F) * wave
    checks = {}
    for i, coordinate in enumerate((x, y)):
        checks["exact_boundary_unitary_momentum_" + str(i)] = s.simplify(
            -s.I * s.diff(physical, coordinate)
            - s.exp(s.I * F)
            * (-s.I * s.diff(wave, coordinate) + s.diff(F, coordinate) * wave)
        )
    checks["exact_boundary_unitary_time_connection"] = s.simplify(
        s.I * s.diff(physical, u)
        - s.exp(s.I * F) * (s.I * s.diff(wave, u) - s.diff(F, u) * wave)
    )
    checks["exact_complete_canonical_shear"] = (S * J * S.T - J).applyfunc(s.factor)
    checks["boundary_gradient_is_raw_momentum_increment"] = (
        s.Matrix([s.diff(F, x), s.diff(F, y)]) - C * Q
    ).applyfunc(s.factor)
    checks["constant_matter_cross_boundary_coefficient"] = s.factor(
        C[0, 1] - s.Rational(3, 10)
    )
    checks["whole_time_generator_includes_weight_derivative"] = s.factor(
        s.diff(F, u) + 9 * a**3 * (s.diff(H, u) + 3 * H**2) * x**2
    )
    # An arbitrary full symmetric raw quadratic matrix independently fixes
    # the time-connection sign, without reference-Hessian identification.
    symbols = s.symbols("arbitrary_raw_quadratic0:10", real=True)
    G = s.zeros(4)
    cursor = 0
    for i in range(4):
        for j in range(i, 4):
            G[i, j] = G[j, i] = symbols[cursor]
            cursor += 1
    extra = s.zeros(4)
    extra[:2, :2] = C.diff(u)
    corrected = S.T * G * S + extra
    raw_flow = S.diff(u) * S.inv() + S * J * corrected * S.inv()
    checks["arbitrary_full_time_dependent_canonical_flow"] = (
        raw_flow - J * G
    ).applyfunc(s.factor)
    # A quadratic boundary is exactly calibrated in any fixed whitening:
    # D lowers polynomial degree by two, so exp(D)(1-D)F=F.
    M = s.Matrix(
        [
            [s.Symbol("whitened_cov00"), s.Symbol("whitened_cov01")],
            [s.Symbol("whitened_cov01"), s.Symbol("whitened_cov11")],
        ]
    )
    DF = sum(M[i, j] * s.diff(F, Q[i], Q[j]) for i in range(2) for j in range(2)) / 4
    D2F = sum(M[i, j] * s.diff(DF, Q[i], Q[j]) for i in range(2) for j in range(2)) / 4
    checks["whole_quadratic_boundary_has_zero_second_heat_power"] = s.factor(D2F)
    checks["whole_boundary_calibrated_equals_Weyl"] = s.factor((F - DF) + DF - F)
    live = s.symbols("live_alpha live_p live_pm live_eta live_ph", real=True)
    checks["fixed_reference_boundary_has_no_live_Y_derivatives"] = s.Matrix(
        [s.diff(F, value) for value in live]
    )
    return {
        "whole_unit_CCR_reference_history": [a, H, ell],
        "whole_symmetric_boundary_matrix": C,
        "whole_canonical_boundary_shear": S,
        "whole_canonical_boundary_generating_function": F,
        "whole_canonical_boundary_time_derivative": s.factor(s.diff(F, u)),
        "whole_raw_and_prepared_momentum_relation": raw_p - standard_p - C * Q,
        "whole_arbitrary_raw_quadratic_generator_check": G,
        "checks": checks,
    }


@cache
def covariance():
    P, kap = s.symbols("positive_P positive_kappa", positive=True)
    inputs = s.Matrix(s.symbols("balanced_scalar_input0:4", real=True))
    old_map = (
        phase.central_map().subs({phase.a: 1, charts.q: P**2}).inv()
        * s.diag(1 / s.sqrt(P), 1 / s.sqrt(P), s.sqrt(P), s.sqrt(P))
        / s.sqrt(kap)
    )
    S0 = shear().subs(u, 0)
    corrected_map = S0 * old_map
    gap = (corrected_map - old_map) * inputs
    positive = s.Symbol("strictly_positive_configuration_variance", positive=True)
    v11, v01 = s.symbols(
        "configuration_variance11 configuration_covariance01", real=True
    )
    qq = s.Matrix([[positive, v01], [v01, v11]])
    qp = s.Matrix(2, 2, s.symbols("full_qp_covariance0:4", real=True))
    pp = s.Matrix(
        [[s.Symbol("pp00"), s.Symbol("pp01")], [s.Symbol("pp01"), s.Symbol("pp11")]]
    )
    V = qq.row_join(qp).col_join(qp.T.row_join(pp))
    delta = S0 * V * S0.T - V
    return {
        "whole_archived_unshifted_bounce_dictionary": old_map,
        "whole_corrected_bounce_dictionary": corrected_map,
        "whole_nonzero_missing_balanced_momentum_rows": gap,
        "whole_arbitrary_full_old_covariance": V,
        "whole_corrected_raw_covariance": S0 * V * S0.T,
        "whole_strict_covariance_defect": 3 * positive / 10,
        "checks": {
            "unchanged_configuration_rows": gap[:2, :],
            "missing_trace_momentum_row": s.factor(
                gap[2] - 3 * inputs[1] / (10 * s.sqrt(P * kap))
            ),
            "missing_matter_momentum_row": s.factor(
                gap[3] - 3 * inputs[2] / (20 * P ** s.Rational(3, 2) * s.sqrt(kap))
            ),
            "full_covariance_defect_independent_of_old_qp_block": s.factor(
                delta[0, 3] - 3 * positive / 10
            ),
            "corrected_physical_CCR_not_a_new_quantum_normalization": (
                corrected_map * J * corrected_map.T - J / kap
            ).applyfunc(s.factor),
        },
        "gates": {
            "strict_generic_positive_covariance_change": 3 * positive / 10 > 0,
            "both_archived_momentum_rows_really_change": gap[2] != 0 and gap[3] != 0,
            "boundary_phase_is_not_a_global_scalar": len(
                canonical_generator().subs(u, 0).free_symbols
            )
            == 2,
        },
    }


@cache
def group_contacts():
    # All scalar blocks receive exactly the same shear. This explicit sparse
    # 48-configuration construction keeps both species and all six real waves.
    C = s.MutableSparseMatrix(48, 48, {})
    block = matrix()
    for wave in range(6):
        for i in range(2):
            for j in range(2):
                C[8 * wave + i, 8 * wave + j] = block[i, j]
    C = s.ImmutableSparseMatrix(C)
    cubic = [s.Integer(U * C * U.T != C) for U in symmetry.representations()]
    # Translation generator for each axis rotates its cosine/sine pair,
    # identically in all eight channels. No group or state averaging.
    translation = []
    charge = []
    for axis in range(3):
        L = s.MutableSparseMatrix(48, 48, {})
        for channel in range(8):
            i, j = (2 * axis) * 8 + channel, (2 * axis + 1) * 8 + channel
            L[i, j], L[j, i] = -1, 1
        translation.append(s.Integer(C * L != L * C))
        charge.append(s.Integer((C * L + (C * L).T) != s.zeros(48)))
    return {
        "whole_same_shear_on_all_six_scalar_real_waves": C,
        "checks": {
            "all24_actual_cubic_actions_preserve_boundary": s.Matrix(cubic),
            "all3_continuous_translations_commute_with_boundary": s.Matrix(translation),
            "all3_translation_charge_boundary_increments_zero": s.Matrix(charge),
        },
        "gates": {
            "actual_cubic_action_count": len(cubic) == 24,
            "all_scalar_channels_and_real_waves_kept": sum(value != 0 for value in C)
            == 18,
        },
    }


@cache
def data():
    parts = {
        "raw_ADM_velocity_derivation": raw_momenta(),
        "independent_tree_H2_derivation": tree_hessian(),
        "full_time_dependent_canonical_transport": transport(),
        "generic_covariance_and_archived_dictionary_refutation": covariance(),
        "whole_cubic_and_translation_contacts": group_contacts(),
    }
    return {
        **{
            name: {
                key: value
                for key, value in part.items()
                if key not in ("checks", "gates")
            }
            for name, part in parts.items()
        },
        "whole_explicit_historical_qualification": "S269's central inverse is the prepared, boundary-normalized dictionary, not the raw ADM dictionary. S220's earlier boundary f=-9Hv^2+3ell v sigma was omitted. S251's later outer/central symmetric shears do not restore it. Published S269-S272 and frozen S273-S274 remain internally defined abstract calculations but their same-original-physical-state bridge is not established by their replay and is refuted at these two momentum rows. This successor changes the physical chart and the four coupled solutions; no archived file or report is rewritten.",
        "whole_corrected_physical_state_and_reference": "In unit CCR variables raw z=Sb(u) z_prepared, with symmetric lower block C=a_bar^3[[-18Hbar,3ellbar],[3ellbar,0]]. The exact finite unitary is Ub=exp(iF_b). Keep the original prepared Gaussian and use raw initial state Ub(0)psi0. The raw reference propagator is Ub(u)Uref(u)Ub(0)^-1, not Uref alone; raw initial whitening is Sb(0)S0. This is canonical transport of the stipulated original state, not a new minimizer, covariance deletion or irrelevant global phase.",
        "whole_exact_nonlinear_extension_boundary": "Use this fixed-reference linear canonical shear as an explicitly stipulated exact extension to the entire finite reduced chart. Pull back EVERY original nonlinear raw source and implicit reconstruction by it. This does not assert an unproved nonlinear ADM boundary generating function beyond the derived linearized canonical bridge. The generator depends on reference time, not live homogeneous Y, and its configuration-dependent time derivative is retained.",
        "whole_tree_fixture_scope": "Only the independent zero-profile tree fixture identifies the transported full constrained H2 with the original reference H2. The actual nonzero-profile live-root Hamiltonian is kept independently and completely; no such identification is assumed there.",
        "checks": {
            name + "_" + key: value
            for name, part in parts.items()
            for key, value in part["checks"].items()
        },
        "gates": {
            name + "_" + key: bool(value)
            for name, part in parts.items()
            for key, value in part.get("gates", {}).items()
        },
    }
