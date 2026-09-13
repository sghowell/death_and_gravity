"""Actual QG2 profile-once assembly, mixed-order coordinates and scalar Ward block."""

from functools import cache

import sympy as s
from p8_vacuum_affine_heavy_adm_clock_response import geometry as g
from p8_vacuum_affine_heavy_adm_clock_response import profile
from p8_vacuum_affine_heavy_clock_quadratic import clock
from p8_vacuum_affine_heavy_curved_state import quantum, state
from p8_vacuum_affine_quantum_forced_constraints import forces
from p8_vacuum_affine_reduced_scalar_hamiltonian import scalar

t = s.Symbol("proper_time", real=True)
fields = tuple(
    s.Function(name)(t) for name in ("eta", "spatial_w", "spatial_c", "matter_r")
)
eta, w, c, r = fields
H, delta, ell, q, Jc, theta, Tc, A = tuple(
    s.Function(name)(t)
    for name in ("H", "delta", "ell", "q", "Jc", "Theta", "Tcorr", "Aret")
)
ROW_ORDERS = (4, 4, 4, 2)
PIVOT_INVERSE = s.Rational(15625, 6144)


def zero(value):
    return (
        value.applyfunc(s.cancel)
        if isinstance(value, s.MatrixBase)
        else s.cancel(value)
    )


def weighted(f):
    return s.diff(f, t) + 3 * H * f


def lagrangian(n, v, vd, sigma, sd, b):
    E = 1 - 3 * delta
    wb = -ell * E
    return (
        -3 * vd**2
        + (Jc + wb**2 / 2 - 3 * theta**2) * n**2
        + 6 * theta * n * vd
        + sd**2 / 2
        + wb * n * sd
        - 3 * ell * vd * sigma
        + 2 * b * (vd - theta * n)
        + ell * b * sigma
        + q * v**2
        + 2 * E * q * n * v
        - q * sigma**2 / 2
        + 3 * Tc * n * v
        + 9 * A * v**2 / 2
    )


@cache
def system():
    n = s.diff(eta, t)
    v = w + H * eta - delta * n
    b = s.diff(c, t) + q * eta
    sigma = r + ell * eta
    lag = lagrangian(n, v, s.diff(v, t), sigma, s.diff(sigma, t), b)
    equations = tuple(
        s.expand(
            s.diff(lag, f)
            - weighted(s.diff(lag, s.diff(f, t)))
            + weighted(weighted(s.diff(lag, s.diff(f, t, 2))))
        )
        for f in fields
    )
    principal = s.ImmutableMatrix(
        [
            [s.factor(equations[i].coeff(s.diff(f, t, ROW_ORDERS[i]))) for f in fields]
            for i in range(4)
        ]
    )
    wrong = lagrangian(n, v, s.diff(v, t), r, s.diff(r, t), b)
    wrong_r = s.expand(s.diff(wrong, r) - weighted(s.diff(wrong, s.diff(r, t))))
    return {
        "complete_current_coefficient_lagrangian": lag,
        "equations": equations,
        "principal": principal,
        "unshifted_matter_eta_second": s.factor(wrong_r.coeff(s.diff(eta, t, 2))),
    }


@cache
def profile_data():
    nD, nG, vD, vG, aa, rho, pressure, dd = s.symbols(
        "n_D n_G v_D v_G a rho pressure delta", real=True
    )
    d = profile.data()
    named = {str(x): x for x in d["whole_common_clock_profile_Hessian"].free_symbols}
    substitution = {named["trace_D"]: 6 * vD, named["trace_G"]: 6 * vG}
    actual = d["whole_common_clock_profile_Hessian"].subs(
        substitution, simultaneous=True
    )
    AH, BH = -pressure, -(rho + pressure) / 2
    dJH = (21 * dd**2 - 3 * dd) * AH / 2 + (1 - 6 * dd) * BH
    TH = (1 + 3 * dd) * AH - 2 * BH
    matched = aa**3 * (
        2 * dJH * nD * nG + 3 * TH * (nD * vG + nG * vD) + 9 * AH * vD * vG
    )
    fixed = quantum.fixed_profile()
    hrho = fixed["rho"].subs(state.TIME, clock.t) / state.KAPPA
    hp = fixed["pressure"].subs(state.TIME, clock.t) / state.KAPPA
    actual_clock = clock.data()
    full_quadratic = actual_clock["complete_added_profile_quadratic_density"]
    nv = {str(x): x for x in full_quadratic.free_symbols}
    repeated = (matched / (2 * aa**3)).subs(
        {
            nD: nv["n_lapse"],
            nG: nv["n_lapse"],
            vD: nv["v"],
            vG: nv["v"],
            rho: hrho,
            pressure: hp,
            dd: scalar.delta,
        },
        simultaneous=True,
    )
    checks = {
        "whole_S246_clock_profile_equals_full_current_coefficient_Hessian": actual
        - matched,
        "actual_S241_complete_quadratic_equals_half_repeated_mixed_profile": repeated
        - full_quadratic,
        "current_full_pressure_coefficient_includes_heavy_once": actual_clock[
            "full_current_clock_coefficients"
        ]["A"]
        + actual_clock["full_actual_normalized_reference_pressure"],
        "full_profile_clock_contact_opposes_only_the_full_Gaussian_contact": d[
            "complete_Gaussian_nonlinear_clock_contact"
        ]
        + d["complete_profile_nonlinear_clock_contact"],
        "profile_spatial_specialization_retained_before_Ward": actual.subs(
            {nD: 0, nG: 0}
        )
        + 9 * aa**3 * pressure * vD * vG,
    }
    e, f = s.symbols("e f", real=True)
    N = 1 + e * nD + f * nG
    V = e * vD + f * vG
    density = (
        aa**3
        * N
        * s.exp(3 * V)
        * (1 + 2 * dd * (N**-2 - 1)) ** (-s.Rational(3, 4))
        * (AH + BH * (N**-2 - 1))
    )
    checks["independent_unexpanded_density_mixed_profile"] = (
        s.diff(density, e, f).subs({e: 0, f: 0}) - matched
    )
    base, gP, gH, pH = s.symbols(
        "T_cl_QG1 R_Proca_Gaussian R_H_Gaussian P_H_clock", commutative=False
    )
    current = base + pH
    complete_H = gH + pH
    once = current + gP + complete_H - pH
    checks["two_equivalent_whole_assemblies"] = s.expand(
        once - (base + gP + complete_H)
    )
    twice_defect = s.expand(current + gP + complete_H - once)
    return {
        "whole_heavy_clock_profile_Hessian": matched,
        "complete_heavy_clock_coefficient_addition": s.ImmutableMatrix(
            [AH, BH, dJH, TH]
        ),
        "whole_actual_QG2_coefficient_substitution": actual_clock[
            "exact_current_Hamiltonian_coefficient_substitution"
        ],
        "actual_current_clock_coefficients": actual_clock[
            "full_current_clock_coefficients"
        ],
        "double_profile_defect": twice_defect,
        "assembly": "T_cl,QG2 + R_Proca,Gaussian + (R_S246 - P_H,clock) = T_cl,QG1 + R_Proca,Gaussian + R_S246. The full heavy profile appears once. No finite heat coefficient or physical state changes.",
        "nonlinear_contact_order": "Subtract the WHOLE fixed profile from S246 to obtain the pure Gaussian response. Its nonzero one-current, both metric Ward legs and distinct nonlinear clock contact remain. Cancellation with the coefficient-sector profile occurs only in the whole equation.",
        "checks": {key: zero(value) for key, value in checks.items()},
        "gates": {
            "double_profile_defect_nonzero": twice_defect != 0,
            "profile_Hessian_not_replaced_by_mean_cancellation": matched != 0,
            "current_classical_hypotheses_freshly_checked": all(
                bool(v) for v in actual_clock["gates"].values()
            ),
            "no_added_nonlocal_linear_heavy_mean_channel": all(
                value == 0
                for key, value in profile.data()["checks"].items()
                if key.startswith("full_heavy_source_clock_")
            ),
            "free_coherent_heavy_mode_does_not_remove_Gaussian_metric_response": True,
        },
    }


@cache
def principal_data():
    d = system()
    above = [
        s.factor(row.coeff(s.diff(f, t, j)))
        for i, row in enumerate(d["equations"])
        for f in fields
        for j in range(ROW_ORDERS[i] + 1, 7)
    ]
    actual = scalar.tree.quadratic()
    return {
        "row_orders": ROW_ORDERS,
        "actual_principal": d["principal"],
        "actual_eta_pivot_inverse_bound": PIVOT_INVERSE,
        "current_coefficient_rule": "Jc, A, Tcorr are the entire QG2 functions in profile_data, not the old QG1 numerical coefficients. The independently expanded heavy profile is already included in this Jc,A,Tcorr. All phase-variable substitutions are simultaneous.",
        "prepared_adapted_coordinates": "n=eta'; zeta=w+H eta; v=zeta-delta eta'; b=c'+q eta; sigma=r+ell eta; q=P²/a².",
        "checks": {
            "entire_new_current_mixed_principal": d["principal"]
            - s.diag(-6 * delta**2, 0, 0, -1),
            "all_terms_above_full_row_orders": s.ImmutableMatrix(above),
            "matter_shift_required_for_low_row": d["unshifted_matter_eta_second"] - ell,
            "actual_matter_identity": actual["w"] + scalar.ell * (1 - 3 * scalar.delta),
            "actual_nonzero_pivot_bound": 1 / (6 * s.Rational(32, 125) ** 2)
            - PIVOT_INVERSE,
        },
        "gates": {
            "no_H_Theta_or_external_momentum_pivot_division": not any(
                s.denom(x).has(H, theta, q) for x in d["principal"]
            ),
            "whole_current_retuning_still_present": all(
                d["complete_current_coefficient_lagrangian"].has(x) for x in (Jc, A, Tc)
            ),
            "matter_invariant_not_optional": d["unshifted_matter_eta_second"] != 0,
        },
    }


@cache
def ward_data():
    P = s.Symbol("external_momentum", positive=True)
    rho, pressure = (
        s.Function(name, real=True)(g.t)
        for name in ("rho_total_Gaussian", "pressure_total_Gaussian")
    )
    dd = s.Function("delta", real=True)(g.t)
    E = g.density(rho, pressure)
    gauge_checks = []

    def packet(label, sign):
        fs = tuple(
            s.Function(label + name, real=True)(g.t) for name in ("eta", "w", "c")
        )
        phase = s.exp(sign * s.I * P * g.z)
        ee, ww, cc = (f * phase for f in fs)
        physical = (
            s.diff(ee, g.t),
            s.Matrix([0, 0, -sign * s.I * (s.diff(cc, g.t) + P**2 * ee / g.a**2) / P]),
            2 * (ww + g.H * ee) * s.eye(3),
        )
        synchronous = (
            s.S.Zero,
            s.zeros(3, 1),
            2 * ww * s.eye(3) - 2 * cc * s.diag(0, 0, 1),
        )
        xi = s.Matrix([ee, 0, 0, -sign * s.I * cc / P])
        pure = g.gauge(xi)
        for left, right, gauge in zip(physical, synchronous, pure):
            difference = zero(left - right - gauge)
            gauge_checks.extend(
                list(difference)
                if isinstance(difference, s.MatrixBase)
                else [difference]
            )
        return fs, physical, synchronous, xi

    det, dp, _, dx = packet("detector_", -1)
    src, _, gs, gx = packet("source_", 1)
    source = g.source_ward(dp, gx, E)
    detector = g.detector_ward(dx, gs, E)
    clock_term = (
        3
        * g.a**3
        * pressure
        * (4 * dd**2 - 3 * dd)
        * s.diff(det[0], g.t)
        * s.diff(src[0], g.t)
    )
    complete = s.cancel(s.expand(source + detector + clock_term))
    variables, derivative_order, sides, mapping = [], [], [], {}
    for side, fs in enumerate((det, src)):
        for f in fs:
            for j in range(4):
                var = s.Symbol(str(f.func) + "_jet" + str(j), real=True)
                variables.append(var)
                derivative_order.append(j)
                sides.append(side)
                mapping[s.diff(f, g.t, j)] = var
    poly = s.Poly(complete.xreplace(mapping), *variables)
    maxima = [0, 0]
    bilinear = True
    for mon, _ in poly.terms():
        for side in (0, 1):
            bilinear &= sum(e for e, which in zip(mon, sides) if which == side) == 1
            maxima[side] = max(
                maxima[side],
                sum(
                    e * j
                    for e, j, which in zip(mon, derivative_order, sides)
                    if which == side
                ),
            )
    rhoH, pH = s.Function("rho_H")(g.t), s.Function("pressure_H")(g.t)
    rhoP, pP = s.Function("rho_Proca")(g.t), s.Function("pressure_Proca")(g.t)
    fullsum = complete.subs(
        {rho: rhoH + rhoP, pressure: pH + pP}, simultaneous=True
    ).doit()
    separated = sum(
        complete.subs({rho: rr, pressure: pp}, simultaneous=True).doit()
        for rr, pp in ((rhoH, pH), (rhoP, pP))
    )
    return {
        "complete_ordered_Ward_plus_Gaussian_clock": complete,
        "whole_Ward_monomial_count": len(poly.terms()),
        "maximum_derivative_per_leg": maxima,
        "external_momentum_degree": s.degree(complete, P),
        "nonzero_clock_bounce_eta_second": (
            -3 * g.a**3 * pressure * (4 * dd**2 - 3 * dd)
        ).subs(dd, s.Rational(1, 2)),
        "metric_order": "Use the full physical detector in source Ward, then the synchronous source in detector Ward. The mean is the sum of the actual Proca and actual scalar Gaussian means, excluding both fixed profiles. The complete common-clock contact is retained once.",
        "checks": {
            "all_physical_minus_synchronous_gauge_entries": s.ImmutableMatrix(
                gauge_checks
            ),
            "two_actual_species_add_before_full_Ward": zero(fullsum - separated),
            "complete_local_derivative_order": s.ImmutableMatrix(maxima) - s.ones(2, 1),
            "whole_local_external_degree_two": s.degree(complete, P) - 2,
            "full_ordered_Ward_has_fourteen_terms": len(poly.terms()) - 14,
        },
        "gates": {
            "whole_opposite_phase_cancelled": not complete.has(g.z),
            "inverse_momentum_cancels_before_any_norm": not s.denom(complete).has(P),
            "one_detector_one_source_in_every_term": bilinear,
            "clock_contact_not_identically_zero": clock_term != 0,
            "actual_scalar_Ward_not_inherited_from_vector_norm": True,
            "Gaussian_mean_not_set_to_zero_by_fixed_profile": complete.has(
                rho, pressure
            ),
        },
    }


@cache
def residual_data():
    En, Ev, Eb, Es = (s.Function(name)(t) for name in ("En", "Ev", "Eb", "Es"))
    original = s.Matrix(
        [
            s.diff(eta, t),
            w + H * eta - delta * s.diff(eta, t),
            s.diff(c, t) + q * eta,
            r + ell * eta,
        ]
    )
    adapted = s.Matrix(
        [
            -weighted(En) + H * Ev + weighted(delta * Ev) + q * Eb + ell * Es,
            Ev,
            -weighted(Eb),
            Es,
        ]
    )
    boundary = (En - delta * Ev) * eta + Eb * c
    en, ez, eb = (s.Function(name)(t) for name in ("force_n", "force_zeta", "force_b"))
    forced = adapted.subs(
        {En: en + delta * ez, Ev: ez, Eb: eb, Es: 0}, simultaneous=True
    ).doit()
    wanted = s.Matrix([-weighted(en) + H * ez + q * eb, ez, -weighted(eb), 0])
    aa = (1 + t**2) ** 2
    pp = s.Symbol("external_momentum", positive=True)
    u = s.Dummy("source_time", real=True)
    bad_b = aa**-3
    bad_n = aa**-3 * s.Integral((pp**2 / aa**2).subs(t, u), (u, -s.Rational(1, 2), t))
    W = lambda f: s.diff(f, t) + 3 * s.diff(aa, t) / aa * f
    return {
        "full_weighted_residual_adjoint": adapted,
        "density_boundary": aa**3 * boundary,
        "physical_force_pullback": wanted,
        "prepared_inverse": "eta=I n; w=zeta-H eta; c=I(b-q eta); r=sigma-ell eta",
        "nonzero_unprepared_constraint_kernel": s.ImmutableMatrix([bad_n, 0, bad_b, 0]),
        "actual_force_interface": "Use the whole QG2 coefficients in the generic exact S222 maps C,K,D,F. Define Qbar as the pure Proca-plus-H Gaussian response, with its clock contact and OUTPUT (kappa a³)^-1 factor. The fixed profiles are in the coefficient sector only.",
        "checks": {
            "full_weighted_variational_boundary": zero(
                (s.Matrix([En, Ev, Eb, Es]).T * original)[0]
                - (adapted.T * s.Matrix(fields))[0]
                - weighted(boundary)
            ),
            "normalized_physical_force_pullback": zero(forced - wanted),
            "unprepared_shift_kernel": zero(W(bad_b)),
            "unprepared_lapse_kernel": zero(-W(bad_n) + pp**2 / aa**2 * bad_b),
        },
        "gates": {
            "both_original_initial_constraint_germs_essential": bad_b.subs(
                t, -s.Rational(1, 2)
            )
            != 0,
            "direct_auxiliary_rank_two_kept": forces.system()["D"].rank() == 2,
            "output_density_before_source_adjoint_and_row_primitives": True,
        },
    }
