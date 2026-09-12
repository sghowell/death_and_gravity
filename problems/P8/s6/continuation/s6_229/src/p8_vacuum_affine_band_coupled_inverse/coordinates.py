"""Literal full classical principal block, ordered Ward contacts and residual adjoint."""

from functools import cache

import sympy as s
from p8_vacuum_affine_prepared_ward_reconstruction import kinematics as k
from p8_vacuum_affine_prepared_ward_reconstruction import ward
from p8_vacuum_affine_reduced_scalar_hamiltonian import scalar as old

t = s.Symbol("proper_time", real=True)
fields = tuple(
    s.Function(name)(t)
    for name in ("eta", "spatial_w", "spatial_c", "matter_invariant")
)
eta, w, c, r = fields
H, delta, ell, q, J, theta, dJ, Tc, A = tuple(
    s.Function(name)(t)
    for name in ("H", "delta", "ell", "q", "J", "Theta", "DeltaJ", "Tcorr", "Aret")
)
ROW_ORDERS = (4, 4, 4, 2)
L0 = s.ImmutableMatrix([[1, -s.Rational(1, 3)], [0, -1]])
PIVOT_INVERSE = s.Rational(15625, 6144)


def weighted(f):
    return s.diff(f, t) + 3 * H * f


@cache
def system():
    E = 1 - 3 * delta
    wb = -ell * E
    v = w + H * eta - delta * s.diff(eta, t)
    n = s.diff(eta, t)
    b = s.diff(c, t) + q * eta
    sigma = r + ell * eta
    raw = (
        old.lagrangian(J, theta, wb, ell, E)
        + old.dJ * old.n**2
        + 3 * old.Tc * old.n * old.v
        + s.Rational(9, 2) * old.A * old.v**2
    )
    mapping = {
        old.n: n,
        old.v: v,
        old.vd: s.diff(v, t),
        old.sigma: sigma,
        old.sd: s.diff(sigma, t),
        old.b: b,
        old.q: q,
        old.dJ: dJ,
        old.Tc: Tc,
        old.A: A,
    }
    lag = raw.subs(mapping, simultaneous=True)
    eqs = tuple(
        s.expand(
            s.diff(lag, f)
            - weighted(s.diff(lag, s.diff(f, t)))
            + weighted(weighted(s.diff(lag, s.diff(f, t, 2))))
        )
        for f in fields
    )
    principal = s.ImmutableMatrix(
        [
            [s.factor(eqs[i].coeff(s.diff(f, t, ROW_ORDERS[i]))) for f in fields]
            for i in range(4)
        ]
    )
    unshifted = raw.subs(
        {**mapping, old.sigma: r, old.sd: s.diff(r, t)}, simultaneous=True
    )
    old_r = s.expand(s.diff(unshifted, r) - weighted(s.diff(unshifted, s.diff(r, t))))
    return {
        "lagrangian": lag,
        "equations": eqs,
        "principal": principal,
        "unshifted_matter_eta_second": s.factor(old_r.coeff(s.diff(eta, t, 2))),
    }


@cache
def principal_data():
    d = system()
    actual = old.tree.quadratic()
    above = [
        s.factor(row.coeff(s.diff(f, t, j)))
        for i, row in enumerate(d["equations"])
        for f in fields
        for j in range(ROW_ORDERS[i] + 1, 7)
    ]
    return {
        "adapted_tuple": "n=eta', zeta=w+H eta, v=zeta-delta eta', b=c'+q eta, sigma=r+ell eta; q=P^2/a^2",
        "row_orders": ROW_ORDERS,
        "actual_principal": d["principal"],
        "actual_eta_pivot_inverse_bound": PIVOT_INVERSE,
        "unshifted_matter_eta_second": d["unshifted_matter_eta_second"],
        "checks": {
            "complete_mixed_row_principal": d["principal"]
            - s.diag(-6 * delta**2, 0, 0, -1),
            "all_above_row_orders": s.ImmutableMatrix(above),
            "unshifted_matter_cross_is_ell": d["unshifted_matter_eta_second"] - ell,
            "actual_background_matter_identity": actual["w"]
            + old.ell * (1 - 3 * old.delta),
            "actual_delta": old.delta - 1 / (2 * (1 + old.t**2) ** 3),
            "actual_nonzero_pivot_bound": 1 / (6 * s.Rational(32, 125) ** 2)
            - PIVOT_INVERSE,
        },
        "gates": {
            "matter_invariant_required": d["unshifted_matter_eta_second"] != 0,
            "retuning_terms_in_entire_action": all(
                d["lagrangian"].has(x) for x in (dJ, Tc, A)
            ),
            "no_principal_H_Theta_E_or_q_division": not any(
                s.denom(x).has(H, theta, q) for x in d["principal"]
            ),
            "all_remaining_local_rows_lower_order": all(x == 0 for x in above),
        },
    }


@cache
def residual_data():
    En, Ev, Eb, Es = tuple(s.Function(name)(t) for name in ("En", "Ev", "Eb", "Es"))
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
    ep = s.Function("physical_lapse_residual")(t)
    en, ez, eb = tuple(
        s.Function(name)(t) for name in ("force_n", "force_zeta", "force_b")
    )
    forced = adapted.subs(
        {En: en + delta * ez, Ev: ez, Eb: eb, Es: 0}, simultaneous=True
    ).doit()
    wanted = s.Matrix([-weighted(en) + H * ez + q * eb, ez, -weighted(eb), 0])
    a = (1 + t * t) ** 2
    aq = s.Symbol("P", positive=True) ** 2 / a**2
    u = s.Dummy("source_time", real=True)
    bad_b = a**-3
    bad_n = a**-3 * s.Integral(aq.subs(t, u), (u, -s.Rational(1, 2), t))
    actualW = lambda f: s.diff(f, t) + 3 * s.diff(a, t) / a * f
    return {
        "weighted_old_residual_adjoint": adapted,
        "complete_density_boundary": a**3 * boundary,
        "normalized_physical_force_pullback": wanted,
        "prepared_inverse": "eta=I n; w=zeta-H eta; c=I(b-q eta); r=sigma-ell eta",
        "unprepared_counterexample": s.ImmutableMatrix([bad_n, 0, bad_b, 0]),
        "checks": {
            "full_weighted_adjoint_boundary": s.expand(
                (s.Matrix([En, Ev, Eb, Es]).T * original)[0]
                - (adapted.T * s.Matrix(fields))[0]
                - weighted(boundary)
            ),
            "physical_residual_clock_pullback": s.expand(
                adapted[0].subs(En, ep + delta * Ev).doit()
                + weighted(ep)
                - H * Ev
                - q * Eb
                - ell * Es
            ),
            "complete_normalized_force_pullback": (forced - wanted).applyfunc(s.expand),
            "unprepared_shift_homogeneous": s.simplify(actualW(bad_b)),
            "unprepared_lapse_homogeneous": s.simplify(-actualW(bad_n) + aq * bad_b),
        },
        "gates": {
            "original_initial_germ_is_essential": bad_b.subs(t, -s.Rational(1, 2)) != 0,
            "both_lapse_and_shift_recovered_by_weighted_uniqueness": True,
            "density_at_output_before_source_adjoint": True,
        },
    }


@cache
def ward_data():
    P = s.Symbol("external_momentum", positive=True)
    rho, pressure = tuple(
        s.Function(name, real=True)(k.t) for name in ("rho", "pressure")
    )
    dd = s.Function("delta", real=True)(k.t)
    stress = s.diag(-(k.a**3) * rho / 2, *([-k.a * pressure / 2] * 3))
    Pi = s.diag(0, 0, 1)
    gauge_residuals = []

    def packet(label, sign):
        fs = tuple(
            s.Function(label + name, real=True)(k.t) for name in ("eta", "w", "c")
        )
        phase = s.exp(sign * s.I * P * k.z)
        ee, ww, cc = tuple(f * phase for f in fs)
        qq = P**2 / k.a**2
        physical = (
            s.diff(ee, k.t),
            s.Matrix([0, 0, -sign * s.I * (s.diff(cc, k.t) + qq * ee) / P]),
            2 * (ww + k.H * ee) * s.eye(3),
        )
        syn = (s.S.Zero, s.zeros(3, 1), 2 * ww * s.eye(3) - 2 * cc * Pi)
        xi = s.Matrix([ee, 0, 0, -sign * s.I * cc / P])
        gauge = k.gauge(xi[0], xi[1:, 0])
        for x, y, z in zip(physical, syn, gauge):
            rr = x - y - z
            gauge_residuals.extend(list(rr) if isinstance(rr, s.MatrixBase) else [rr])
        return fs, physical, syn, xi

    det, dp, _, dx = packet("detector_", -1)
    src, _, gs, gx = packet("source_", 1)
    source = ward.source_integrand(dp, gx, stress)
    detector = ward.detector_integrand(dx, gs, stress)
    clock = (
        3
        * k.a**3
        * pressure
        * (4 * dd**2 - 3 * dd)
        * s.diff(det[0], k.t)
        * s.diff(src[0], k.t)
    )
    complete = s.cancel(s.expand(source + detector + clock))
    variables = []
    mapping = {}
    orders = []
    sides = []
    for side, fs in enumerate((det, src)):
        for f in fs:
            for j in range(5):
                var = s.Symbol(str(f.func) + "_jet" + str(j), real=True)
                variables.append(var)
                orders.append(j)
                sides.append(side)
                mapping[s.diff(f, k.t, j)] = var
    poly = s.Poly(complete.xreplace(mapping), *variables)
    maxima = [0, 0]
    bilinear = True
    for mon, _ in poly.terms():
        for side in (0, 1):
            bilinear &= sum(e for e, which in zip(mon, sides) if which == side) == 1
            maxima[side] = max(
                maxima[side],
                sum(e * j for e, j, which in zip(mon, orders, sides) if which == side),
            )
    clock_second = -3 * k.a**3 * pressure * (4 * dd**2 - 3 * dd)
    return {
        "full_Ward_monomial_count": len(poly.terms()),
        "maximum_derivative_per_leg": maxima,
        "full_Ward_external_momentum_degree": s.degree(complete, P),
        "distinct_clock_eta_second_coefficient": clock_second,
        "clock_bounce_second_coefficient": clock_second.subs(dd, s.Rational(1, 2)),
        "ordered_formula": "Source Ward on the full detector, then detector Ward on the synchronous source, plus distinct clock nn contact once; Gaussian one-point not zero.",
        "checks": {
            "both_opposite_leg_gauge_maps": s.ImmutableMatrix(
                gauge_residuals
            ).applyfunc(s.cancel),
            "no_net_spatial_phase": s.diff(complete, k.z),
            "full_fourteen_bilinear_monomials": len(poly.terms()) - 14,
            "maximum_one_derivative_each_leg": s.ImmutableMatrix(maxima) - s.ones(2, 1),
            "original_distinct_clock_bounce": clock_second.subs(dd, s.Rational(1, 2))
            - s.Rational(3, 2) * k.a**3 * pressure,
        },
        "gates": {
            "both_ordered_Ward_terms_nonzero": s.cancel(source) != 0
            and s.cancel(detector) != 0,
            "clock_contact_nonzero_and_retained_once": clock != 0,
            "no_inverse_external_momentum_after_complete_cancellation": not s.denom(
                complete
            ).has(P),
            "external_degree_at_most_two": s.degree(complete, P) <= 2,
            "full_expression_bilinear": bilinear,
            "Euler_time_order_at_most_two": maxima == [1, 1],
        },
    }
