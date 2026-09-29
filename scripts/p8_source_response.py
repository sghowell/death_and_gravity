"""Actual source-tagged constrained vertices and conditional Regge pole matching.

Read-only exact diagnostic. Inherits the frozen nonlinear constraint reduction
and transports its vertices to the ORIGINAL prepared canonical picture. The
finite-regulator Wick kernel is not a renormalized curved response bound.
Pole matching is conditional on the explicitly stated Regge/dispersion scope;
neither its hypotheses nor a finite-angle envelope are inferred from kappa.
"""

import hashlib
import json
from fractions import Fraction as F

import sympy as s
from p8_match1_rate_input import REPO, exact
from p8_rate4_candidate import original_borns
from p8_reference_regge import inputs as previous_inputs

EXTRA_SOURCES = {
    "scripts/p8_reference_regge.py": "1f0c84bddbc2bbe28f42cb01933126980e38d94a2ca627e63009e7f89ede3621",
    "docs/assessment-2026-09-24-p8-reference-regge.md": "bde4bc6454851fbedf59c882aaa3ff3dc88550de9bec99f19c66f555dd599ecc",
}
REPORTS = {
    "problems/P8/s6/continuation/s6_180/certificates/polynomial-vacuum-affine-coupled-response.json": "ad283ebc8d7159256c851b51a1f46056ec5ccb1924e34e6f6f02a28b0866aa3d",
    "problems/P8/s6/continuation/s6_253/certificates/polynomial-vacuum-affine-physical-background-vertices.json": "9b7e382f5ceb147457a8ee5d6eb5a5265834a2d7499ff5ec04ab5462195e927f",
    "problems/P8/s6/continuation/s6_257/certificates/polynomial-vacuum-affine-nonlinear-auxiliary-measure.json": "ae63ae3c213a5a5cd8b02174ec9eacfd4d19b3d0fb103c05aa221649ead4e05f",
    "problems/P8/s6/continuation/s6_275/certificates/polynomial-vacuum-affine-canonical-boundary-corrected-hybrid.json": "4951c18e25f7e51811592467e67e060231f9de7fd09f90722d08b96658fcc537",
}


def inputs():
    manifest = previous_inputs()
    for name, expected in {**EXTRA_SOURCES, **REPORTS}.items():
        assert hashlib.sha256((REPO / name).read_bytes()).hexdigest() == expected
        assert name not in manifest or manifest[name] == expected
        manifest[name] = expected
    for name in REPORTS:
        root = (REPO / name).parent.parent
        for relative, expected in json.loads((REPO / name).read_text())[
            "source_sha256"
        ].items():
            path = (root / relative).resolve()
            assert path.is_relative_to(root)
            assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
            key = str(path.relative_to(REPO))
            assert key not in manifest or manifest[key] == expected
            manifest[key] = expected
    return manifest


def source_rows(*, time, comoving_k_squared, lapse_pivot, profile_current):
    """Rows on y=(v,sigma,pv,ps); momenta divided by kappa*a^3.

    J and T are explicit inputs, not assigned their tree values by default.
    The new interacting-vertex proof uses the formal grade-zero clock root;
    inserting the fixed grade-one profiles here is a higher-order resummation,
    not a proof of a nonlinear root with those nonzero profiles.
    n is the linear lapse; P=pv+3*ell*sigma. Spatial derivatives on the
    external and internal legs must use their OWN momenta in a convolution.
    """
    u, k2, J, T = map(exact, (time, comoving_k_squared, lapse_pivot, profile_current))
    if k2 < 0 or J <= 0:
        raise ValueError("Require k_squared>=0 and a positive regular lapse pivot")
    a, h = (1 + u * u) ** 2, (1 + u * u) ** 3
    ell, E = 1 / (10 * a**3), 1 - 3 / (2 * h)
    theta = 4 * u / (1 + u * u) - u / (1 + u * u) ** 4
    q = k2 / a**2
    n = s.Matrix([[-(2 * E * q + 3 * T), 3 * theta * ell, theta, ell * E]]) / (2 * J)
    P = s.Matrix([[0, 3 * ell, 1, 0]])
    return n, P, s.Rational(theta), s.Rational(h)


def source_matrices(n, P, theta):
    """H3=-G*y^T F y/(2h); H4_tag2=G^2*y^T K y/h^2, before 1/J.

    Return F and the row A=P/2-6*theta*n. The contact matrix is
    A.T*A/J-3*n.T*n. Products are Weyl ordered at the same regulator.
    """
    return n.T * P + P.T * n - 12 * theta * n.T * n, P / 2 - 6 * theta * n


def polarized_source_matrices(n_left, P_left, n_right, P_right, theta, J):
    """Different leg momenta for a convolution, NOT a diagonal-mode shortcut.

    F_lr=F_rl.T and K_lr=K_rl.T; the matrices need not be symmetric for
    unequal leg momenta. H3=-G*y_left.T*F_lr*y_right/(2h), summed over
    ordered momenta, and H4_tag2=G^2*y_left.T*K_lr*y_right/h^2.
    """
    F_lr = n_left.T * P_right + P_left.T * n_right - 12 * theta * n_left.T * n_right
    A_left, A_right = P_left / 2 - 6 * theta * n_left, P_right / 2 - 6 * theta * n_right
    K_lr = A_left.T * A_right / J - 3 * n_left.T * n_right
    return F_lr, K_lr


def pole_matched_ratio_upper(*, kappa, amplitude_factor, slope_fraction):
    """C_envelope/alpha_min < (11/7)*B/(kappa*a_rel), CONDITIONALLY.

    Requires the written single leading Regge-pole matching hypotheses,
    |C(t)|<=B*C(0) and (2-j(t))/(-t)>=a_rel*j'(0) on the ENTIRE detector
    interval. No default B or a_rel is a physical bound. pi/2<11/7.
    """
    kappa, B, a_rel = map(exact, (kappa, amplitude_factor, slope_fraction))
    if kappa <= 0 or B <= 0 or a_rel <= 0:
        raise ValueError("Require positive coupling normalization and shape bounds")
    return F(11, 7) * B / (kappa * a_rel)


def omega(n):
    return s.zeros(n).row_join(s.eye(n)).col_join((-s.eye(n)).row_join(s.zeros(n)))


def mixed_hessian(row, g):
    return (
        s.zeros(row.cols)
        .row_join(-row.T * g)
        .col_join((-g.T * row).row_join(s.zeros(g.cols)))
    )


def audit():
    before = inputs()
    # Independent joint trace/temporal Legendre solve, not a Proca-only solve.
    a, B, U, p, G, r, c, xi, K, T = s.symbols("a B U p G r c xi K T", nonzero=True)
    lag = a * K**2 + B * K + U * (T - xi * (r * K + c)) ** 2 / 2
    solved = s.solve((p - s.diff(lag, K), -s.diff(lag, T) - G), (K, T))
    reduced = (p - B - xi * r * G) ** 2 / (4 * a) + G**2 / (2 * U) - xi * c * G
    assert s.cancel((p * K - lag - T * G).subs(solved) - reduced) == 0
    assert s.cancel(solved[K] - (p - B - xi * r * G) / (2 * a)) == 0

    # Literal clock expansion. Only needed jets are used, with their frozen
    # full-source factorization inherited; no global polynomial replacement.
    e, n, H, M1, b0, b1, v, raw, ell, sig, h = s.symbols(
        "e n H M1 b0 b1 v raw ell sig h", nonzero=True
    )
    N, M, b = 1 + e * n, 1 + M1 * e * n, b0 + b1 * e * n
    pbar = b0 - 2 * H
    trace_p = (pbar + e * raw / 3) * s.exp(-3 * e * v)
    clock_r = ((1 + e * n) ** -2 - 1) / h
    # xi-linear Hamiltonian coefficient after T elimination, divided by G.
    literal = clock_r * (3 * N * (trace_p - b) / (2 * M) + 3 * H)
    jet = s.series(literal, e, 0, 3).removeO().expand()
    assert jet.coeff(e, 0) == jet.coeff(e, 1) == 0
    theta = -H * (M1 - 1) + b1 / 2
    prepared = s.expand(jet.coeff(e, 2).subs(raw, p + 9 * pbar * v + 3 * ell * sig))
    P = p + 3 * ell * sig
    cubic = -n * (P - 6 * theta * n) / h
    assert s.factor(prepared - cubic) == 0
    direct = s.series(-3 * N * clock_r**2 / (4 * M), e, 0, 3).removeO()
    assert s.expand(direct).coeff(e, 2) == -3 * n**2 / h**2
    # Dropping the original cross-boundary shift changes the actual vertex.
    omitted = s.expand(jet.coeff(e, 2).subs(raw, p + 9 * pbar * v))
    assert s.factor(prepared - omitted) == -3 * ell * n * sig / h

    # Actual R/M/B clock jets give the inherited Theta, at arbitrary time.
    u, lapse = s.symbols("u lapse", real=True, positive=True)
    clock_h = (1 + u**2) ** 3
    R = 1 + (lapse**-2 - 1) / clock_h
    M_clock = R ** s.Rational(1, 4)
    B_clock = -(R ** -s.Rational(3, 4)) * s.diff(R, u) / (2 * lapse)
    # I_N(u,1)=0 is inherited and pinned; I_NN is NOT zero and is retained in J.
    Hu = 4 * u / (1 + u**2)
    theta_clock = (
        -Hu * (s.diff(M_clock, lapse).subs(lapse, 1) - 1)
        + s.diff(B_clock, lapse).subs(lapse, 1) / 2
    )
    assert s.factor(theta_clock - Hu + u / (1 + u**2) ** 4) == 0

    # Independent exact lapse elimination in the degree-relevant polynomial.
    J, L, th, PP = s.symbols("J L th PP", nonzero=True)
    Haux = -J * n**2 + L * n - e * xi * G * n * (PP - 6 * th * n) / h
    Haux -= 3 * e**2 * xi**2 * G**2 * n**2 / h**2
    nstar = (L - e * xi * G * PP / h) / (
        2 * (J - 6 * e * xi * th * G / h + 3 * e**2 * xi**2 * G**2 / h**2)
    )
    assert s.cancel(s.diff(Haux, n).subs(n, nstar)) == 0
    effective = s.cancel(Haux.subs(n, nstar))
    jet = s.series(effective, e, 0, 3).removeO().expand()
    n1 = L / (2 * J)
    expected4 = xi**2 * G**2 / h**2 * ((PP / 2 - 6 * th * n1) ** 2 / J - 3 * n1**2)
    assert s.factor(jet.coeff(e, 1) + xi * G * n1 * (PP - 6 * th * n1) / h) == 0
    assert s.factor(jet.coeff(e, 2) - expected4) == 0
    assert s.factor(expected4 + 3 * xi**2 * G**2 * n1**2 / h**2) != 0

    # Bind the bounce J to the pinned FULL primitive-sensitive tree expression.
    report = json.loads((REPO / next(iter(REPORTS))).read_text())

    def find_J(value):
        if isinstance(value, dict):
            if "J" in value and isinstance(value["J"], str):
                yield value["J"]
            for item in value.values():
                yield from find_J(item)
        elif isinstance(value, list):
            for item in value:
                yield from find_J(item)

    candidates = list(find_J(report))
    assert len(candidates) == 1
    tree_J = s.sympify(candidates[0], locals={"u": u})
    assert tree_J.subs(u, 0) == s.Rational(243, 160)
    nrow, Prow, trow, hrow = source_rows(
        time=0, comoving_k_squared=1, lapse_pivot=F(243, 160), profile_current=0
    )
    assert nrow == s.Matrix([[s.Rational(80, 243), 0, 0, -s.Rational(4, 243)]])
    assert Prow == s.Matrix([[0, s.Rational(3, 10), 1, 0]])
    assert trow == 0 and hrow == 1
    v0, sigma0, pv0, ps0 = s.symbols("v0 sigma0 pv0 ps0", real=True)
    y = s.Matrix([v0, sigma0, pv0, ps0])
    n_bounce = (nrow * y)[0]
    assert n_bounce == (80 * v0 - 4 * ps0) / 243
    # Both raw momenta must be transported; not a fresh Gaussian preparation.
    raw_ps = ps0 + s.Rational(3, 10) * v0
    n_raw = ((1 + s.Rational(3, 200)) * v0 - raw_ps / 20) / (2 * s.Rational(243, 160))
    assert s.expand(n_raw - n_bounce) == 0
    f3, ar = source_matrices(nrow, Prow, trow)
    contact = ar.T * ar / s.Rational(243, 160) - 3 * nrow.T * nrow
    assert (y.T * contact * y)[0].subs({v0: 0, sigma0: 0, pv0: 1, ps0: 0}) > 0
    assert (y.T * contact * y)[0].subs({v0: 1, sigma0: 0, pv0: 0, ps0: 0}) < 0
    assert s.expand((y.T * f3 * y)[0] / 2 - n_bounce * (Prow * y)[0]) == 0

    # Exact Wick controls with correlated physical Gaussian fixtures. These
    # are algebra tests, NOT replacements for the inherited prepared state.
    o4, o2 = omega(2), omega(1)
    shear = (
        s.eye(2)
        .row_join(s.zeros(2))
        .col_join(
            s.Matrix([[1, s.Rational(1, 3)], [s.Rational(1, 3), -2]]).row_join(s.eye(2))
        )
    )
    assert shear * o4 * shear.T == o4
    V = (
        shear
        * s.diag(s.Rational(1, 2), s.Rational(1, 6), s.Rational(1, 2), s.Rational(3, 2))
        * shear.T
    )
    assert V * o4 * V == o4 / 4
    S = s.Matrix(
        [
            [s.Rational(3, 5), 0, s.Rational(4, 5), 0],
            [0, s.Rational(5, 13), 0, s.Rational(12, 13)],
            [-s.Rational(4, 5), 0, s.Rational(3, 5), 0],
            [0, -s.Rational(12, 13), 0, s.Rational(5, 13)],
        ]
    )
    Rv = s.Matrix(
        [
            [s.Rational(8, 17), s.Rational(15, 17)],
            [-s.Rational(15, 17), s.Rational(8, 17)],
        ]
    )
    assert S * o4 * S.T == o4 and Rv * o2 * Rv.T == o2
    W = S * (V + s.I * o4 / 2)
    WV = Rv * (s.diag(s.Rational(1, 4), 1) + s.I * o2 / 2)
    fullW = s.diag(W, WV)
    g = s.Matrix([[0, 2]])
    wick_checks = 0
    for probe_left, probe_right in (
        (s.eye(4)[:, 0], s.eye(4)[:, 2]),
        (s.Matrix([1, 2, 3, 4]), s.Matrix([2, -1, 4, 3])),
        (s.eye(4)[:, 1], s.eye(4)[:, 3]),
    ):
        lx, ly = probe_left.T * f3, probe_right.T * f3
        A, B = mixed_hessian(lx, g), mixed_hessian(ly, g)
        connected = s.expand(s.trace(A * fullW * B * fullW.T) / 2)
        product = s.expand((g * WV * g.T)[0] * (lx * W * ly.T)[0])
        assert s.expand(connected - product) == 0
        assert (
            s.expand(-s.I * (product - s.conjugate(product)) - 2 * s.im(product)) == 0
        )
        assert s.im(product) != 0
        # Ordinary transpose is essential for this Wightman ordering.
        wrong = s.expand(s.trace(A * fullW * B * fullW.conjugate().T) / 2)
        assert wrong != connected
        wick_checks += 1
    # The longitudinal external-G bubble has two internal scalar lines. Its
    # factor 1/2 is independently checked by the indexed Wick contractions.
    scalar_bubble = s.trace(f3 * W * f3 * W.T) / 2
    indexed = sum(
        f3[i, j] * f3[k, l] * (W[i, k] * W[j, l] + W[i, l] * W[j, k]) / 4
        for i in range(4)
        for j in range(4)
        for k in range(4)
        for l in range(4)
    )
    assert s.expand(scalar_bubble - indexed) == 0
    assert s.im(s.expand(scalar_bubble)) != 0
    # The quartic's two local diagonal blocks are retained, with mixed mean
    # zero in the unchanged reference product state (not a factorized scalar state).
    assert s.hessian((y.T * contact * y)[0], list(y)) == 2 * contact
    g_external = s.Symbol("g_external", real=True)
    assert s.diff(g_external**2 * s.trace(contact * V), g_external, 2) == 2 * s.trace(
        contact * V
    )
    row_checks = 0
    for time in (F(0), F(1, 8), F(-1, 8)):
        for k2 in (F(0), F(1), F(9)):
            nr, pr, tr, hr = source_rows(
                time=time,
                comoving_k_squared=k2,
                lapse_pivot=F(3, 2),
                profile_current=F(1, 100),
            )
            fm, aa = source_matrices(nr, pr, tr)
            assert fm == fm.T and hr > 0
            assert (
                s.expand(
                    (y.T * fm * y)[0] / 2
                    - (nr * y)[0] * ((pr * y)[0] - 6 * tr * (nr * y)[0])
                )
                == 0
            )
            km = aa.T * aa / s.Rational(3, 2) - 3 * nr.T * nr
            assert km == km.T
            row_checks += 1
    polarized_checks = 0
    for time in (F(0), F(1, 8), F(-1, 8)):
        for left_k, right_k in ((0, 1), (1, 9)):
            rows = [
                source_rows(
                    time=time,
                    comoving_k_squared=kk,
                    lapse_pivot=F(3, 2),
                    profile_current=0,
                )
                for kk in (left_k, right_k)
            ]
            nl, pl, theta_lr, _ = rows[0]
            nr, pr, _, _ = rows[1]
            flr, klr = polarized_source_matrices(
                nl, pl, nr, pr, theta_lr, s.Rational(3, 2)
            )
            frl, krl = polarized_source_matrices(
                nr, pr, nl, pl, theta_lr, s.Rational(3, 2)
            )
            assert flr == frl.T and klr == krl.T
            assert flr != source_matrices(nl, pl, theta_lr)[0]
            probe = s.Matrix([1, 2, 3, 4])
            leg_row = (nl * probe)[0] * pr + (
                (pl * probe)[0] - 12 * theta_lr * (nl * probe)[0]
            ) * nr
            assert probe.T * flr == leg_row
            polarized_checks += 1

    # Frozen flat vacuum longitudinal contact control, not the curved state's
    # finite part. G_k covariance=k^2/(2*kappa*sqrt(k^2+m^2)).
    k, mass = s.symbols("k mass", positive=True)
    primitive = (
        k * s.sqrt(k**2 + mass**2) * (2 * k**2 - 3 * mass**2)
        + 3 * mass**4 * s.asinh(k / mass)
    ) / 8
    assert s.simplify(s.diff(primitive, k) - k**4 / s.sqrt(k**2 + mass**2)) == 0
    assert primitive.subs(k, 0) == 0
    assert s.limit(primitive / k**4, k, s.oo) == s.Rational(1, 4)

    # Our massive s0=2 two-cut convention, independently of massless FESRs.
    t, C0, C1, slope, curvature, kap = s.symbols(
        "t C0 C1 slope curvature kap", nonzero=True
    )
    singular_tail = 2 * (C0 + C1 * t) / (s.pi * (-slope * t - curvature * t**2 / 2))
    pole_coefficient = s.limit(t * singular_tail, t, 0)
    assert pole_coefficient == -2 * C0 / (s.pi * slope)
    assert pole_coefficient.subs(C0, s.pi * slope / (2 * kap)) == -1 / kap
    finite = s.limit(singular_tail - pole_coefficient / t, t, 0)
    assert (
        s.factor(finite + 2 * C1 / (s.pi * slope) - C0 * curvature / (s.pi * slope**2))
        == 0
    )
    split, energy = s.symbols("split energy", positive=True)
    denominator_remainder = (
        s.log(split / (split - 2)) + 4 / (split - 2) + 2 / (split - 2) ** 2
    )
    integrand = energy**2 / (energy - 2) ** 3 - 1 / energy
    assert (
        s.cancel(s.diff(denominator_remainder, split) + integrand.subs(energy, split))
        == 0
    )
    assert s.limit(denominator_remainder, split, s.oo) == 0
    # At n=0,1 a leading Regge FESR determines a trajectory only when its
    # finite-energy moments and complex-arc remainder have actually been given.
    j, f = s.symbols("j f", nonzero=True)
    m1, m3 = f / (j + 2), f / (j + 4)
    assert s.cancel((4 * m3 - 2 * m1) / (m1 - m3) - j) == 0
    _, _, kappa, lam, _ = original_borns()
    ratio = pole_matched_ratio_upper(
        kappa=kappa, amplitude_factor=10**196, slope_fraction=1
    )
    assert 60 * ratio < lam / 100
    assert 60 * F(11, 7) < 95
    # Increasing the finite split does not automatically reduce the residue:
    # C_T2(t)=C_T1(t)*(T2/T1)^(j(t)-2); at t0 the slope gains j'(0)*log(T2/T1).
    scale = s.Symbol("scale", positive=True)
    transported = (C0 + C1 * t) * scale ** (slope * t + curvature * t**2 / 2)
    assert (
        s.simplify(s.diff(transported, t).subs(t, 0) - C1 - C0 * slope * s.log(scale))
        == 0
    )

    rejected = 0
    good_rows = {
        "time": F(0),
        "comoving_k_squared": F(1),
        "lapse_pivot": F(3, 2),
        "profile_current": F(0),
    }
    good_pole = {"kappa": F(2), "amplitude_factor": F(3), "slope_fraction": F(1, 2)}
    for function, kwargs in (
        (source_rows, good_rows),
        (pole_matched_ratio_upper, good_pole),
    ):
        for name in kwargs:
            try:
                function(**{key: value for key, value in kwargs.items() if key != name})
            except TypeError:
                rejected += 1
            else:
                raise AssertionError("A missing physical input was defaulted")
        for name in kwargs:
            for bad in (True, 0.5):
                try:
                    function(**{**kwargs, name: bad})
                except (TypeError, ValueError):
                    rejected += 1
                else:
                    raise AssertionError("Inexact or boolean input accepted")
    for function, kwargs, update in (
        (source_rows, good_rows, {"comoving_k_squared": -1}),
        (source_rows, good_rows, {"lapse_pivot": 0}),
        (source_rows, good_rows, {"lapse_pivot": -1}),
        (pole_matched_ratio_upper, good_pole, {"kappa": 0}),
        (pole_matched_ratio_upper, good_pole, {"amplitude_factor": 0}),
        (pole_matched_ratio_upper, good_pole, {"slope_fraction": 0}),
    ):
        try:
            function(**{**kwargs, **update})
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("Invalid domain accepted")
    assert inputs() == before
    return {
        "milestone": "P8.PRESCRIPTION-2.CONSTRAINED-SOURCE + GRAVITY.POLE-MATCH",
        "outcome": "ACTUAL_SOURCE_VERTICES_AND_FINITE_REGULATOR_RESPONSE_KERNEL_WITH_CONDITIONAL_POLE_NORMALIZATION",
        "protected_input_files": len(before),
        "physical_model_adopted": False,
        "tracks_complete": False,
        "inherited_not_new": "S250 heavy filtration; S253 joint Gaussian reduction/state; S257 finite-regulator constraint measure; S275 both original canonical boundary shifts",
        "source_cubic": "h3_tag1=-xi*G*n*(P-6*Theta*n)/h; P=pv+3*ell*sigma",
        "source_quartic": "h4_tag2=xi^2*G^2*((P/2-6*Theta*n)^2/J-3*n^2)/h^2",
        "source_scope": "source-squared one-loop quadratic blocks: mixed scalar-vector fish, two-scalar longitudinal bubble and both lapse-generated contacts; formal tree-clock vertices and fixed-state finite-regulator kernels, not a renormalized continuum bound",
        "bounce_tree_J": "243/160",
        "bounce_contact_sign": "indefinite; no sign or smallness claimed for the full response",
        "arbitrary_time_momentum_row_controls": row_checks,
        "unequal_momentum_polarization_controls": polarized_checks,
        "correlated_wick_controls": wick_checks,
        "longitudinal_bubble_and_both_contact_blocks": True,
        "contact_uv_control": "frozen flat vacuum <G^2> has cutoff^4/(16*pi^2*kappa) leading term; not a physical cutoff or the curved renormalized answer",
        "conditional_regge_matching": "C(0)/j'(0)=pi/(2*kappa) only for the specified single leading-pole dispersion scope",
        "conditional_finite_angle_tail": "<95*B/(kappa*a_rel) plus supplied residual error; B and a_rel still need physical interval bounds",
        "fesr_scope": "published massless FESRs require finite-energy moments and a Regge complex-arc approximation; neither follows from four RATE4 samples",
        "invalid_or_missing_inputs_rejected": rejected,
        "remaining": "renormalize and bound the same-state source kernel with other reference sectors; establish actual Regge/finite-window/contour/observable inputs; M/V/G/B/R/P8 OPEN",
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
