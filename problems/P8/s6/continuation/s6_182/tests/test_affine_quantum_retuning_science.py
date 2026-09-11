"""Independent full-coefficient, covariant-current and nonadaptive controls."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_proca_gaussian import stress
from p8_vacuum_affine_quantum_retuning import audit, bounds, profile, response, ward


def mpq(value):
    value = s.Rational(value)
    return mp.mpf(int(value.p)) / int(value.q)


def T(x):
    return x**1024 / (x**1024 + (1 - x) ** 1024)


def T_minus_one(x):
    q = (1 - x) / x
    return -(q**1024) / (1 + q**1024)


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_all_exact_residuals(name):
    v = audit.residuals()[name]
    assert all(x == 0 for x in (list(v) if isinstance(v, s.MatrixBase) else [v])), name


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_unsupported_scope_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("name", list(audit.packets()))
def test_every_quantitative_packet(name):
    assert all(audit.packets()[name]["gates"].values())


@pytest.mark.parametrize("order", range(6))
def test_reference_bound_is_actual_source_pinned_stress(order):
    row = stress.data()["bounds_by_time_derivative_order"][order]
    for kind in ("energy", "pressure"):
        assert (
            0
            < row["absolute_" + kind + "_derivative_upper"] / profile.KAPPA
            < profile.EPS
        )
    assert "complete" in profile.data()["fixed_profile_definition"]
    assert (
        "radial integral"
        in profile.data()["source_defined_energy_and_pressure"]["definition"]
    )


@pytest.mark.parametrize(
    "x",
    (
        -s.Rational(1, 8192),
        0,
        s.Rational(1, 4),
        s.Rational(1, 2),
        s.Rational(3, 4),
        1,
        s.Rational(119, 100),
    ),
)
def test_literal_complex_discs_cover_full_original_strip(x):
    with mp.workdps(90):
        xx = mpq(x)
        radius = mpq(bounds.data()["full_X_strip_complex_disc_radius"])
        for k in range(16):
            z = radius * mp.exp(2j * mp.pi * k / 16)
            numerator = (xx + z) ** 1024
            other = (1 - xx - z) ** 1024
            assert numerator + other != 0
            assert abs(numerator / (numerator + other)) < 2


@pytest.mark.parametrize("x", (s.Rational(1, 4), s.Rational(1, 2), s.Rational(3, 4)))
def test_literal_middle_disc_phase_bound(x):
    with mp.workdps(90):
        xx = mpq(x)
        d = mpq(bounds.data()["full_X_strip_complex_disc_radius"])
        for k in range(16):
            z = d * mp.exp(2j * mp.pi * k / 16)
            q = (1 - xx - z) / (xx + z)
            angle = 1024 * mp.arg(q)
            assert abs(angle) <= mp.mpf(1) / 4
            assert mp.re(q**1024) > 0
            assert abs(1 / (1 + q**1024)) <= 1


@pytest.mark.parametrize("x", (s.Rational(3, 4), s.Rational(4, 5), s.Rational(6, 5)))
def test_nonzero_complete_clock_remainder_derivatives(x):
    with mp.workdps(110):
        xx = mpq(x)
        for order, upper in bounds.data()[
            "Cauchy_T_minus_one_derivative_envelopes"
        ].items():
            actual = abs(mp.diff(T_minus_one, xx, order))
            assert 0 < actual < mpq(upper)
        # An independent literal denominator evaluation at enough precision
        # confirms that complement evaluation is not a new switch.
        with mp.workdps(1100):
            assert abs((T(xx) - 1) - T_minus_one(xx)) < abs(T_minus_one(xx)) * mp.mpf(
                "1e-80"
            )


@pytest.mark.parametrize(
    "time_order,X_order", [(j, k) for j in range(6) for k in range(6)]
)
def test_all_mixed_error_budget_entries(time_order, X_order):
    d = bounds.data()
    assert d["full_X_strip_mixed_correction_derivative_upper"][
        (time_order, X_order)
    ] < s.Rational(1, 10**742)
    assert (
        d["mixed_correction_derivative_upper"][(time_order, X_order)] <= 2 * profile.EPS
    )
    if X_order >= 2:
        assert d["mixed_correction_derivative_upper"][
            (time_order, X_order)
        ] <= s.Rational(1, 10**1159)


@pytest.mark.parametrize("order", range(6))
def test_literal_transition_derivatives_with_adversarial_bounded_profiles(order):
    # These are coefficient-bound fixtures, NOT replacement quantum states.
    with mp.workdps(90):
        for x in (mp.mpf("0.25"), mp.mpf("0.5"), mp.mpf("0.75"), mp.mpf("1.19")):
            for t in (mp.mpf("-0.5"), mp.mpf(0), mp.mpf("0.5")):
                r, p = mp.sin(t), mp.cos(t)
                pv = mp.mpf("1e-19")
                f = lambda y, pv=pv, p=p, r=r: (
                    -pv + T(y) * (-p + pv - (r + p) * (y - 1) / 2)
                )
                actual = abs(mp.diff(f, x, order))
                envelope = (
                    bounds.data()["full_X_strip_mixed_correction_derivative_upper"][
                        (0, order)
                    ]
                    / profile.EPS
                )
                assert actual < mpq(envelope)


@pytest.mark.parametrize("r,p", ((2, 3), (-1, 1), (5, -2), (0, 0)))
def test_literal_lapse_volume_current_at_fixed_reference_coefficients(r, p):
    N, v = s.symbols("N v", positive=True)
    F = -p - s.Rational(r + p, 2) * (N**-2 - 1)
    L = N * s.exp(3 * v) * F
    point = {N: 1, v: 0}
    mean = s.Matrix([s.diff(L, N).subs(point), s.diff(L, v).subs(point)])
    Hess = s.hessian(L, (N, v)).subs(point)
    assert mean == s.Matrix([r, -3 * p])
    assert Hess == s.Matrix([[-r - p, 3 * r], [3 * r, -9 * p]])
    # This is not the derivative of a freshly reselected stress.
    assert (
        s.factor(s.diff((N**-2 + 1) * (-s.Rational(r + p, 2)) + p, N).subs(N, 1))
        == r + p
    )
    assert s.diff(F, N).subs(N, 1) == r + p


@pytest.mark.parametrize(
    "time",
    (-s.Rational(1, 2), -s.Rational(1, 4), 0, s.Rational(1, 4), s.Rational(1, 2)),
)
def test_complete_chart_lapse_pivot_by_independent_high_precision(time):
    with mp.workdps(90):
        tt = mpq(time)
        delta = 1 / (2 * (1 + tt * tt) ** 3)
        r, p = mp.mpf(2) / 7, mp.mpf(3) / 11
        pv = mp.mpf("1e-19")
        correction = lambda X: -pv + T(X) * (-p + pv - (r + p) * (X - 1) / 2)
        density = lambda N: (
            N * (1 + 2 * delta * (N**-2 - 1)) ** (-mp.mpf(3) / 4) * correction(N**-2)
        )
        literal = mp.diff(density, mp.mpf(1), 2) / 2
        expected = (21 * delta * delta - 3 * delta) * (-p) / 2 + (1 - 6 * delta) * (
            -(r + p) / 2
        )
        assert abs(literal - expected) < mp.mpf("1e-70")
        assert abs(literal) < 4
    assert bounds.data()["new_classical_J_lower_on_slab"] > s.Rational(1, 100)


@pytest.mark.parametrize("order", range(6))
def test_full_switch_endpoint_jets_independent_literal_derivative(order):
    X = profile.X
    actual = profile.original.data()["T"]
    assert s.diff(actual, X, order).subs(X, 0) == 0
    assert s.diff(actual - 1, X, order).subs(X, 1) == 0


def test_unconserved_arbitrary_profile_does_not_solve_scalar_mean_equation():
    u = profile.u
    expression = ward.data()["literal_clock_profile_scalar_Euler"]
    changed = expression.subs(
        {profile.rho: 0, profile.P: u, ward.H: 4 * u / (1 + u * u)}, simultaneous=True
    ).doit()
    assert s.factor(changed - 12 * u * u / (1 + u * u)) == 0
    assert changed.subs(u, s.Rational(1, 4)) != 0


@pytest.mark.parametrize("c", (-3, 1, 5))
def test_same_constant_vacuum_density_cancels_contacts_not_only_mean(c):
    new = response.data()["new_complete_stationary_adapted_response"]
    mapped = new.subs(
        {profile.rho: -c, profile.P: c, response.R: 0, response.Q: 0}, simultaneous=True
    ).doit()
    assert mapped.applyfunc(s.simplify) == s.zeros(2, 1)
    individual = ward.data()["profile_current_Hessian"].subs(
        {profile.rho: -c, profile.P: c}
    )
    assert individual != s.zeros(2, 2)
    assert individual == s.Matrix([[0, -3 * c], [-3 * c, -9 * c]])


def test_actual_nonlocal_response_survives_and_profile_is_not_adaptive():
    d = response.data()
    current = d["new_complete_stationary_adapted_response"]
    assert s.diff(current[1], response.Q) == 3
    assert (
        s.diff(current[0], s.diff(response.eta, profile.u, 2))
        == profile.rho + profile.P
    )
    fixed = profile.data()["full_coefficient_correction"]
    assert s.diff(fixed, response.Q) == 0 and s.diff(fixed, response.R) == 0


def test_exact_new_local_lapse_recovery_and_no_small_full_inverse_claim():
    d = bounds.data()
    assert d["new_first_row_quantum_added_upper"] == s.Rational(59, 2) * profile.EPS
    assert (
        d["new_retained_first_row_C1_upper"]
        == bounds.old_inverse.first_row()["physical_lapse_reconstruction_C1_upper"]
        + 19 * profile.EPS
    )
    assert s.Rational(14976) < d["new_retained_first_row_C1_upper"] < 15000
    assert "remain unevaluated" in response.data()["uncomputed"]


@pytest.mark.parametrize("ratio", (1, 2, 7))
def test_complete_canonical_retuning_is_fixed_not_just_a_quartic_germ(ratio):
    d = profile.data()["full_coefficient_correction"]
    fixed_physical_at_base_coordinates = d / s.Integer(ratio)
    assert (
        s.simplify(
            (profile.KAPPA * ratio) * fixed_physical_at_base_coordinates
            - profile.KAPPA * d
        )
        == 0
    )
    assert profile.data()["vacuum_nonconstant_first_possible_field_degree"] == 2048


def test_nonconstant_vacuum_factor_does_not_require_time_reflection_symmetry():
    X = profile.X
    denominator = X**1024 + (1 - X) ** 1024
    correction = profile.data()["full_coefficient_correction"]
    bracket = -profile.P + profile.PV - (profile.rho + profile.P) * (X - 1) / 2
    assert s.factor((correction + profile.PV) * denominator - X**1024 * bracket) == 0
    # An odd coefficient can occur, but its first field degree is2049.
    e, Y, phi = s.symbols("epsilon Y phi", real=True)
    odd = e * phi * (e * e * Y) ** 1024
    assert s.degree(odd, e) == 2049


def test_new_action_and_mean_boundary_are_explicit():
    assert profile.data()["candidate"] == "CD-REG-AFFINE-ISO-QG1"
    assert "No old certificate" in profile.data()["parent"]
    assert ward.data()["gates"]["only_retained_conditional_mean_stationarity"] is True
    assert "not a rewrite" in audit.observable()["new_candidate"]
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 38
    assert len(audit.residuals()) == 63 and audit.scalar_entry_count() == 74
    assert len(audit.gates()) == 36 and audit.rejected_inputs() == 104
