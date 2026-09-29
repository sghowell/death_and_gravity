"""Source-kernel extension and physical-completion obstruction audit.

The continuum proof is written in the companion note. Exact controls bind its
power counting and two covariant finite-functional witnesses to the actual
prepared response. Witnesses are NOT adopted changes of the parent. They show
that the retained normalizations cannot choose a common quantum completion.
No O(G) estimate is promoted to a certified physical gravitational allowance.
"""

import hashlib
import json
from fractions import Fraction as F
from itertools import product
from math import comb

import sympy as s
from p8_match1_rate_input import REPO, exact
from p8_rate4_candidate import original_borns
from p8_source_response import inputs as previous_inputs
from p8_source_response import source_rows

EXTRA_SOURCES = {
    "scripts/p8_source_response.py": "3286d7f91f0d4abbcb3c2b78128585b51cddac3f4ba9fbf19bc73a64a44b19b3",
    "docs/assessment-2026-09-24-p8-source-response-pole-match.md": "91f61824c06cbc99eb4b242cc568cc82ca2f6aa06b73e481fa5eddaf25698dce",
}
REPORTS = {
    "problems/P8/s6/continuation/s6_176/certificates/polynomial-vacuum-affine-proca-gaussian.json": "1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75",
    "problems/P8/s6/continuation/s6_241/certificates/polynomial-vacuum-affine-heavy-clock-quadratic.json": "f7d3b2ddfaf2d8fb2e09dcaf28ca942ddcda8c0ac36ad9e95711b2db9a1e5a04",
    "problems/P8/s6/continuation/s6_251/certificates/polynomial-vacuum-affine-coupled-gaussian-state.json": "e60b97e9a4634795ea788ac6a3d35618c938bad000cb57cc5a85a40bd43b1f5c",
    "problems/P8/s6/continuation/s6_252/certificates/polynomial-vacuum-affine-gaussian-measure-response.json": "504045fc426a73d2de1d22d5532cdf0238b6659cbbe7c24e162abf8ba2a5e3da",
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


def normal_indices(order):
    """Four NORMAL coordinates to the diagonal; not eight independent ones."""
    if type(order) is not int or not 0 <= order <= 10:
        raise ValueError("Use the proved finite normal-jet orders 0 through 10")
    return tuple(
        index for index in product(range(order + 1), repeat=4) if sum(index) <= order
    )


def polynomial_jet_remainder(terms, order):
    """Algebra control for the local Taylor projector, NOT a distribution solver.

    The actual test-function projector has a smooth compact cutoff equal to
    one near the diagonal. This exact polynomial control checks its local jet.
    """
    allowed = set(normal_indices(order))
    checked = {}
    for index, value in terms.items():
        if (
            not isinstance(index, tuple)
            or len(index) != 4
            or any(type(i) is not int or i < 0 for i in index)
        ):
            raise ValueError("Require four nonnegative integer normal exponents")
        checked[index] = exact(value)
    return {index: value for index, value in checked.items() if index not in allowed}


def audit():
    before = inputs()
    # Source-specific symbol counting uses the inherited per-field orders,
    # rather than adding two derivatives to a worst-entry covariance bound.
    position, momentum = F(-1, 2), F(3, 2)
    lapse_order = max(position + 2, momentum)
    P_order = max(position, momentum)
    longitudinal_momentum = F(-1, 2)
    G_order = longitudinal_momentum + 1
    scalar_internal = 6 + 2 * max(lapse_order, P_order) + 2 * G_order
    long_internal = 6 + 4 * max(lapse_order, P_order)
    assert (scalar_internal, long_internal) == (10, 12)
    phase_scalar, phase_long = scalar_internal + 4, long_internal + 2
    assert phase_scalar == phase_long == 14
    derivative_cap = int(phase_scalar - 4)
    assert derivative_cap == 10
    assert 4 + (derivative_cap + 1) - phase_scalar == 1
    jet_counts = {}
    for order in (6, 8, 10):
        indices = normal_indices(order)
        assert len(indices) == comb(order + 4, 4)
        poly = {index: F((-1) ** sum(index), 1 + sum(index)) for index in indices}
        poly[(order + 1, 0, 0, 0)] = F(7, 13)
        assert polynomial_jet_remainder(poly, order) == {(order + 1, 0, 0, 0): F(7, 13)}
        jet_counts[str(order)] = len(indices)

    # Actual longitudinal normalization on the reference FLRW clock. This is
    # the leading symbol of the SAME all-order Proca state, not a new vacuum.
    k, a, zeta, kappa = s.symbols("k a zeta kappa", positive=True)
    kinetic = a / (1 + zeta * k**2 / a**2)
    frequency = s.sqrt(k**2 / a**2 + 1 / zeta)
    G_covariance_symbol = zeta * k**2 * kinetic * frequency / (2 * kappa * a**6)
    assert s.limit(G_covariance_symbol / k, k, s.oo) == 1 / (2 * kappa * a**4)
    assert s.limit(G_covariance_symbol / k**2, k, 0) == s.sqrt(zeta) / (
        2 * kappa * a**5
    )

    # Two regular local COVARIANT finite functionals, evaluated as witnesses:
    # O_N=X^4(X-1)^2;
    # O_G=zeta^2 X^4 [grad_nu(u) nabla_mu F^{mu nu}(W)]^2.
    X, e, x2, j2, n, v, h, G = s.symbols("X e x2 j2 n v h G", real=True)
    polynomial = X**4 * (X - 1) ** 2
    assert [s.diff(polynomial, X, i).subs(X, 1) for i in range(3)] == [0, 0, 2]
    vacuum_N = s.expand(polynomial.subs(X, e**2 * x2))
    vacuum_G = s.expand(zeta**2 * (e**2 * x2) ** 4 * (e**2 * j2) ** 2)
    assert min(power[0] for power, _ in s.Poly(vacuum_N, e).terms()) == 8
    assert min(power[0] for power, _ in s.Poly(vacuum_G, e).terms()) == 12
    for expression in (vacuum_N, vacuum_G):
        assert all(expression.coeff(e, order) == 0 for order in range(5))
    # The loop-graded graph identity places their first possible E4 effects
    # at grades >=3 and >=5, not at the first quantum order being matched.
    first_four_point_grades = [(degree + 2 - 4) // 2 for degree in (8, 12)]
    assert first_four_point_grades == [3, 5]

    N = 1 + e * n
    clock_X = N**-2
    clock_R = 1 + (clock_X - 1) / h
    physical_volume = N * s.exp(3 * e * v) * clock_R ** -s.Rational(3, 4)
    density_N = (
        s.series(physical_volume * polynomial.subs(X, clock_X), e, 0, 3)
        .removeO()
        .expand()
    )
    # At the reference, grad(u).div(F)=G/zeta, from Pi_W=kappa*zeta*a*E_i.
    density_G = (
        s.series(physical_volume * zeta**2 * clock_X**4 * (e * G / zeta) ** 2, e, 0, 3)
        .removeO()
        .expand()
    )
    assert density_N == 4 * e**2 * n**2
    assert density_G == e**2 * G**2
    assert density_N.coeff(e, 1) == density_G.coeff(e, 1) == 0

    # Clock first variations vanish, but the reduced physical Hessians do not.
    J0, F0 = s.Rational(243, 160), s.Rational(1199, 800)
    u = s.Symbol("u", real=True)
    H_clock = 4 * u / (1 + u**2)
    theta_clock = H_clock - u / (1 + u**2) ** 4
    E_clock = 1 - s.Rational(3, 2) / (1 + u**2) ** 3
    ell_clock = 1 / (10 * (1 + u**2) ** 6)
    F_clock = (
        theta_clock * s.diff(E_clock, u)
        - E_clock * s.diff(theta_clock, u)
        + H_clock * E_clock * theta_clock
        - theta_clock**2
        - ell_clock**2 * E_clock**2 / 2
    )
    assert F_clock.subs(u, 0) == F0 == J0 - s.Rational(1, 50)
    row_n, _, _, _ = source_rows(
        time=0, comoving_k_squared=1, lapse_pivot=F(243, 160), profile_current=0
    )
    scalar_direction = 8 * row_n.T * row_n
    response_map = s.diag(scalar_direction[0, 0], 2)
    assert response_map == s.diag(s.Rational(51200, 59049), 2)
    assert response_map.det() == s.Rational(102400, 59049) != 0
    cN, cG, loop, L = s.symbols("cN cG loop L", real=True)
    reduced_H = L**2 / (4 * (J0 + 4 * loop * cN))
    assert s.diff(reduced_H, loop).subs(loop, 0) == -4 * cN * (L / (2 * J0)) ** 2
    # Exact characteristic pencil shows that the scalar witness is not just
    # a disposable c-number or an untransported canonical phase change.
    E0, ell0, speed = s.Rational(-1, 2), s.Rational(1, 10), s.Symbol("speed")
    K = s.Matrix([[2 * (J0 + 4 * loop * cN) / E0**2 + ell0**2, ell0], [ell0, 1]])
    spatial = s.Matrix([[2 * F0 / E0**2 + ell0**2, ell0], [ell0, 1]])
    assert (
        s.factor(
            (speed * K - spatial).det()
            - K.det() * (speed - 1) * (speed - F0 / (J0 + 4 * loop * cN))
        )
        == 0
    )
    speed_sensitivity = s.diff(F0 / (J0 + 4 * loop * cN), loop).subs(loop, 0)
    assert speed_sensitivity == -s.Rational(153472, 59049) * cN
    for comparison in (s.Rational(-1, 10**6), s.Rational(1, 10**6)):
        assert 0 < F0 / (J0 + 4 * comparison) < 1
        assert s.Rational(1, 1000) - 4 * abs(comparison) > 0

    # Independent temporal-vector solve for the second witness. This is a
    # frozen-coefficient, first-order finite-k control, NOT a stationary
    # dispersion law for the evolving clock or a new UV front-velocity claim.
    coordinate, velocity, temporal, momentum_W = s.symbols(
        "coordinate velocity temporal momentum_W", real=True
    )
    zeta_eff = zeta + 2 * loop * cG * zeta**2 * k**2
    lag_vector = (
        zeta_eff * (velocity - k * temporal) ** 2 / 2
        + temporal**2 / 2
        - coordinate**2 / 2
    )
    temporal_star = zeta_eff * k * velocity / (1 + zeta_eff * k**2)
    assert s.cancel(s.diff(lag_vector, temporal).subs(temporal, temporal_star)) == 0
    reduced_lag = s.cancel(lag_vector.subs(temporal, temporal_star))
    velocity_star = (1 / zeta_eff + k**2) * momentum_W
    assert (
        s.cancel(
            s.diff(reduced_lag, velocity).subs(velocity, velocity_star) - momentum_W
        )
        == 0
    )
    reduced_vector_H = s.cancel(
        (momentum_W * velocity - reduced_lag).subs(velocity, velocity_star)
    )
    assert (
        s.simplify(
            s.diff(reduced_vector_H, loop).subs(loop, 0) + cG * k**2 * momentum_W**2
        )
        == 0
    )
    vector_frequency_squared = s.diff(reduced_vector_H, momentum_W, 2)
    assert (
        s.simplify(s.diff(vector_frequency_squared, loop).subs(loop, 0))
        == -2 * cG * k**2
    )

    # Finite reference transport preserves the physical response; resetting
    # the complementary boundary coordinates does not. No coordinate is set.
    eta_N, eta_G, ref_N, ref_G = s.symbols("eta_N eta_G ref_N ref_G", real=True)
    ref, boundary, shift = (
        s.Matrix([ref_N, ref_G]),
        s.Matrix([cN, cG]),
        s.Matrix([eta_N, eta_G]),
    )
    physical = response_map * (ref + boundary)
    assert response_map * ((ref + shift) + (boundary - shift)) == physical
    reset_defect = response_map * (ref + shift) - response_map * ref
    assert reset_defect == response_map * shift
    assert reset_defect != s.zeros(2, 1)

    # Literature O(G) is asymptotic, not a physical constant/radius enclosure.
    # Both controls below are ENTIRE functions K*g; not alleged UV amplitudes.
    _, _, actual_kappa, lam, _ = original_borns()
    g_phys, allocation = F(1, actual_kappa), lam / 100
    constant_limit = actual_kappa * allocation
    assert constant_limit == 10**198
    small_error, large_error = g_phys, (2 * actual_kappa * lam) * g_phys
    assert small_error < allocation and large_error == 2 * lam > allocation
    gravity_certificate_fields = (
        "same_observable_and_low_cut_definition",
        "positive_reference_conversion_bound",
        "finite_angle_or_finite_contour_spectral_bound",
        "finite_coupling_remainder_constant",
        "validity_domain_containing_actual_parameters",
    )

    rejected = 0
    for bad in (-1, 11, True, 2.0):
        try:
            normal_indices(bad)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("Unsupported jet order accepted")
    for terms in (
        {(0, 0, 0): F(1)},
        {(0, 0, 0, -1): F(1)},
        {(0, 0, 0, True): F(1)},
        {(0, 0, 0, 0): 0.5},
        {(0, 0, 0, 0): True},
    ):
        try:
            polynomial_jet_remainder(terms, 6)
        except (TypeError, ValueError):
            rejected += 1
        else:
            raise AssertionError("Malformed exact normal jet accepted")
    assert inputs() == before
    return {
        "milestone": "P8.PRESCRIPTION-2.RESPONSE-EXTENSION-AND-CLOSURE-AUDIT",
        "outcome": "CONTINUUM_EXTENSION_FAMILY_WITH_EXPLICIT_PHYSICAL_INPUT_BLOCKERS",
        "protected_input_files": len(before),
        "physical_model_adopted": False,
        "tracks_complete": False,
        "prescription2_status": "BLOCKED_FOR_UNIQUE_PHYSICAL_COMPLETION_ON_FINITE_BOUNDARY_SPECIFICATION",
        "gravity_status": "CURRENT_ROUTE_BLOCKED_ON_QUANTITATIVE_SAME_OBSERVABLE_PHYSICAL_INPUT",
        "continuum_source_result": "oriented connected scalar/longitudinal noise products exist; retarded source-squared kernels have extensions with local ambiguities; not a chosen covariant physical subtraction",
        "scalar_loop_scaling_bound_before_external_rows": int(scalar_internal),
        "longitudinal_loop_scaling_bound_before_external_rows": int(long_internal),
        "prepared_phase_scaling_bound": int(phase_scalar),
        "normal_derivative_cap": derivative_cap,
        "normal_jet_counts_not_physical_operator_count": jet_counts,
        "covariant_witness_vacuum_field_degrees": [8, 12],
        "covariant_witness_clock_quadratic_densities": ["4*cN*n^2", "cG*G^2"],
        "response_map_determinant": str(response_map.det()),
        "scalar_squared_speed_sensitivity": "-153472*cN/59049 at the formal tree bounce",
        "longitudinal_squared_frequency_sensitivity": "-2*cG*k^2 in the frozen-coefficient control at fixed finite momentum and first quantum order; not a stationary curved-clock dispersion or UV front-velocity result",
        "first_possible_four_point_loop_grades": first_four_point_grades,
        "witness_scope": "different finite quantum boundary completions with the same retained first-order flat data and reference mean; not adopted parents, not all-orders equality or complete healthy UV theories",
        "finite_reference_transport": "boundary shifts oppositely to reference; reset-to-zero changes the physical response",
        "gravity_absolute_error_allocation_control": "For delta=C_err/kappa<lambda/100 one needs C_err<1e198 AND a proved validity domain; no physical C_err is supplied by an O(G) statement",
        "gravity_required_certificate_fields": gravity_certificate_fields,
        "invalid_inputs_rejected": rejected,
        "remaining": "select a physically explicit common quantum boundary prescription or parent and establish same-observable gravity hypotheses/errors; no original M/V/G/B/R/P8 gate closed",
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
