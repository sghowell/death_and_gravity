"""Independent full-source, physical-helicity, angular, spectral and scope tests."""

import pytest
import sympy as s
from p8_vacuum_affine_massive_graviton_cut_isolated_replay import (
    amplitude,
    audit,
    sewing,
    source,
    window,
)

ROWS = audit.residuals()
GATES = audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_identity(name):
    value = ROWS[name]
    assert all(
        x == 0 for x in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not value.atoms(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[x[0] for x in audit.bad_cases()]
)
def test_all_unsupported_inputs(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("spin", range(13))
def test_spin_validation_precedes_cached_values(spin):
    assert sewing.legendre_q(spin) == sewing.legendre_q(s.Integer(spin))
    for bad in (bool(spin), float(spin), s.Float(spin), str(spin)):
        with pytest.raises((TypeError, ValueError)):
            sewing.legendre_q(bad)


@pytest.mark.parametrize("sign1,sign2", ((1, 1), (1, -1), (-1, 1), (-1, -1)))
@pytest.mark.parametrize(
    "energy,momentum,cosine",
    ((5, 3, s.Rational(3, 5)), (13, 5, s.Integer(0)), (17, 8, s.Rational(-4, 5))),
)
def test_complete_four_graph_amplitude_at_independent_massive_points(
    sign1, sign2, energy, momentum, cosine
):
    label = (
        ("plus" if sign1 == 1 else "minus") + "_" + ("plus" if sign2 == 1 else "minus")
    )
    row = amplitude.data()["literal_four_graph_helicity_rows"][label]
    point = {amplitude.E: energy, amplitude.P: momentum, amplitude.X: cosine}
    direct = sum(
        row[k]
        for k in (
            "contact",
            "first_scalar_exchange",
            "second_scalar_exchange",
            "entire_EH_exchange",
        )
    ).subs(point)
    mass2 = energy * energy - momentum * momentum
    den = energy * energy - momentum * momentum * cosine * cosine
    wanted = (
        mass2 * mass2 if sign1 == sign2 else momentum**4 * (1 - cosine * cosine) ** 2
    ) / den
    assert s.factor(direct - wanted) == 0


@pytest.mark.parametrize("mass2", (s.Rational(1, 3), s.Integer(1), s.Integer(7)))
def test_massive_threshold_not_replaced_by_massless_scalar(mass2):
    same, opposite = amplitude.helicity_trees(
        4 * mass2, s.Symbol("x"), mass2, s.Integer(5)
    )
    assert s.factor(same - mass2 / 5) == 0
    assert opposite == 0
    assert same != 0


def test_full_TT_projector_is_idempotent_and_rank_two():
    P = amplitude.data()["canonical_TT_projector"]
    assert P * P == P and s.trace(P) == 2 and P.rank() == 2
    assert P.T == P
    trace = s.Matrix(
        [1 if i // 4 == i % 4 and i // 4 in (1, 2) else 0 for i in range(16)]
    )
    assert P * trace == s.zeros(16, 1)


def test_opposite_helicity_phase_cannot_be_replaced_by_modulus():
    a = s.Rational(3, 5)
    b = s.Rational(4, 5)
    assert s.re((a + s.I * b) ** 4) == -s.Rational(527, 625)
    assert (a * a + b * b) ** 2 == 1
    # A negative nonforward angular integrand is compatible with positive spin weights.
    value = sewing.full_helicity_sewing(
        s.Integer(100), s.Integer(1), s.Integer(1), s.Integer(0), s.Integer(0), a
    )
    assert value < 0


@pytest.mark.parametrize(
    "r", (s.Rational(1, 9), s.Rational(1, 4), s.Rational(9, 16), s.Rational(25, 36))
)
@pytest.mark.parametrize(
    "y", (s.Integer(0), s.Rational(1, 5), s.Rational(4, 9), s.Integer(1))
)
def test_full_pointwise_helicity_bound(r, y):
    A, B = 1 - r, r * (1 - y)
    den = A + B
    u, v = A * A / den, B * B / den
    assert 0 <= u + v <= den <= 1 and 0 < u * u + v * v <= 1


@pytest.mark.parametrize("beta", (s.Rational(1, 3), s.Rational(3, 5), s.Rational(4, 5)))
def test_forward_integral_by_independent_exact_rational_integration(beta):
    x = s.Symbol("x", real=True)
    r = beta * beta
    actual = (
        s.integrate(
            ((1 - r) ** 4 + r**4 * (1 - x * x) ** 4) / (1 - r * x * x) ** 2, (x, -1, 1)
        )
        / 2
    )
    expected = sewing.forward_shape(beta).rewrite(s.log)
    assert s.simplify(s.expand_log(actual - expected, force=True)) == 0


@pytest.mark.parametrize("spin", (0, 2, 4, 6, 8))
@pytest.mark.parametrize("w", (s.Integer(2), s.Integer(3)))
def test_Q_recurrence_against_independent_fixed_real_sheet_integral(spin, w):
    x = s.Symbol("x", real=True)
    actual = s.integrate(s.legendre(spin, x) / (w - x), (x, -1, 1)) / 2
    expected = sewing.legendre_q(spin, w)
    assert s.simplify(s.expand_log(actual - expected, force=True)) == 0


@pytest.mark.parametrize("spin", (0, 2, 4, 6, 8, 10, 12))
def test_every_sampled_even_weight_is_a_nonnegative_real_square(spin):
    value = sewing.spin_weight(spin).subs({sewing.W: 2, source.S: 13, source.K: 7})
    assert value.is_real is True and value.is_nonnegative is True
    f0, f4, norm2 = sewing.spin_components(spin)
    assert (f4 == 0 and norm2 == 0) if spin < 4 else norm2 > 0
    assert f0 != 0


@pytest.mark.parametrize("spin", (1, 3, 5, 7, 9, 11))
def test_odd_spin_weights_vanish(spin):
    assert sewing.spin_weight(spin) == 0


@pytest.mark.parametrize("pair", ((4, 6), (4, 8), (6, 8)))
def test_independent_spin4_cross_orthogonality(pair):
    a, b = pair
    x = s.Symbol("x", real=True)
    product = (
        (1 - x * x) ** 4
        * s.diff(s.legendre(a, x), x, 4)
        * s.diff(s.legendre(b, x), x, 4)
    )
    assert s.integrate(product, (x, -1, 1)) == 0


def test_independent_low_velocity_Parseval_includes_spin4():
    r, x = s.symbols("r x", real=True)
    f0 = s.expand((1 - r) ** 2 * sum(r**j * x ** (2 * j) for j in range(5)))
    f4 = s.expand(r * r * (1 - x * x) ** 2 * sum(r**j * x ** (2 * j) for j in range(3)))
    total = 0
    for ell in range(0, 9, 2):
        a0 = (2 * ell + 1) * s.integrate(f0 * s.legendre(ell, x), (x, -1, 1)) / 2
        total += a0 * a0 / (2 * ell + 1)
        if ell >= 4:
            norm2 = s.factorial(ell - 4) / s.factorial(ell + 4)
            a4 = (
                (2 * ell + 1)
                * s.integrate(
                    f4 * (1 - x * x) ** 2 * s.diff(s.legendre(ell, x), x, 4), (x, -1, 1)
                )
                / 2
            )
            total += norm2 * a4 * a4 / (2 * ell + 1)
    expected = (
        1
        - s.Rational(10, 3) * r
        + s.Rational(59, 15) * r * r
        - s.Rational(64, 35) * r**3
        + s.Rational(64, 105) * r**4
    )
    assert s.expand(s.series(total, r, 0, 5).removeO() - expected) == 0
    assert (
        s.expand(
            s.series(sewing.forward_shape(), sewing.BETA, 0, 10).removeO()
            - expected.subs(r, sewing.BETA**2)
        )
        == 0
    )


@pytest.mark.parametrize(
    "lower,upper,mass2,kappa",
    (
        (4, 5, 1, 1),
        (8, 17, 2, 7),
        (s.Rational(4, 3), 2, s.Rational(1, 3), 11),
        (4, 10**196, 1, 10**800),
    ),
)
def test_exact_physical_window_and_all_arguments_validated(lower, upper, mass2, kappa):
    assert window.require_window(lower, upper, mass2, kappa) == (
        lower,
        upper,
        mass2,
        kappa,
    )
    for i in range(4):
        args = [lower, upper, mass2, kappa]
        args[i] = True
        with pytest.raises((TypeError, ValueError)):
            window.require_window(*args)


def test_window_majorant_is_an_integral_not_a_subtraction_constant():
    z, mu = s.symbols("z mu", positive=True)
    derivative = s.diff(window.moment_primitive(z, mu), z)
    assert s.factor(derivative - z * z / (z - 2 * mu) ** 3) == 0
    assert "NOT the complete b20" in window.data()["not_established"]
    assert "other intermediate states" in window.data()["not_established"]


def test_named_window_bound_is_small_only_relative_to_the_old_reference():
    d = window.data()["named_window"]
    assert d["comparison_to_unchanged_reference_4lambda"] < s.Rational(1, 10**990)
    assert d["strict_upper_bound"] == s.Rational(591, 128 * 10**1600)
    assert "bare term canceled" in window.data()["window_proof"]


def test_entire_old_action_and_all_three_fixed_constants_return_at_marker_one():
    d = source.data()
    R, F = d["entire_original_coefficient_functions"]
    marked = d["entire_formal_loop_marked_scalar_coefficient"]
    assert s.expand(marked.subs(source.LOOP, 1) - F) == 0
    assert R == source.current.heavy.coefficients()["R"]
    assert (
        d["all_three_retained_vacuum_constants"]
        == source.current.fixed_functions()["vacuum_pressure_total"]
    )
    assert marked.has(source.LOOP)


def test_all_previous_rows_and_qualifications_are_unchanged():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 136
    assert audit.matching()[:-2] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert audit.validate_scope(audit.frontier(), audit.matching())
    assert "unphysical" in audit.observable()["not_established"].lower()


def test_selected_predecessor_contract_does_not_alias_parent_dictionary():
    full, selected = source.retained_vacuum_checks()
    parent_full, parent_rows = source.previous.full_vacuum()
    assert full is parent_full
    assert selected is not parent_rows
    assert tuple(selected) == source.VACUUM_CHECK_NAMES
    selected.pop(source.VACUUM_CHECK_NAMES[0])
    assert source.VACUUM_CHECK_NAMES[0] in parent_rows


def test_parent_packet_warmup_cannot_add_checks_to_this_contract():
    full, before = source.retained_vacuum_checks()
    source.previous.data()
    after_full, after = source.retained_vacuum_checks()
    assert after_full is full and after == before
    assert len(after) == 10
    assert len(source.previous.full_vacuum()[1]) == 19


def test_entire_source_packet_identical_after_parent_warmup_and_own_cache_reset():
    before = source.data()
    source.previous.data()
    source.data.cache_clear()
    after = source.data()
    assert before == after
    assert len(after["checks"]) == 21


def test_rejected_native_packet_is_pinned_without_adopting_its_verdict():
    from p8_vacuum_affine_massive_graviton_cut_isolated_replay import verify

    data = verify.rejected_packet_evidence()
    assert data["status"] == "REJECTED_NOT_AN_ACCEPTED_PREDECESSOR"
    assert data["native_named_and_scalar"] == [139, 409]
    assert data["cold_preflight_named_and_scalar"] == [130, 391]
    assert data["unchanged_source_files"] == 18


def test_rejected_and_corrected_matching_records_are_both_explicit():
    assert audit.matching()[-2] == audit.REJECTED_ITEM
    assert audit.matching()[-1] == audit.ITEM
    assert audit.REJECTED_ITEM["status"].startswith("REJECTED_")
    assert audit.matching()[:-2] == audit.previous.matching()
