"""Rigorous ball corroboration with independent phase/normalization controls."""

import json
from fractions import Fraction

import pytest
from flint import acb, arb, ctx, fmpq
from p8_own_scattering import connection_audit as audit


@pytest.fixture(scope="module")
def replay():
    return audit.report()


@pytest.mark.parametrize("precision", [128, 256])
def test_every_ball_residual_is_finite_contains_zero_and_has_rational_width(replay, precision):
    run = next(item for item in replay["runs"] if item["precision_bits"] == precision)
    assert run["residual_count"] == 46
    assert Fraction(run["residual_absolute_cap"]) == Fraction(1, 2**(precision // 2))
    assert all(all(item.values()) for item in run["residual_checks"].values())
    assert all(run["checks"].values())


def test_report_is_json_native_and_replay_is_deterministic(replay):
    assert json.loads(json.dumps(replay, sort_keys=True)) == replay
    assert audit.report() == replay
    assert all(replay["checks"].values())
    assert all(replay["cross_precision_coefficient_overlap"].values())


def test_rational_B_enclosure_itself_proves_the_positive_margin(replay):
    for run in replay["runs"]:
        lower, upper = map(Fraction, run["outer_rational_enclosures"]["B_abs"])
        assert Fraction(1, 40) < lower < upper < Fraction(1, 30)
        assert lower == Fraction(31349917, 10**9)
        assert upper == Fraction(15674959, 500000000)
        assert run["checks"]["B_strict_lower_1_over_40"]


def test_no_decimal_ball_strings_are_used_as_assertions(replay):
    for run in replay["runs"]:
        for key, enclosure in run["outer_rational_enclosures"].items():
            intervals = [enclosure] if key == "B_abs" else list(enclosure.values())
            for lower, upper in intervals:
                assert "+/-" not in lower + upper
                assert Fraction(lower) <= Fraction(upper)
                assert (Fraction(lower) * audit.ENCLOSURE_DENOMINATOR).denominator == 1
                assert (Fraction(upper) * audit.ENCLOSURE_DENOMINATOR).denominator == 1


@pytest.mark.parametrize("precision", [128, 256])
def test_generic_Gamma_B_and_flux_normalization(precision):
    with ctx.workprec(precision):
        data = audit.balls(precision)
        rho = arb(7).sqrt() / 2
        aa, bb = data["A"], data["B_gamma"]
        assert bb.overlaps(acb(0, 1 / (arb.pi() * rho).sinh()))
        assert (aa * aa.conjugate() - bb * bb.conjugate()).contains(1)
        assert (aa * aa.conjugate()).real > 1
        # A reflected flux is |B/A|^2, not |B|^2.
        difference = abs(bb)**2 - abs(bb / aa)**2
        assert difference > 0
        assert not bb.overlaps(-data["B_closed"])


@pytest.mark.parametrize("precision", [128, 256])
def test_parity_Wronskian_and_phase_are_independently_normalized(precision):
    with ctx.workprec(precision):
        data = audit.balls(precision)
        rho, ee, oo = data["rho"], data["E_plus"], data["O_plus"]
        assert (4 * rho * (ee * oo.conjugate()).imag).contains(1)
        assert not (-4 * rho * (ee * oo.conjugate()).imag).contains(1)
        for fixture in data["radial"].values():
            for row in range(2):
                for col in range(2):
                    assert fixture["parity"][row][col].overlaps(fixture["Jost"][row][col])
        first = data["radial"]["1"]["parity"][0][1]
        assert not first.overlaps(data["A"].conjugate())  # Missing 4*sqrt(2) phase.
        assert not first.overlaps(data["radial"]["1/1000000"]["parity"][0][1])


@pytest.mark.parametrize("x", ["-1", "-1/4", "0", "1/4", "1"])
def test_finite_point_Jost_Cauchy_data_and_both_left_Gauss_branches(replay, x):
    keys = (
        "parity_Wronskian", "Jost_parity_value", "Jost_parity_x_derivative",
        "Jost_time_Wronskian", "left_right_Jost_value", "left_right_Jost_t_derivative",
    )
    for run in replay["runs"]:
        for key in keys:
            assert all(run["residual_checks"][f"x_{x}_{key}"].values())


def test_parity_initial_values_and_Wronskian_sign_directly():
    with ctx.workprec(128):
        data = audit.balls(128)
        point = data["points"]["0"]
        assert point["even"].contains(1)
        assert point["even_x"].contains(0)
        assert point["odd"].contains(0)
        assert point["odd_x"].contains(1)
        assert point["Wronskian_t"].overlaps(-2 * acb(0, data["rho"]))
        assert not point["Wronskian_t"].overlaps(2 * acb(0, data["rho"]))


def test_all_deliberate_mistakes_are_rigorously_separated(replay):
    names = (
        "wrong_B_sign_disjoint", "missing_radial_phase_disjoint",
        "wrong_Wronskian_orientation_disjoint", "reflected_flux_not_B_squared",
        "wrong_unit_Jost_x_amplitude_disjoint",
    )
    for run in replay["runs"]:
        assert all(run["checks"][name] for name in names)


def test_precision_context_is_restored_even_after_report():
    saved = ctx.prec
    try:
        ctx.prec = 73
        audit.balls(128)
        assert ctx.prec == 73
        audit.report()
        assert ctx.prec == 73
    finally:
        ctx.prec = saved


@pytest.mark.parametrize("bad", [True, False, 128.0, Fraction(128), "128", None])
def test_precision_type_guard_is_strict(bad):
    with pytest.raises(TypeError):
        audit.balls(bad)


@pytest.mark.parametrize("bad", [-1, 0, 64, 127, 129, 512])
def test_only_declared_precision_domains_are_admitted(bad):
    with pytest.raises(ValueError):
        audit.balls(bad)


def test_report_does_not_promote_finite_overlap_to_theorem(replay):
    assert replay["role"] == "numerical corroboration, not proof of the analytic connection"
    assert "not a uniform ODE" in replay["boundaries"][1]
    assert "no physical S-matrix" in replay["boundaries"][2]


def test_exact_ball_operations_are_not_plain_floating_tolerances():
    # The actual audit rejects an enclosure that straddles the threshold even
    # if its midpoint would pass a numerical comparison.
    with ctx.workprec(128):
        straddling = arb(fmpq(1, 40), fmpq(1, 10**9))
        assert not straddling > arb(fmpq(1, 40))
        assert straddling.contains(arb(fmpq(1, 40)))
