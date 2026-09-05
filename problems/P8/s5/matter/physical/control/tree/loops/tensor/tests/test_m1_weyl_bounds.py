import pytest
import sympy as sp
from p8_m1_weyl import bounds as b
from p8_m1_weyl import reduction as r


def test_all_sufficient_window_energy_and_Duhamel_inequalities():
    assert all(value > 0 for value in b.validate_constants().values())
    assert b.EPS_MAX*b.Q_MAX == sp.Rational(4, 62500)


def test_replay_includes_complete_critical_points_not_only_sampled_bounds():
    extrema = b.compact_extrema()
    assert extrema["friction_F"]["absolute_sup"] == "320"
    assert len(extrema["friction_F"]["stationary_roots"]) == 2
    assert extrema["canonical_G_non_q"]["stationary_roots"] == ["7/24"]
    assert extrema["canonical_G_non_q"]["candidate_values"] == ["0", "-1920", "392"]


def test_original_equation_residual_majorants_match_direct_polynomials():
    assert all(value == 0 for value in b.residual_majorant_checks().values())
    for power, row in r.residual_coefficients().items():
        for component, polynomial in row.items():
            for x in (-1, -sp.Rational(1, 3), 0, sp.Rational(2, 3), 1):
                for q in (1, 17):
                    substitutions = {**dict(zip(r.H_JETS, r.compact_hubble_jets())), r.Q: q}
                    direct = abs(polynomial.subs(substitutions, simultaneous=True).subs(r.X, x))
                    assert direct <= b.residual_majorants()[power][component].subs(r.Q, q)


def test_matching_bound_is_an_input_with_explicit_homogeneous_scaling():
    record = b.estimates(1, b.MTAU_EXAMPLE)
    assert record["relative_energy_norm_evolution_error"] == sp.Rational(10)**(-633)
    assert record["relative_first_Born_remainder"] == sp.Rational(10)**(-1266)
    assert record["relative_original_equation_residual"] == sp.Rational(10)**(-1258)
    larger = b.estimates(2, b.MTAU_EXAMPLE)
    assert larger["relative_energy_norm_evolution_error"] == 2*record["relative_energy_norm_evolution_error"]
    assert larger["relative_first_Born_remainder"] == 4*record["relative_first_Born_remainder"]


@pytest.mark.parametrize("bad", [1.0, sp.Float("0.1"), sp.oo, sp.nan])
def test_bound_interface_rejects_rounded_or_nonfinite_inputs(bad):
    with pytest.raises(ValueError, match="exact finite"):
        b.estimates(bad, b.MTAU_EXAMPLE)


def test_unsafe_coefficient_band_and_zero_mass_fail_closed():
    with pytest.raises(ValueError, match=r"10\^-30"):
        b.estimates(1, 1)
    with pytest.raises(ValueError, match="positive"):
        b.estimates(1, 0)
    assert all(value == 0 for value in b.estimates(0, 1).values())


def test_frequency_enhancement_cannot_be_replaced_by_coefficient_smallness():
    k, epsilon, time = sp.symbols("k epsilon t", real=True)
    solution = sp.sin(k*(1+epsilon)*time)
    leading_change = sp.diff(solution, epsilon).subs(epsilon, 0)
    assert leading_change == k*time*sp.cos(k*time)
    assert sp.simplify(leading_change.subs(time, 1).subs(k, 2*sp.pi)-2*sp.pi) == 0


def test_exponential_enclosure_is_rational_and_rejects_invalid_tail_ratio():
    assert b.exponential_upper(sp.Rational(4, 11)).is_Rational
    with pytest.raises(ValueError):
        b.exponential_upper(20)
