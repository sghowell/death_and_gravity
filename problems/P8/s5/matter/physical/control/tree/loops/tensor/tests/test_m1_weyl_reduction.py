import sympy as sp
from p8_m1_weyl import reduction as r


def test_off_shell_fourth_order_factorization_and_baseline_substitution():
    assert all(value == 0 for value in r.operator_checks().values())


def test_off_shell_action_map_and_physical_branch_map_are_distinct():
    assert all(value == 0 for value in r.field_map_checks().values())
    assert r.controls()["drop_off_shell_field_map_E0_term_error"] != 0
    assert r.controls()["naive_y_speed_used_as_physical_at_bounce_error"] == 64


def test_full_old_tensor_normalization_retains_volume_mass_correction():
    assert all(value == 0 for value in r.canonical_checks().values())
    assert r.controls()["omit_canonical_volume_mass_error_at_x_half"] != 0


def test_exact_fourth_order_residual_starts_at_beta_squared():
    assert all(value == 0 for value in r.residual_checks().values())
    assert r.residual_coefficients()[0]["velocity"] != 0


def test_compact_Chebyshev_derivative_jets_and_coefficients():
    assert all(value == 0 for value in r.compact_checks().values())
    for n in range(5):
        assert r.compact_hubble_jets()[n].subs(r.X, 1) == 4*(-1)**n*sp.factorial(n)


def test_independent_cosmic_CD_jets_at_rational_times():
    t, tau = sp.symbols("t tau", positive=True)
    u = t/tau
    ell = tau*sp.sqrt(1+u**2)
    h = 4*u/(tau*(1+u**2))
    for n, compact in enumerate(r.compact_hubble_jets()):
        assert sp.simplify(sp.diff(h, t, n)-compact.subs(r.X, u/sp.sqrt(1+u**2))/ell**(n+1)) == 0


def test_sign_changing_physical_two_derivative_speed_is_not_a_front_claim():
    controls = r.controls()
    assert controls["physical_speed_shift_at_bounce_per_beta_over_tau2"] == -64
    assert controls["physical_speed_shift_at_u_sqrt3_per_beta_over_ell2"] == 32
    assert controls["speed_shift_zero_at_u_one"] == 0
    assert controls["flat_background_reduced_correction"] == 0


def test_exact_compact_friction_and_mass_extrema():
    x = r.X
    # The cubic B has stationary points +/-sqrt(7)/6 and endpoints +/-1.
    cubic = r.compact_coefficients()["B"]
    assert cubic.subs(x, 1) == 40
    assert abs(cubic.subs(x, sp.sqrt(7)/6)) < 40
    assert sp.diff(cubic, x).subs(x, sp.sqrt(7)/6) == 0
    v = sp.Symbol("v")
    mass = 2688*v-4608*v**2
    assert mass.subs(v, sp.Rational(7, 24)) == 392
    assert mass.subs(v, 1) == -1920


def test_flat_resummed_extra_branch_is_not_kept_in_reduced_equation():
    omega, momentum, beta = sp.symbols("omega k beta")
    d = momentum**2-omega**2
    fourth = sp.expand(d-4*beta*d**2)
    assert sp.factor(fourth) == sp.factor(d*(1-4*beta*d))
    assert sp.expand(fourth.subs(omega**2, momentum**2-1/(4*beta))) == 0
    assert d.subs(omega**2, momentum**2-1/(4*beta)) != 0
