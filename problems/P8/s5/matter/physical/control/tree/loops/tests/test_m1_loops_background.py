import pytest
import sympy as sp
from p8_m1_loops import background as b


def test_cosmic_compact_variation_trace_and_conservation_identities():
    assert all(value == 0 for value in b.algebra_checks().values())


def test_independent_hubble_contractions_at_rational_times_and_two_scales():
    t = sp.Symbol("t", real=True)
    for tau in (sp.Rational(2, 3), sp.Integer(7)):
        a = (1+(t/tau)**2)**2
        hubble = sp.diff(a, t)/a
        hdot = sp.diff(hubble, t)
        r00 = -3*(hdot+hubble**2)
        rp = hdot+3*hubble**2
        ricci = r00-3*rp
        ricci_squared = r00**2+3*rp**2
        riemann_squared = 12*((sp.diff(a, t, 2)/a)**2+hubble**4)
        for u in (-3, -sp.Rational(1, 2), 0, sp.Rational(2, 3), 5):
            substitutions = {b.U: u, b.TAU: tau}
            for key, expected in (("R", ricci), ("Ric2", ricci_squared), ("Riem2", riemann_squared)):
                assert sp.simplify(b.cosmic()[key].subs(substitutions)-expected.subs(t, tau*u)) == 0


def test_exact_all_time_extrema_include_interior_density_minimum():
    bounds = b.bounds()
    assert bounds["I_R2_00"]["min"] == -224
    assert bounds["I_R2_00"]["absolute_sup"] == 6048
    assert sp.Rational(2, 9) in bounds["I_R2_00"]["candidate_v"]
    assert bounds["I_R2_pressure"]["min"] == -4032
    assert bounds["I_R2_pressure"]["max"] == 576
    assert bounds["I_A2_00"]["absolute_sup"] == 84
    assert bounds["I_A2_pressure"]["absolute_sup"] == 56
    assert bounds["BoxR"]["min"] == -sp.Rational(728, 3)
    assert bounds["BoxR"]["absolute_sup"] == 3024
    assert bounds["A2_chosen_representative"]["absolute_sup"] == 480
    assert bounds["R"]["absolute_sup"] == 168


def test_extrema_checker_rejects_unproved_domains_and_degrees():
    with pytest.raises(ValueError, match="even"):
        b.polynomial_range(b.X)
    with pytest.raises(ValueError, match="quadratics"):
        b.polynomial_range(b.X**6)


def test_conditional_log_bound_restores_M_tau_hbar_and_never_divides_by_G00():
    mtau = sp.Symbol("M_tau", positive=True)
    log, hbar = sp.symbols("L hbar", nonnegative=True)
    exact = b.constant_log_source_bound(mtau, log, hbar, rational=False)
    assert sp.simplify(exact-2*6048*hbar*log/(72*(4*sp.pi)**2*mtau**2)) == 0
    assert b.constant_log_source_bound(mtau, log, hbar) == 7*hbar*log/(6*mtau**2)
    assert b.constant_log_source_bound(b.MTAU_S58, 1) == sp.Rational(7, 6)*sp.Rational(10)**(-648)
    assert 168/sp.Integer(12)**2 == sp.Rational(7, 6)  # pi>3 supplies the upper bound.
    for args in ((0, 1), (1, -1)):
        with pytest.raises(ValueError):
            b.constant_log_source_bound(*args)


def test_finite_matching_coefficients_are_not_bounded_by_the_log_coefficient():
    finite = sp.Symbol("C", nonnegative=True)
    bound = b.finite_local_source_bound(b.MTAU_S58, finite)
    assert sp.diff(bound, finite) != 0
    assert b.constant_log_source_bound(1, 0) == 0
    assert b.finite_local_source_bound(1, 1) == 12096
    with pytest.raises(ValueError):
        b.finite_local_source_bound(1, -1)


@pytest.mark.parametrize("bad", [1.0, sp.Float("0.1"), sp.oo, -sp.oo, sp.zoo, sp.nan,
                                 sp.Symbol("unproved_finite")])
def test_bound_interface_rejects_rounded_nonfinite_and_unproved_inputs(bad):
    with pytest.raises(ValueError, match="exact finite"):
        b.constant_log_source_bound(bad, 1)
    with pytest.raises(ValueError, match="exact finite"):
        b.constant_log_source_bound(1, bad)
    with pytest.raises(ValueError, match="exact finite"):
        b.constant_log_source_bound(1, 1, bad)
    with pytest.raises(ValueError, match="exact finite"):
        b.finite_local_source_bound(1, bad)


def test_bound_interface_accepts_exact_algebraic_inputs_without_rounding():
    assert b.constant_log_source_bound(sp.sqrt(2), sp.Rational(1, 3)) == sp.Rational(7, 36)
    assert b.finite_local_source_bound(sp.sqrt(2), sp.Rational(1, 3)) == 2016


def test_R_flat_variation_and_misleading_pointwise_density_controls():
    controls = b.controls()
    assert all(controls[key] == 0 for key in ("R_flat_radiation_R", "R_flat_radiation_I_R2_00",
                                             "R_flat_radiation_I_R2_pressure", "full_A2_density_bounce"))
    assert controls["bulk_A2_variation_bounce_not_zero"] == 4
    for key in ("drop_curvature_derivatives_I00_error_at_x_half", "wrong_R_sign_curvature_error_at_bounce",
                "drop_ell_derivative_Rdot_error", "finite_local_coefficient_remains_free"):
        assert sp.expand(controls[key]) != 0


def test_flat_tensor_second_variation_survives_zero_background_weyl():
    assert b.compact()["C2"] == 0
    assert all(value == 0 for value in b.tensor_principal_checks().values())
    # A timelike Fourier TT trial is a coefficient control, not an on-shell mode.
    tensor_polynomial = sp.Rational(1, 2)*(sp.Integer(3)**2-sp.Integer(2)**2)**2
    assert tensor_polynomial == sp.Rational(25, 2)
