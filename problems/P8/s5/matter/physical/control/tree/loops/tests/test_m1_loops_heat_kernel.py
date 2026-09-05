import sympy as sp
from p8_m1_loops import heat_kernel as h


def test_general_laplace_coefficient_specialization_and_decompositions():
    assert all(value == 0 for value in h.decomposition_checks().values())
    direct = (2*h.RIEM2-2*h.RIC2+5*h.R**2+12*h.BOX_R)/360
    assert sp.expand(h.scalar_a2()-direct) == 0
    assert h.bulk_coefficients() == {"C2": sp.Rational(1, 120), "E4": -sp.Rational(1, 360),
                                     "R2": sp.Rational(1, 72), "BoxR": sp.Rational(1, 30)}


def test_real_determinant_dimensional_pole_and_wick_signs():
    assert all(value == 0 for value in h.convention_checks().values())
    table = h.convention_table()
    euclidean = sp.sympify(table["Euclidean_pole_coefficient_per_A2_over_epsilon"])
    lorentzian = sp.sympify(table["P8_Lorentzian_pole_coefficient_per_A2_over_epsilon"])
    assert sp.simplify(euclidean+lorentzian) == 0
    assert euclidean != lorentzian


def test_conformal_and_r_flat_controls_do_not_erase_weyl():
    controls = h.controls()
    assert controls["conformal_R2_coefficient"] == 0
    assert controls["conformal_C2_coefficient"] == sp.Rational(1, 120)
    assert controls["R_flat_bulk_Weyl_term"] != 0
    assert controls["Ricci_flat_local_is_Euler"] == 0
    # Wrong sign for the endomorphism E changes the conformal R2 coefficient.
    wrong = h.universal_a2().subs({h.E: h.R/6, h.BOX_E: h.BOX_R/6, h.OMEGA2: 0})
    assert sp.expand(wrong).coeff(h.R, 2) == controls["flip_E_sign_at_conformal_coupling_error"]


def test_off_shell_minisuperspace_nonclosure_has_nonzero_fourth_order_symbol():
    witness = h.nonclosure_witness()
    a, ad, add, a3, a4 = (witness[name] for name in ("a", "adot", "addot", "a3", "a4"))
    expected = a*a4+2*ad*a3+3*add**2/2-6*ad**2*add/a+3*ad**4/(2*a**2)
    assert witness["acceleration_hessian"] == a
    assert witness["fourth_derivative_coefficient"] == a
    assert witness["strict_class_acceleration_hessian"] == 0
    assert sp.expand(witness["euler_lagrange"]-expected) == 0


def test_total_derivative_cannot_erase_fourth_order_euler_lagrange_symbol():
    t = sp.Symbol("t")
    a = sp.Function("a")(t)
    boundary = sp.diff(a*sp.diff(a, t)**3, t)
    el_boundary = sum((-1)**j*sp.diff(sp.diff(boundary, sp.diff(a, t, j)), t, j) for j in range(3))
    assert sp.expand(el_boundary) == 0


def test_derivative_metric_redefinition_trades_for_matter_not_free_frame():
    assert all(value == 0 for value in h.field_redefinition_checks().values())
    m, y = sp.symbols("M Y", positive=True)
    traded_on_shell = y**2/(40*m**4)
    assert traded_on_shell != 0
    # This comparison is Einstein plus chi only, not an invocation of the CD EOM.


def test_real_complex_and_euclidean_lorentzian_omission_controls_fire():
    controls = h.controls()
    for key in ("omit_real_scalar_half_residue_error", "omit_R2_acceleration_hessian_error",
                "using_Euclidean_as_Lorentzian_beta_error_per_R2"):
        assert sp.expand(controls[key]) != 0
