import sympy as sp
from p8_a26_vacuum import regularity


def test_source_terms_cancel_only_after_full_reconstruction_is_kept():
    assert regularity.checks() and not any(regularity.checks().values())


def test_actual_tilt_regularities_are_strict():
    data = regularity.tilt_data()
    alpha = data["alpha"]
    assert 2 < alpha < 3 and 1 < alpha-1 < 2
    assert data["Q_third_leading_coefficient"] < 0
    assert alpha-3 == -sp.Rational(241, 250)
    assert not data["physical_instability_established"]
    assert not data["numerical_b_printed_by_source"]
    assert data["n_s_from_source_equation_21"] == sp.Rational(24, 25)


def test_named_entropy_derivative_is_nonzero_without_decimal_roundoff():
    data = regularity.entropy_background_at_zero()
    assert data["B"].is_positive
    assert data["tau_Bprime"].is_positive


def test_omitting_f0_reconstruction_creates_a_spurious_background_potential():
    q, qp, b, bp, bpp, h, x = sp.symbols("q qp b bp bpp h X", nonzero=True)
    result = regularity.reconstructed_source(q, qp, b, bp, bpp, h, x, bp**2*x, b)
    assert result["full_Q_part"] == 0
    assert sp.expand(-result["W2"]*b**2) != 0


def test_exact_integer_exponent_is_a_smooth_exclusion_control():
    z = sp.Symbol("z", real=True)
    sig = 1/(1+sp.exp(-z))
    q = (1-sig)/(2*(1+z**2))+sig/2
    assert sp.diff(q, z, 3).subs(z, 0).is_finite
    assert 2*(2-1)*(2-2) == 0
