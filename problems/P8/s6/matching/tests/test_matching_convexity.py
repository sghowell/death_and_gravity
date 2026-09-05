import sympy as sp
from p8_matching import convexity as c


def test_stationary_and_target_identities():
    result = c.derive()
    assert set(result["exact_residuals"].values()) == {"0"}
    assert result["target_bounce_F_Fx_Fxx"] == ["-4", "10", "-10"]
    assert result["target_minus_Fxx_sign"]["numerator"]["roots"] == 0


def test_two_heavy_positive_and_zero_coupling_controls():
    j1, j2, k11, k12, k22 = sp.symbols("j1 j2 k11 k12 k22", real=True)
    curvature = c.stationary_checks()["two_heavy_Pxx"]
    assert curvature.subs({j1: 2, j2: -1, k11: 2, k12: 1, k22: 3}) == sp.Rational(18, 5)
    assert curvature.subs({j1: 0, j2: 0}) == 0


def test_indefinite_stationary_hessian_is_not_positive_branch():
    j1, j2, k11, k12, k22 = sp.symbols("j1 j2 k11 k12 k22", real=True)
    curvature = c.stationary_checks()["two_heavy_Pxx"]
    assert curvature.subs({j1: 1, j2: 0, k11: -1, k12: 0, k22: 1}) == -1


def test_nonlinear_stationary_gap_and_sign():
    X, mu, b, source = sp.symbols("X mu b j", real=True)
    result = c.stationary_checks()
    assert result["nonlinear_Pxx"].subs({X: 1, mu: 2, b: 1, source: 3}) == sp.Rational(16, 3)
    assert result["nonlinear_Hessian"].subs({X: 4, mu: 2, b: 1}) == 0


def test_mixed_gap_does_not_imply_stationary_stability():
    result = c.mixing_control()
    assert result["negative_control"]["gap_squared"] == 1
    assert result["negative_control"]["lower_root_negative"] is True
    assert result["residuals"]["bare_velocity_Hessian"] == sp.zeros(2)
    omega, k2, kappa, mixing = sp.symbols("omega k2 kappa c", real=True)
    y = sp.Symbol("y", real=True)
    polynomial = result["determinant"].subs({kappa: -1, mixing: sp.sqrt(2), k2: 2, omega**2: y})
    assert all(root > 0 for root in sp.solve(polynomial, y))


def test_physical_second_jet_error_has_no_M_tau_suppression():
    result = c.target_checks()
    M, tau = sp.symbols("M tau", positive=True)
    assert sp.cancel(result["physical_remainder_second_jet_magnitude_lower_bound"]*tau**2/M**2) == 10
