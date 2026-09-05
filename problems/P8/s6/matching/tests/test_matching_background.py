import sympy as sp
from p8_matching import background as b


def test_metric_ricci_and_time_boundary():
    t, a, N = b.t, b.a, b.N
    H = sp.diff(a, t)/a
    cosmic_R = b.ricci_from_metric().subs({N: 1, sp.diff(N, t): 0})
    assert sp.simplify(cosmic_R+6*(sp.diff(H, t)+2*H**2)) == 0
    assert set(b.derive()["exact_residuals"].values()) == {"0"}


def test_minimal_raychaudhuri_is_potential_independent():
    M = sp.Symbol("M", positive=True)
    ray = b.equations()["Raychaudhuri_equation"]
    minimal = ray.subs({b.Q: M**2, sp.diff(b.Q, b.t): 0, sp.diff(b.Q, b.t, 2): 0})
    assert sp.simplify(minimal-2*M**2*sp.diff(sp.diff(b.a, b.t)/b.a, b.t)-b.kinetic) == 0
    assert not minimal.has(b.potential)


def test_target_requires_negative_null_stress():
    M, tau = sp.symbols("M tau", positive=True)
    target = b.equations()["target_null_stress"]
    assert target.subs(b.t, 0) == -4*M**2/tau**2
    assert b.equations()["bounce_stress_remainder_minimum"] == 4*M**2/tau**2


def test_mixed_field_metric_stays_positive_quadratic_form():
    residual, squares = b.field_metric_control()
    assert residual == 0
    symbols = {str(symbol): symbol for symbol in squares.free_symbols}
    point = {symbols[name]: value for name, value in
             {"v1": 1, "v2": -2, "v3": 3, "d1": 2, "d2": 3, "d3": 5,
              "l21": 1, "l31": -1, "l32": 2}.items()}
    assert squares.subs(point) > 0
