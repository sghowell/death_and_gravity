import sympy as sp
from p8_matching import escape


def test_nonminimal_escape_solves_both_frames():
    result = escape.derive()
    assert set(result["exact_residuals"].values()) == {"0"}
    assert result["bounce_Planck_jets"] == ["1", "0", "-6"]
    assert result["past_q_limit"] == result["future_q_limit"] == "0"


def test_escape_never_becomes_an_accepted_D_witness():
    result = escape.derive()
    assert "violates P8 tensor tails" in result["warning"]
    assert "no vacuum/UV/perturbation verdict" in result["warning"]
    assert all(value > 0 for value in result["Einstein_kinetic_positive_polynomial_coefficients_in_u_squared"])


def test_escape_Einstein_frame_is_a_turnaround_not_a_bounce():
    result = escape.derive()
    # Parsing only expressions generated locally by this module.
    u = sp.Symbol("u", real=True)
    M, tau = sp.symbols("M tau", positive=True)
    HdE = sp.sympify(result["Einstein_H_dot"], locals={"u": u, "M": M, "tau": tau})
    assert HdE.subs(u, 0) == -1/tau**2
