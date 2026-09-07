"""Coefficient-16 connection, endpoint normalization, and rational budgets."""

from fractions import Fraction

import mpmath as mp
import pytest
import sympy as sp
from p8_own_scattering import connection as c


def test_exact_identity_replay():
    assert c.identities()
    assert all(value == 0 for value in c.identities().values())


def test_all_rational_margins():
    assert all(c.checks().values())


def test_even_odd_initial_normalization():
    data = c.even_odd()
    x = data["x"]
    assert sp.simplify(data["even"].subs(x, 0)) == 1
    assert sp.simplify(sp.diff(data["even"], x).subs(x, 0)) == 0
    assert sp.simplify(data["odd"].subs(x, 0)) == 0
    assert sp.simplify(sp.diff(data["odd"], x).subs(x, 0)) == 1


def test_coefficient_80_is_a_real_omission_control():
    # A/(1+8x^2) becomes (A/8-1/4)+(3/4)sech^2(t).
    assert c.RHO**2 == sp.Rational(7, 4)
    assert sp.Rational(80, 8) - sp.Rational(1, 4) - c.RHO**2 == 8


def test_reference_transfer_numeric_parity_crosscheck():
    # Corroboration only: the certificate uses the exact parameter identities
    # and source-audited Gamma connection, not numerical integration.
    with mp.workdps(60):
        rho = mp.sqrt(7) / 2
        aa = mp.gamma(1 - 1j * rho) * mp.gamma(-1j * rho) / (
            mp.gamma(mp.mpf("1.5") - 1j * rho)
            * mp.gamma(-mp.mpf(".5") - 1j * rho)
        )
        bb = 1j / mp.sinh(mp.pi * rho)
        ee = mp.sqrt(mp.pi) * mp.gamma(1j * rho) * 8**(mp.mpf(".25") + 1j * rho / 2) / (
            mp.gamma(-mp.mpf(".25") + 1j * rho / 2)
            * mp.gamma(mp.mpf(".75") + 1j * rho / 2)
        )
        oo = mp.gamma(mp.mpf("1.5")) * mp.gamma(1j * rho) * 8**(-mp.mpf(".25") + 1j * rho / 2) / (
            mp.gamma(mp.mpf(".25") + 1j * rho / 2)
            * mp.gamma(mp.mpf("1.25") + 1j * rho / 2)
        )
        parity_basis = mp.matrix([[ee, oo], [mp.conj(ee), mp.conj(oo)]])
        transfer = parity_basis * mp.matrix([[1, 0], [0, -1]]) * parity_basis**-1
        phase = (4 * mp.sqrt(2))**(2j * rho)
        expected = mp.matrix([[-mp.conj(bb), mp.conj(aa) * phase], [aa / phase, -bb]])
        assert mp.norm(transfer - expected) < mp.mpf("1e-55")
        assert abs(4 * rho * mp.im(ee * mp.conj(oo)) - 1) < mp.mpf("1e-55")
        assert abs(abs(aa)**2 - abs(bb)**2 - 1) < mp.mpf("1e-55")


def test_nonzero_reflection_is_not_reflected_transfer_modulus():
    data = c.coefficients()
    assert data["B_abs_squared"].is_positive is True
    assert sp.simplify(data["B_abs_squared"] / data["A_abs_squared"] - data["reflected_flux"]) == 0
    assert sp.simplify(data["transmitted_flux"] + data["reflected_flux"]) == 1


def test_two_kinds_of_determinant_and_real_radial_basis():
    ar, ai, bb = sp.symbols("ar ai bb", real=True)
    aa = ar + sp.I * ai
    time = sp.Matrix([[sp.conjugate(aa), sp.I * bb], [-sp.I * bb, aa]])
    radial = sp.Matrix([[sp.I * bb, sp.conjugate(aa)], [aa, -sp.I * bb]])
    assert sp.expand(time.det() + radial.det()) == 0
    real = c.unitary_power_coordinates().H * radial * c.unitary_power_coordinates()
    assert sp.simplify(real - sp.Matrix([[ar, ai + bb], [ai - bb, -ar]])) == sp.zeros(2)


def test_physical_delta_phase_must_not_be_dropped():
    delta = sp.symbols("delta", positive=True)
    phase = sp.exp(2 * sp.I * c.RHO * sp.log(c.PHASE_SCALE / sp.sqrt(delta)))
    assert sp.simplify(sp.diff(phase, delta) / phase + sp.I * c.RHO / delta) == 0
    assert c.power_transfer(Fraction(1, 10**6)).shape == (2, 2)


@pytest.mark.parametrize("u", [Fraction(-1, 100), Fraction(1, 100)])
def test_actual_endpoint_map_from_direct_chain_rule(u):
    delta = sp.Rational(1, 10**6)
    uu = sp.symbols("u", real=True)
    q = sp.Function("Q")(uu)
    k = 3 + (uu + 1)**2
    radius = sp.sqrt(uu**2 + delta / 8)
    psi = sp.sqrt(k) * q / sp.sqrt(radius)
    direct = sp.Matrix([psi, radius * sp.diff(psi, uu) / c.RHO])
    actual = c.endpoint_map(uu, delta, k, sp.diff(k, uu)) * sp.Matrix([q, sp.diff(q, uu)])
    assert all(sp.simplify(value.subs(uu, u)) == 0 for value in direct - actual)


def test_kinetic_derivative_omission_changes_actual_clock_readout():
    full = c.endpoint_map(Fraction(1, 100), Fraction(1, 10**6), 4, 1)
    omitted = c.endpoint_map(Fraction(1, 100), Fraction(1, 10**6), 4, 0)
    assert (full - omitted)[1, 0].is_positive is True
    assert (full - omitted)[0, 0] == 0


def test_wave_coefficient_extraction_has_no_hidden_factor():
    t = sp.symbols("t", real=True)
    frame = c.plane_wave_frame(t)
    assert sp.simplify(frame.H * frame) == sp.eye(2)
    assert sp.simplify(frame * frame.H) == sp.eye(2)


def test_energy_factor_is_one_half_not_a_transplanted_mu_bound():
    potential = sp.symbols("V", positive=True)
    generator = sp.Matrix([[0, c.RHO], [-c.RHO - potential / c.RHO, 0]])
    symmetric = (generator + generator.T) / 2
    assert symmetric.eigenvals() == {-potential / (2 * c.RHO): 1, potential / (2 * c.RHO): 1}
    assert c.RHO < 3  # The old response's mu>3 shortcut is inapplicable.


def test_independent_fraction_exponential_majorants():
    f = Fraction
    # e tail begins at n=3; successive ratios are <=1/4.
    e_upper = 1 + 1 + f(1, 2) + f(1, 6) / (1 - f(1, 4))
    exp_small_upper = 1 + f(3, 5) + f(3, 5)**2 / 2 + f(3, 5)**3 / 6 / (1 - f(3, 20))
    assert e_upper == f(49, 18) < f(11, 4)
    assert exp_small_upper == f(1549, 850) < 2
    assert f(11, 4)**4 / (1 - f(1, 4)) == f(14641, 192) < 80


def test_independent_fraction_full_window_budget():
    f = Fraction
    length, delta, remainder = f(1, 100), f(1, 10**6), f(44)
    radius_upper = f(1, 2000)
    assert delta / 8 < radius_upper**2
    central = remainder * (length**2 + 2 * length * radius_upper)
    tails = 3 * delta / (32 * length**2)
    budget = f(4, 5) * (central + tails)
    assert (central, tails, budget) == (f(121, 25000), f(3, 3200), f(2311, 500000))
    assert budget < f(1, 200)
    error = 4 * (1 / (1 - f(1, 400)) - 1)
    assert error == f(4, 399) < f(1, 80)
    assert f(1, 40) - error > f(1, 80)


def test_broader_delta_domain_does_not_pass_this_budget():
    # The exact branch's .01 domain is NOT the fixed-window scattering gate.
    length = Fraction(1, 100)
    tails = 3 * Fraction(1, 100) / (32 * length**2)
    assert tails * Fraction(4, 5) > Fraction(1, 200)


def test_scope_is_conditional_and_homogeneous():
    data = c.calibration()
    assert data["conditional_premise"].startswith("actual canonical remainder")
    assert "no physical vacuum" in data["exterior_scope"]
    assert "not a zero-data" in data["data_scope"]


@pytest.mark.parametrize("bad", [True, False, sp.true, sp.false, 1.0, sp.Float("0.01"), sp.oo, sp.nan, -1, 0])
def test_delta_guard_rejects_invalid_exact_domain(bad):
    c.power_transfer(sp.Rational(1, 10**6))  # Warm unrelated valid calls first.
    with pytest.raises((TypeError, ValueError)):
        c.power_transfer(bad)


@pytest.mark.parametrize("bad", [True, 1.0, sp.oo, sp.nan, -1, 0])
def test_kinetic_guard(bad):
    with pytest.raises((TypeError, ValueError)):
        c.endpoint_map(0, Fraction(1, 10**6), bad, 0)


def test_real_clock_guard():
    with pytest.raises(ValueError):
        c.plane_wave_frame(sp.I)
    with pytest.raises(ValueError):
        c.canonical_endpoint_map(sp.Float(".01"), Fraction(1, 10**6))
