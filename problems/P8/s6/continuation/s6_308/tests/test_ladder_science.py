"""Independent raw-rung limits, Fourier integrals, phases and scope checks."""

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_leading_ladder_coulomb_phase import (
    audit,
    bounds,
    dimensional,
    phase,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()


@pytest.mark.parametrize("name", tuple(ROWS))
def test_all_exact_residuals(name):
    value = ROWS[name]
    assert all(
        z == 0 for z in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_original_scope_rejections(name, call, args):
    assert name
    with pytest.raises((ValueError, TypeError)):
        call(*args)


@pytest.mark.parametrize("order", (2, 3, 4))
def test_general_symbolic_ordered_denominator_identity(order):
    variables = s.symbols("x0:" + str(order), positive=True)
    assert dimensional.ordered_identity(variables) == 0


@pytest.mark.parametrize("order", range(1, 7))
def test_raw_full_D_Gamma_rungs_after_exact_same_IR_division(order):
    with mp.workdps(120):
        V, tau = mp.mpf(47), mp.mpf("0.2")
        c = mp.log(tau) - mp.log(4 * mp.pi) + mp.euler + 2 / V

        def factor(strength):
            return (
                mp.exp(mp.j * strength * (c - 2 * mp.euler))
                * mp.gamma(1 - mp.j * strength)
                / mp.gamma(1 + mp.j * strength)
            )

        expected = mp.taylor(factor, 0, order)[order]
        errors = []
        for ep in (mp.mpf("1e-6"), mp.mpf("1e-9")):
            ratio = 1 + 2 * ep / (V * (1 + ep))
            total = 0
            for loops in range(order + 1):
                H = (
                    mp.gamma(1 + ep) ** (loops + 1)
                    * mp.gamma(1 - loops * ep)
                    / mp.gamma(1 + (loops + 1) * ep)
                    * (tau / (4 * mp.pi)) ** (loops * ep)
                )
                total += (
                    (-1) ** (order - loops)
                    * ratio**loops
                    * H
                    / (mp.factorial(loops) * mp.factorial(order - loops))
                )
            got = (mp.j / ep) ** order * total
            errors.append(abs(got - expected))
        assert errors[1] < errors[0] / 100
        assert errors[1] < mp.mpf("1e-5")


@pytest.mark.parametrize(
    "strength,q,h", (("0.1", "1", "0.125"), ("0.4", "2", "0.0625"))
)
def test_literal_radial_Fourier_integral_against_Gaussian_Abel_transform(
    strength, q, h
):
    with mp.workdps(55):
        a, q, h = map(mp.mpf, (strength, q, h))
        upper = mp.sqrt(150 / h)
        count = int(mp.ceil(upper * q / mp.pi))
        knots = [upper * j / count for j in range(count + 1)]
        literal = (
            2
            * mp.pi
            * mp.quad(
                lambda b: (
                    b ** (1 - 2 * mp.j * a) * mp.besselj(0, q * b) * mp.exp(-h * b * b)
                ),
                knots,
            )
        )
        hyper = (
            mp.pi
            * h ** (mp.j * a - 1)
            * mp.gamma(1 - mp.j * a)
            * mp.hyp1f1(1 - mp.j * a, 1, -q * q / (4 * h))
        )
        assert abs(literal - hyper) < mp.mpf("1e-45")


def test_Abel_removal_at_nonzero_transfer_matches_Riesz_distribution():
    with mp.workdps(60):
        a, q = mp.mpf("0.1"), mp.mpf(1)
        riesz = (
            4
            * mp.pi
            / (q * q)
            * mp.gamma(1 - mp.j * a)
            / mp.gamma(mp.j * a)
            * (q * q / 4) ** (mp.j * a)
        )
        errors = []
        for h in (mp.mpf("0.01"), mp.mpf("0.001"), mp.mpf("0.0001")):
            hyper = (
                mp.pi
                * h ** (mp.j * a - 1)
                * mp.gamma(1 - mp.j * a)
                * mp.hyp1f1(1 - mp.j * a, 1, -q * q / (4 * h))
            )
            errors.append(abs(hyper - riesz))
        assert errors[1] < errors[0] / 5 and errors[2] < errors[1] / 5


@pytest.mark.parametrize("strength", ("0.0001", "0.1", "0.4", "0.5"))
@pytest.mark.parametrize("transfer", ("0.2", "1e-30", "1e-1000"))
def test_whole_phase_modulus_and_complete_gamma_error(strength, transfer):
    with mp.workdps(65):
        a, t = map(mp.mpf, (strength, transfer))
        Q = (
            mp.exp(-2 * mp.j * mp.euler * a)
            * mp.gamma(1 - mp.j * a)
            / mp.gamma(1 + mp.j * a)
        )
        c = mp.log(t) - mp.log(4 * mp.pi) + mp.euler + mp.mpf(2) / 47
        whole = mp.exp(mp.j * a * c) * Q
        assert abs(abs(whole) - 1) < mp.mpf("1e-60")
        assert abs(Q - 1) < 2 * a**3
        assert abs(whole - mp.exp(mp.j * a * c)) < 2 * a**3


@pytest.mark.parametrize("excursion", (1, 2, 10, 100))
def test_nonuniform_loop_square_and_complex_monodromy_controls(excursion):
    with mp.workdps(65):
        a = mp.mpf("0.01")
        c = -mp.mpf(excursion) / a
        Q = (
            mp.exp(-2 * mp.j * mp.euler * a)
            * mp.gamma(1 - mp.j * a)
            / mp.gamma(1 + mp.j * a)
        )
        whole = mp.exp(mp.j * a * c) * Q
        truncated = 1 + mp.j * a * c
        assert abs(abs(whole) - 1) < mp.mpf("1e-60")
        assert abs(abs(truncated) ** 2 - 1 - excursion**2) < mp.mpf("1e-55")
        continued = mp.exp(mp.j * a * (c + 2 * mp.j * mp.pi)) * Q
        assert abs(abs(continued) - mp.exp(-2 * mp.pi * a)) < mp.mpf("1e-60")
        assert abs(continued) < 1


@pytest.mark.parametrize("energy", (s.Rational(25, 4), 9, 16))
@pytest.mark.parametrize("tau", (s.Rational(1, 10**400), s.Rational(1, 100), 1))
def test_original_parameter_and_positive_full_Born_embedding(energy, tau):
    assert bounds.original_gamma_remainder_bound(energy, tau) == s.Rational(2, 10**2400)
    assert bounds.original_full_Born_embedding_bound(energy, tau) == s.Rational(
        4, 10**2400
    )
    eta = phase.eta(energy, source.KAPPA)
    assert 0 < eta < 1 / source.KAPPA
    z = 1 - 2 * tau / (energy - 4)
    am = source.born.matter_born(energy, z, source.HEAVY_MASS2, source.CUBIC)
    ag = source.born.gravity_born(energy, z, source.KAPPA)
    pole = source.old_phase.eikonal(energy) / (source.KAPPA * tau)
    assert am > 0 and ag > pole / 2
    assert 0 < pole / (am + ag) < 2


@pytest.mark.parametrize(
    "bad", (True, False, 1.0, s.Float(1), "1", None, s.oo, s.nan, s.I, s.Symbol("x"))
)
def test_nonexact_physical_or_strength_input_rejected(bad):
    with pytest.raises((TypeError, ValueError)):
        bounds.require_domain(9, bad)
    with pytest.raises((TypeError, ValueError)):
        bounds.gamma_remainder_bound(bad)


@pytest.mark.parametrize("energy,tau", ((4, 1), (17, 1), (9, 0), (9, -1), (9, 2)))
def test_outside_stated_physical_domain_rejected(energy, tau):
    with pytest.raises(ValueError):
        bounds.original_gamma_remainder_bound(energy, tau)


@pytest.mark.parametrize("bad", (-1, s.Rational(3, 4), 1))
def test_outside_Gamma_bound_domain_rejected(bad):
    with pytest.raises(ValueError):
        bounds.gamma_remainder_bound(bad)


@pytest.mark.parametrize("bad", (True, False, 1.0, s.Float(1), s.Rational(1, 2), -1))
def test_noninteger_loop_orders_rejected(bad):
    with pytest.raises((ValueError, TypeError)):
        dimensional.normalized_loop(bad)


def test_zero_coupling_and_full_source_frontiers():
    assert phase.gamma_phase(0) == 1 and phase.from_log(0, 7) == 1
    assert bounds.gamma_remainder_bound(0) == 0
    assert dimensional.normalized_loop(0) == 1
    assert len(audit.frontier()) == 9
    assert len(audit.matching()) == 164
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "full interacting amplitude minus" in audit.observable()["not_established"]
    assert "not_full_P8" in audit.MODEL
