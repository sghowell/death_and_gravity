"""Independent physical forward asymptotics and all-angle known-rate checks."""

from functools import cache

import mpmath as mp
import pytest
import sympy as s
from p8_vacuum_affine_finite_forward_phase import (
    audit,
    phase,
    rate,
    regularized,
    source,
)

ROWS, GATES = audit.residuals(), audit.gates()
a, b, ell = s.symbols("a b ell")
master_rows = [source.masters.symbols("m" + str(i)) for i in range(3)]
formal = source.assembly.known_finite(a, b, 1, ell, master_rows)
arguments = (
    a,
    b,
    ell,
    *(row[key] for row in master_rows for key in ("J", "Q", "L", "C", "la")),
)
evaluate = s.lambdify(arguments, formal, "mpmath")
tree = s.lambdify((a, b), source.assembly.tree_jets(a, b, 1)[0], "mpmath")
ENERGIES = ("6.25", "9", "16")
TRANSFERS = ("1e-4", "1e-8", "1e-16")


def lower_log(z):
    z = mp.mpc(z)
    return mp.log(-z.real) - mp.j * mp.pi if z.imag == 0 and z.real < 0 else mp.log(z)


def numeric_triangle(channel):
    def integrand(x):
        z = abs(channel) * x * (1 - x)
        if channel < 0:
            return -(mp.log(z) + mp.pi / mp.sqrt(z)) / (2 * (1 + z))
        h = z - 1
        real = (
            -(1 - h / 2 + h * h / 3)
            if abs(h) < mp.mpf("1e-25")
            else mp.log(z) / (1 - z)
        )
        return -(real + mp.j * mp.pi / (mp.sqrt(z) * (1 + mp.sqrt(z)))) / 2

    return mp.quad(integrand, [0, mp.mpf(".25"), mp.mpf(".5"), mp.mpf(".75"), 1])


def numeric_masters(channel):
    ell = mp.log(4 * mp.pi) - mp.euler
    la = ell - lower_log(-channel)
    if channel < 0:
        J = mp.quad(lambda x: 1 / (1 - channel * x * (1 - x)), [0, mp.mpf(".5"), 1])
        Q = mp.quad(
            lambda x: mp.log(1 - channel * x * (1 - x)) / (1 - channel * x * (1 - x)),
            [0, mp.mpf(".5"), 1],
        )
        L = mp.quad(lambda x: mp.log(1 - channel * x * (1 - x)), [0, mp.mpf(".5"), 1])
    else:

        def state(r):
            x = r + mp.j * r * (1 - r) * (1 - 2 * r) / 4
            xp = 1 + mp.j * (1 - 6 * r + 6 * r * r) / 4
            return xp, 1 - channel * x * (1 - x)

        chart = [0, mp.mpf(".25"), mp.mpf(".5"), mp.mpf(".75"), 1]
        J = mp.quad(lambda r: state(r)[0] / state(r)[1], chart)
        Q = mp.quad(lambda r: state(r)[0] * lower_log(state(r)[1]) / state(r)[1], chart)
        L = mp.quad(lambda r: state(r)[0] * lower_log(state(r)[1]), chart)
    C = numeric_triangle(channel)
    return {
        "J": J,
        "Q": Q,
        "L": L,
        "C": C,
        "la": la,
        "cm": (ell * J - Q) / 2,
        "b00": 2 + la,
        "bmm": ell - L,
    }


@cache
def physical(energy, tau):
    with mp.workdps(110):
        energy, tau = mp.mpf(energy), mp.mpf(tau)
        vals = (energy, -tau, 4 - energy + tau)
        ms = [numeric_masters(x) for x in vals]
        rawlog = mp.log(4 * mp.pi) - mp.euler
        finite = evaluate(
            energy,
            -tau,
            rawlog,
            *(row[key] for row in ms for key in ("J", "Q", "L", "C", "la")),
        )
        v = energy * energy - 4 * energy + 2
        d = mp.sqrt(energy * (energy - 4))
        h = mp.atanh(mp.sqrt((energy - 4) / energy))
        pole = -2 * mp.j * mp.pi * v * v / (d * tau) * (rawlog - mp.log(tau) - 2 / v)
        classical = 3 * mp.pi**2 * (5 * v + 6) / 4
        quantum = (
            4 * v * (energy - 2) * (8 - 3 * v) * h / d**3
            - 2 * v * v / d**2
            + (69 * v + 76) / 30
        )
        # This is the stronger gravity-Born-only denominator, with kappa scaled out.
        relative = 2 * finite.real / (16 * mp.pi**2 * (-tree(energy, -tau)))
        limit = 3 * (5 * v + 6) / (32 * v)
        return tau, finite, pole, classical, quantum, relative, limit


@pytest.mark.parametrize("name", tuple(ROWS))
def test_all_exact_residuals(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_all_written_proof_gates(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_original_scope_and_promotion_rejections(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize("energy", ENERGIES)
@pytest.mark.parametrize("transfer", TRANSFERS)
def test_complete_finite_phase_classical_term_and_known_rate_endpoint(energy, transfer):
    tau, finite, pole, classical, _, relative, limit = physical(energy, transfer)
    with mp.workdps(110):
        assert abs(tau * (finite - pole)) < 2 * classical * mp.sqrt(tau)
        assert abs(mp.sqrt(tau) * finite.real / classical - 1) < 3 * mp.sqrt(tau) * (
            1 + abs(mp.log(tau))
        )
        assert abs(relative / mp.sqrt(tau) - limit) < 3 * mp.sqrt(tau) * (
            1 + abs(mp.log(tau))
        )
        assert abs(finite.real) < mp.mpf(10) ** 14 * (
            1 + abs(mp.log(tau)) + 1 / mp.sqrt(tau)
        )
        assert abs(relative) < mp.mpf(10) ** 12


@pytest.mark.parametrize("energy", ENERGIES)
def test_complete_forward_remainders_decrease(energy):
    rows = [physical(energy, tau) for tau in TRANSFERS]
    with mp.workdps(110):
        classical = [abs(mp.sqrt(row[0]) * row[1].real - row[3]) for row in rows]
        phase_error = [abs(row[0] * (row[1] - row[2])) for row in rows]
        assert classical[2] < classical[1] < classical[0]
        assert phase_error[2] < phase_error[1] < phase_error[0]


@pytest.mark.parametrize("energy", ENERGIES)
def test_real_quantum_log_coefficient_by_two_scale_difference(energy):
    one = physical(energy, "1e-16")
    two = physical(energy, "1e-18")
    with mp.workdps(110):
        got = (
            (one[1].real - one[3] / mp.sqrt(one[0]))
            - (two[1].real - two[3] / mp.sqrt(two[0]))
        ) / mp.log(100)
        assert abs(got - one[4]) < mp.mpf("1e-3")


@pytest.mark.parametrize(
    "energy,transfer",
    (
        ("6.25", ".1"),
        ("6.25", "1"),
        ("6.25", "1.125"),
        ("9", ".1"),
        ("9", "1"),
        ("9", "2.5"),
        ("16", ".1"),
        ("16", "1"),
        ("16", "6"),
    ),
)
def test_whole_near_and_interior_known_real_rate(energy, transfer):
    tau, finite, _, _, _, relative, _ = physical(energy, transfer)
    with mp.workdps(110):
        if tau <= 1:
            assert abs(finite.real) < mp.mpf(10) ** 14 * (
                1 + abs(mp.log(tau)) + 1 / mp.sqrt(tau)
            )
        else:
            assert abs(finite) < mp.mpf(10) ** 11
        assert abs(relative) < mp.mpf(10) ** 12


@pytest.mark.parametrize("tau", ("1e-12", "1e-8", "1e-4", ".1", "1"))
def test_independent_spacelike_triangle_explicit_remainder(tau):
    with mp.workdps(110):
        t = mp.mpf(tau)
        got = numeric_triangle(-t)
        leading = -(mp.pi**2) / (2 * mp.sqrt(t)) - mp.log(t) / 2 + 1
        bound = (
            mp.mpf(5) / 8 * mp.sqrt(t) + t * abs(mp.log(t)) / 12 + mp.mpf(5) / 36 * t
        )
        assert abs(got - leading) < bound


@pytest.mark.parametrize("name", ("J", "Q", "L", "one"))
@pytest.mark.parametrize("order", (0, 1, 2))
def test_exact_rational_derivative_majorant_calibration(name, order):
    expr = s.diff(regularized.coefficients()[name], regularized.TRANSFER, order)
    bound = regularized.majorants()[name][order]
    for energy in (s.Rational(25, 4), s.Integer(16)):
        for transfer in (-s.Integer(1), -s.Rational(1, 2), s.Integer(0)):
            for rawlog in (s.Integer(0), s.Integer(3)):
                value = expr.subs(
                    {
                        regularized.ENERGY: energy,
                        regularized.TRANSFER: transfer,
                        regularized.ELL: rawlog,
                    }
                )
                assert abs(value) <= bound


@pytest.mark.parametrize(
    "denominator",
    (
        regularized.TRANSFER,
        regularized.TRANSFER + 2,
        regularized.ENERGY - 1,
    ),
)
def test_unproved_rational_denominator_rejected(denominator):
    with pytest.raises(ValueError):
        regularized.rational_majorant(1 / denominator)


@pytest.mark.parametrize(
    "energy,transfer",
    (
        (6, -1),
        (17, -1),
        (9, 0),
        (9, -5),
        (9, -6),
        (9, -1.0),
        (9, True),
        (9, None),
        (9, "-1"),
        (9, s.oo),
        (9, s.I),
        (9, s.Symbol("t")),
    ),
)
def test_invalid_energy_angle_domain_rejected(energy, transfer):
    with pytest.raises((TypeError, ValueError)):
        rate.require_physical(energy, transfer)


@pytest.mark.parametrize(
    "energy,transfer",
    (
        (s.Rational(25, 4), -1),
        (16, -6),
        (9, -s.Rational(1, 10**400)),
        (9, -5 + s.Rational(1, 10**400)),
    ),
)
def test_no_artificial_endpoint_window_in_domain(energy, transfer):
    a, b, c, tau = rate.require_physical(energy, transfer)
    assert a + b + c == 4 and tau > 0


@pytest.mark.parametrize("kappa", (1, 10, source.KAPPA))
def test_uniform_rate_bound_exact_scaling(kappa):
    assert rate.all_angle_known_rate_bound(kappa) * kappa == 10**12


def test_full_known_phase_is_not_deleted_and_matching_is_not_assigned():
    assert phase.phase_principal_part().has(s.I, s.log(phase.TAU))
    assert rate.data()["independent_matching_endpoint"].has(
        s.Symbol("finite_delta_kappa")
    )
    assert "three unassigned" in audit.observable()["domain"]
    family = source.assembly.data()["whole_matched_finite_representative"]
    assert family.has(source.assembly.ALPHA, source.assembly.BETA, source.assembly.DK)


def test_original_frontiers_and_historical_qualifications():
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.matching()) == 159 and len(audit.frontier()) == 9
    assert audit.qualifications() == audit.previous.qualifications()
    assert len(audit.qualifications()) == 6
    assert any(row["status"].startswith("REJECTED_") for row in audit.matching())
