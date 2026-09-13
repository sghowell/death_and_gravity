"""Full exact identities and independent constrained classical diagnostics."""

from functools import cache

import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_classical_principal_realization import (
    audit,
    homogeneous,
    realization,
    rolling,
)
from p8_vacuum_affine_coupled_principal_obstruction import background as old
from p8_vacuum_affine_coupled_principal_obstruction import coupled as c
from p8_vacuum_affine_coupled_principal_obstruction import principal
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_complete_exact_identity(name):
    value = audit.residuals()[name]
    assert all(
        entry == 0
        for entry in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_every_domain_and_scope_rejection(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_counts_scope_and_full_domain():
    assert (
        len(audit.residuals()),
        audit.scalar_entry_count(),
        len(audit.gates()),
        len(audit.controls()),
        audit.rejected_inputs(),
        len(audit.matching()),
    ) == (56, 104, 27, 9, 317, 111)
    assert all(value is True for value in audit.gates().values())
    assert audit.require_domain(audit.example_domain()) == audit.example_domain()
    assert audit.STATE != audit.previous.STATE


@pytest.mark.parametrize("i", (-1, 0, 1))
@pytest.mark.parametrize("j", (-1, 0, 1))
def test_independent_exact_current_profile_root_and_margins(i, j):
    rho, pressure = i * realization.PROFILE_BOUND, j * realization.PROFILE_BOUND
    rate = audit.central_rate(rho, pressure)
    T = rho - 3 * pressure / 2
    assert rate.is_positive is True
    assert s.factor(-s.Rational(1, 400) + T + rate**2 / 4) == 0
    assert s.factor(rate**2 - s.Rational(1, 100) + 4 * T) == 0
    J = s.Rational(243, 160) + rho - 7 * pressure / 8 + 3 * T / 4
    Cnn = s.Rational(38, 25) + rho - 7 * pressure / 8 + T / 4
    assert J > s.Rational(3, 2) and Cnn > s.Rational(3, 2)
    # Rationalize rather than numerically subtract two roots separated by
    # less than10^-397. This retains the actual finite profile correction.
    assert abs(-4 * T) / (rate + s.Rational(1, 10)) < s.Integer(1) / 10**397


@pytest.mark.parametrize("case", range(3))
def test_independent_entire_homogeneous_lapse_legendre_block(case):
    D, Z, J = s.Rational(5, 4), s.Rational(6, 5), s.Rational(7, 4)
    theta = s.Rational(case - 1, 10)
    w1, w2 = s.Rational(1, 7), -s.Rational(2, 9)
    Cnn = J - 3 * theta**2 / D + (w1**2 + w2**2) / (2 * Z)
    v = s.Matrix(s.symbols("independent_velocity0:3"))
    n = s.Symbol("independent_lapse")
    cross = s.Matrix([6 * theta, w1, w2])
    raw = (
        -3 * D * v[0] ** 2
        + Z * (v[1] ** 2 + v[2] ** 2) / 2
        + n * cross.dot(v)
        + Cnn * n**2
    )
    solution = s.solve(s.diff(raw, n), n)[0]
    reduced = s.hessian(raw.subs(n, solution), v)
    full = s.hessian(raw, [*v, n])
    assert (full.inv()[:3, :3] * reduced - s.eye(3)).applyfunc(s.factor) == s.zeros(3)
    assert s.factor(reduced.det() + 6 * D * Z**2 * J / Cnn) == 0
    assert reduced[1, 2] != 0


def toy_packet(time, field, velocity, lapse):
    """A finite-parameter diagnostic, NOT the giant-mass physical solution.

    R=N^-2, D=N^-3/2,Z=N^1/2,B0; the homogeneous geometric potential
    is -8+t/100-(N-1)/400+2(N-1)^2. The minimally coupled heavy potential
    is N^(5/2)[-9h^2/2+j(t,N)h], j=1/20+t/30+(N-1)/40.
    """
    H, mc, mh = velocity
    h = field[2]
    N = lapse
    D, DN, DNN = N**-1.5, -1.5 * N**-2.5, 3.75 * N**-3.5
    Z, ZN, ZNN = N**0.5, 0.5 * N**-0.5, -0.25 * N**-1.5
    W, WN, WNN = N**2.5, 2.5 * N**1.5, 3.75 * N**0.5
    source, sourceN, sourceT = 0.05 + time / 30 + (N - 1) / 40, 0.025, 1 / 30
    potential = -4.5 * h**2 + source * h
    geometric = -8 + time / 100 - (N - 1) / 400 + 2 * (N - 1) ** 2
    density = -3 * D * H**2 + Z * (mc**2 + mh**2) / 2 + geometric + W * potential
    constraint = (
        -3 * DN * H**2
        + ZN * (mc**2 + mh**2) / 2
        - 1 / 400
        + 4 * (N - 1)
        + WN * potential
        + W * sourceN * h
    )
    Cnn = (
        -3 * DNN * H**2
        + ZNN * (mc**2 + mh**2) / 2
        + 4
        + WNN * potential
        + 2 * WN * sourceN * h
    ) / 2
    cross = np.array([-6 * DN * H, ZN * mc, ZN * mh])
    K0 = np.diag([-6 * D, Z, Z])
    momentum = np.array([-6 * D * H, Z * mc, Z * mh])
    force = np.array([3 * density, 0, W * (-9 * h + source)])
    Qt = WN * sourceT * h
    Qh = WN * (-9 * h + source) + W * sourceN
    return density, constraint, Cnn, cross, K0, momentum, force, Qt + Qh * mh


def toy_lapse(time, field, velocity):
    return brentq(
        lambda lapse: toy_packet(time, field, velocity, lapse)[1], 0.8, 1.2, xtol=2e-14
    )


@pytest.mark.parametrize("case", range(3))
def test_independent_nonlinear_implicit_lapse_full_Hessian(case):
    time = 0.01 * case
    field = np.array([0.01, 0.02, 0.005])
    velocity = np.array([0.05 + 0.01 * case, 0.1, 0.012])
    lapse = toy_lapse(time, field, velocity)
    _, _, Cnn, cross, K0, _, _, _ = toy_packet(time, field, velocity, lapse)
    expected = K0 - np.outer(cross, cross) / (2 * Cnn)

    def reduced(v):
        return toy_packet(time, field, v, toy_lapse(time, field, v))[0]

    step = 2e-4
    actual = np.zeros((3, 3))
    base = reduced(velocity)
    for i in range(3):
        ei = step * np.eye(3)[i]
        actual[i, i] = (
            reduced(velocity + ei) - 2 * base + reduced(velocity - ei)
        ) / step**2
        for j in range(i):
            ej = step * np.eye(3)[j]
            actual[i, j] = actual[j, i] = (
                reduced(velocity + ei + ej)
                - reduced(velocity + ei - ej)
                - reduced(velocity - ei + ej)
                + reduced(velocity - ei - ej)
            ) / (4 * step**2)
    assert np.linalg.norm(actual - expected) < 2e-5
    assert np.linalg.norm(actual - K0) > 1e-3


@cache
def full_evolution_diagnostic():
    initial = np.array([0, 0, 0, 0, 0.1, 0], dtype=float)
    N0 = toy_lapse(0, initial[:3], initial[3:])

    def reduced_rhs(time, state):
        field, velocity = state[:3], state[3:]
        lapse = toy_lapse(time, field, velocity)
        _, _, Cnn, cross, K0, momentum, force, rest = toy_packet(
            time, field, velocity, lapse
        )
        kinetic = K0 - np.outer(cross, cross) / (2 * Cnn)
        acceleration = np.linalg.solve(
            kinetic, force - 3 * velocity[0] * momentum + cross * rest / (2 * Cnn)
        )
        return np.r_[velocity, acceleration]

    def full_rhs(time, state):
        field, velocity, lapse = state[:3], state[3:6], state[6]
        _, _, Cnn, cross, K0, momentum, force, rest = toy_packet(
            time, field, velocity, lapse
        )
        joint = np.block(
            [[K0, cross[:, None]], [cross[None, :], np.array([[2 * Cnn]])]]
        )
        derivative = np.linalg.solve(
            joint, np.r_[force - 3 * velocity[0] * momentum, -rest]
        )
        return np.r_[velocity, derivative]

    times = np.linspace(0, 0.01, 17)
    reduced = solve_ivp(
        reduced_rhs,
        (0, 0.01),
        initial,
        t_eval=times,
        method="DOP853",
        rtol=2e-11,
        atol=2e-13,
    )
    full = solve_ivp(
        full_rhs,
        (0, 0.01),
        np.r_[initial, N0],
        t_eval=times,
        method="DOP853",
        rtol=2e-11,
        atol=2e-13,
    )
    assert reduced.success and full.success
    return times, reduced.y, full.y


def test_independent_full_and_reduced_unforced_evolution_agree():
    times, reduced, full = full_evolution_diagnostic()
    assert np.max(np.abs(reduced - full[:6])) < 3e-10
    assert np.max(np.abs(full[2])) > 1e-7
    assert np.max(np.abs(full[5])) > 1e-4
    for slot, time in enumerate(times):
        field, velocity, lapse = full[:3, slot], full[3:6, slot], full[6, slot]
        _, constraint, Cnn, _, _, _, _, _ = toy_packet(time, field, velocity, lapse)
        assert abs(constraint) < 2e-10 and Cnn > 1.8
        assert abs(toy_lapse(time, field, velocity) - lapse) < 2e-10


def test_independent_full_evolution_retains_M1_charge():
    _, _, full = full_evolution_diagnostic()
    charge = np.exp(3 * full[0]) * np.sqrt(full[6]) * full[4]
    assert np.max(np.abs(charge - charge[0])) < 2e-10


@pytest.mark.parametrize("case", range(3))
def test_independent_full_nonlinear_field_map_onepoint_contact(case):
    x, y, e = s.symbols("x y e")
    point = s.Rational(case, 5)
    original = (x - point) ** 2 / 2 + (case + 1) * (x - point) ** 3 / 3
    mapping = point + y + e * y**2
    assert s.diff(original.subs(x, mapping), y, 2).subs(y, 0) == 1
    offshell = original + x
    assert s.diff(offshell.subs(x, mapping), y, 2).subs(y, 0) - 1 == 2 * e


def test_holding_heavy_zero_is_not_unforced_initial_acceleration():
    packet = rolling.coefficient_map()
    heavy_euler = packet["entire_homogeneous_matter_Euler_expressions"][1]
    held = {homogeneous.hbar: 0, homogeneous.mh: 0, rolling.mhd: 0}
    residual = heavy_euler.subs(held, simultaneous=True)
    assert residual != 0 and residual.has(old.JH)
    acceleration = s.solve(
        heavy_euler.subs({homogeneous.hbar: 0, homogeneous.mh: 0}, simultaneous=True),
        rolling.mhd,
    )[0]
    assert s.factor(acceleration - homogeneous.N**2 * old.JH / s.sqrt(old.kappa)) == 0


def test_full_rolling_map_has_no_old_scalar_rate_placeholder():
    for value in rolling.coefficient_map()[
        "entire_rolling_four_mode_parameter_map"
    ].values():
        assert not s.sympify(value).has(old.m, old.md)


def test_on_shell_lower_contacts_do_not_remove_fast_branch():
    fast = principal.fast_data()["whole_eight_phase_fast_symbol"]
    source = {c.L0: 0, c.Vvv: 0, c.vs[0]: 0, c.vs[1]: 0}
    assert fast.subs(source, simultaneous=True) == fast
    assert fast[0, 4] != 0 and fast[3, 0] != 0 and fast[4, 7] != 0 and fast[7, 3] != 0
