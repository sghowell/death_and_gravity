"""Exact full-block checks and independently re-entered finite diagnostics."""

from functools import cache

import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_coupled_principal_obstruction import audit, background, principal
from p8_vacuum_affine_coupled_principal_obstruction import coupled as c
from p8_vacuum_affine_physical_background_vertices import coupled as prior
from scipy.integrate import solve_ivp


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_full_exact_identity(name):
    value = audit.residuals()[name]
    assert all(
        entry == 0
        for entry in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_every_invalid_domain_or_scope_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_all_counts_and_complete_domain():
    assert (
        len(audit.residuals()),
        audit.scalar_entry_count(),
        len(audit.gates()),
        len(audit.controls()),
        audit.rejected_inputs(),
        len(audit.matching()),
    ) == (51, 813, 32, 9, 276, 110)
    assert all(value is True for value in audit.gates().values())
    assert audit.require_domain(audit.example_domain()) == audit.example_domain()


def example(index=0):
    values = {
        c.D: s.Rational(5, 4),
        c.Z: s.Rational(4, 5),
        c.J: s.Rational(3, 2),
        c.K: s.Rational(7, 3),
        c.a: s.Rational(6, 5),
        c.zeta: s.Rational(1, 7),
        c.P: s.Integer(3 + index),
        c.Y: s.Rational(8, 9),
        c.Yv: s.Rational(5, 8),
        c.r: s.Rational(index - 1, 5),
        c.rN: s.Rational(2, 7),
        c.th: s.Rational(3, 8),
        c.dH: s.Rational(1, 9),
        c.H: s.Rational(1, 10),
        c.L0: s.Rational(2, 9),
        c.L2: -s.Rational(7, 8),
        c.C: s.Rational(3, 5),
        c.Vvv: -s.Rational(2, 11),
        c.cs[0]: s.Rational(1, 10),
        c.cs[1]: s.Rational(1, 11),
        c.ws[0]: s.Rational(1, 12),
        c.ws[1]: -s.Rational(1, 13),
        c.ds[0]: s.Rational(2, 13),
        c.ds[1]: -s.Rational(2, 15),
        c.vs[0]: s.Rational(1, 17),
        c.vs[1]: -s.Rational(2, 19),
        c.m00: -s.Rational(1, 2),
        c.m01: s.Rational(2, 21),
        c.m11: -s.Integer(9),
    }
    for name in ("theta_dot", "D_dot", "L0_dot", "L2_dot"):
        values[s.Symbol(name, real=True)] = s.Rational(1, 10 + index)
    return values


@pytest.mark.parametrize("index", range(3))
def test_independent_joint_legendre_and_auxiliary_reduction(index):
    lag = s.expand(c.lagrangian().subs(example(index), simultaneous=True))
    # Solve the FOUR momentum definitions and THREE auxiliary Euler equations
    # simultaneously from the raw Lagrangian, independently of stored solutions.
    equations = [
        s.diff(lag, rate) - momentum for rate, momentum in zip(c.RATES, c.MOMENTA)
    ]
    equations += [s.diff(lag, x) for x in (c.n, c.b, c.A0)]
    unknowns = [*c.RATES, c.n, c.b, c.A0]
    matrix, rhs = s.linear_eq_to_matrix(equations, unknowns)
    solution = dict(zip(unknowns, matrix.inv() * rhs))
    actual = s.expand((c.MOMENTA.dot(c.RATES) - lag).subs(solution, simultaneous=True))
    expected = c.hamiltonian().subs(example(index), simultaneous=True)
    assert s.factor(actual - expected) == 0


def test_complete_previous_three_mode_block_is_exact_subrestriction():
    disappear = {
        c.sigma[1]: 0,
        c.ps[1]: 0,
        c.cs[1]: 0,
        c.ws[1]: 0,
        c.ds[0]: 0,
        c.ds[1]: 0,
        c.m00: 0,
        c.m01: 0,
        c.m11: 0,
    }
    old_map = {
        prior.D: c.D,
        prior.Z: c.Z,
        prior.J: c.J,
        prior.K: c.K,
        prior.P: c.P,
        prior.th: c.th,
        prior.w: c.ws[0],
        prior.c: c.cs[0],
        prior.Lnv: c.Lnv,
        prior.Vvv: c.Vvv,
        prior.Vvs: c.vs[0],
        prior.C: c.C,
        prior.Y: c.Y,
        prior.Yv: c.Yv,
        prior.r: c.r,
        prior.rN: c.rN,
        prior.dH: c.dH,
        prior.v: c.v,
        prior.sigma: c.sigma[0],
        prior.longitudinal: c.A,
        prior.pv: c.pv,
        prior.ps: c.ps[0],
        prior.pA: c.pA,
        prior.q: c.q,
    }
    assert (
        s.factor(
            c.hamiltonian().subs(disappear, simultaneous=True)
            - prior.reduced_hamiltonian().subs(old_map, simultaneous=True)
        )
        == 0
    )


@pytest.mark.parametrize("index", range(3))
def test_independent_polynomial_heavy_source_density(index):
    n, v, h, hd, e = s.symbols("n v h hd e")
    lapse = s.Rational(7 + index, 6)
    x = lapse + e * n
    U, source, Cchi = 1 + x + x**3, 2 - x + x**4, 1 + x**2
    volume = 1 + 3 * e * v + s.Rational(9, 2) * e**2 * v**2
    density = (
        U * (e * hd) ** 2 / (2 * x)
        - x * U * 9 * (e * h) ** 2 / 2
        + x * U * source * e * h / 5
    ) * volume
    density -= x * Cchi * (1 + e * v) * 4 * (e * h) ** 2 / 2
    actual = s.diff(density, e, 2).subs(e, 0) / 2
    y = s.Symbol("y")
    source_density = y * (1 + y + y**3) * (2 - y + y**4) / 5
    expected = (
        (1 + lapse + lapse**3) * hd**2 / (2 * lapse)
        - lapse * (1 + lapse + lapse**3) * 9 * h**2 / 2
        - lapse * (1 + lapse**2) * 4 * h**2 / 2
        + (
            source_density.diff(y).subs(y, lapse) * n
            + 3 * source_density.subs(y, lapse) * v
        )
        * h
    )
    assert s.factor(actual - expected) == 0
    assert s.diff(actual, n, h) != 0 and s.diff(actual, v, h) != 0


@pytest.mark.parametrize("index", range(3))
def test_whole_exact_time_connection_chain_independent_path(index):
    time = s.Symbol("diagnostic_time", real=True)
    values = example(index)
    path = {
        c.a: values[c.a] * s.exp(time / 10),
        c.H: s.Rational(1, 10),
        c.D: values[c.D] + time / 7,
        c.th: values[c.th] + time / 9,
        c.L0: values[c.L0] + time / 11,
        c.L2: values[c.L2] + time / 13,
    }
    R = principal.fast_data()["whole_exact_old_from_new_shear"]
    derivative = (
        R.subs(path, simultaneous=True)
        .diff(time)
        .subs(time, 0)
        .subs(values, simultaneous=True)
    )
    expected_values = {
        **values,
        s.Symbol("D_dot", real=True): s.Rational(1, 7),
        s.Symbol("theta_dot", real=True): s.Rational(1, 9),
        s.Symbol("L0_dot", real=True): s.Rational(1, 11),
        s.Symbol("L2_dot", real=True): s.Rational(1, 13),
    }
    expected = principal.fast_data()["whole_shear_time_derivative"].subs(
        expected_values, simultaneous=True
    )
    assert (derivative - expected).applyfunc(s.factor) == s.zeros(8)
    assert derivative != s.zeros(8)


@pytest.mark.parametrize("index", range(3))
def test_full_finite_momentum_remainder_not_only_determinant(index):
    packet = principal.fast_data()
    values = example(index)
    del values[c.P]
    normalized = packet["whole_normalized_generator"].subs(values, simultaneous=True)
    leading = np.array(
        packet["whole_eight_phase_fast_symbol"].subs(values, simultaneous=True),
        dtype=float,
    )
    evaluate = s.lambdify(c.P, normalized, modules="numpy", cse=True)
    previous = None
    for P in (256, 1024, 4096):
        error = np.linalg.norm(np.asarray(evaluate(P), dtype=float) - leading)
        assert np.isfinite(error)
        if previous is not None:
            assert error < 0.65 * previous
        previous = error


@cache
def dynamic_diagnostic():
    t = s.Symbol("time", real=True)
    values = example(2)
    del values[c.P]
    values.update(
        {
            c.a: 1,
            c.H: 0,
            c.D: 1,
            c.Z: 1,
            c.J: 1,
            c.zeta: 1,
            c.Y: 1,
            c.Yv: 1,
            c.K: 1,
            c.L2: -1,
            c.L0: 0,
            c.r: s.Rational(1, 8),
            c.th: t / 10,
            c.dH: t / 11,
        }
    )
    values.update(
        {
            s.Symbol("theta_dot", real=True): s.Rational(1, 10),
            s.Symbol("D_dot", real=True): 0,
            s.Symbol("L0_dot", real=True): 0,
            s.Symbol("L2_dot", real=True): 0,
        }
    )
    whole = principal.fast_data()["whole_normalized_generator"].subs(
        values, simultaneous=True
    )
    fast = np.array(
        principal.fast_data()["whole_eight_phase_fast_symbol"].subs(
            values, simultaneous=True
        ),
        dtype=float,
    )
    eigenvalues, eigenvectors = np.linalg.eig(fast)
    slot = int(np.argmax(eigenvalues.real))
    rho = float(eigenvalues[slot].real)
    vector = eigenvectors[:, slot].real
    vector /= np.linalg.norm(vector)
    evaluate = s.lambdify((t, c.P), whole, modules="numpy", cse=True)
    logs = []
    for P in (64, 256, 1024):

        def rhs(time, state, momentum=P):
            return (
                momentum**1.5
                * np.asarray(evaluate(time, momentum), dtype=float)
                @ state
            )

        result = solve_ivp(
            rhs, (0, 0.02), vector, method="DOP853", rtol=2e-10, atol=2e-12
        )
        assert result.success
        rate = np.log(np.linalg.norm(result.y[:, -1])) / (P**1.5 * 0.02)
        logs.append(float(rate))
    return rho, logs


def test_full_time_dependent_evolution_has_fast_growth_diagnostic():
    rho, rates = dynamic_diagnostic()
    assert rates[-1] > rho / 2
    assert abs(rates[-1] - rho) < 0.08 * rho
    assert abs(rates[-1] - rho) < abs(rates[0] - rho)


@pytest.mark.parametrize("index", range(3))
def test_entire_finite_matrix_fast_spectrum_converges(index):
    values = example(index)
    values[c.r] = s.Rational(1 + index, 10)
    del values[c.P]
    packet = principal.fast_data()
    fast_rate = (
        float(packet["positive_fast_rate_fourth_power"].subs(values, simultaneous=True))
        ** 0.25
    )
    normalized = packet["whole_normalized_generator"].subs(values, simultaneous=True)
    evaluate = s.lambdify(c.P, normalized, modules="numpy", cse=True)
    eigenvalues = np.linalg.eigvals(np.asarray(evaluate(65536), dtype=float))
    for target in (fast_rate, -fast_rate, 1j * fast_rate, -1j * fast_rate):
        assert np.min(np.abs(eigenvalues - target)) < 0.015 * fast_rate


def test_negative_source_omission_controls():
    full = c.lagrangian()
    assert s.diff(full, c.n, c.sigma[1]).has(c.ds[1])
    assert s.diff(full, c.v, c.sigma[1]).has(c.vs[1])
    assert s.diff(c.hamiltonian(), c.ds[1]) != 0
    assert s.diff(c.hamiltonian(), c.vs[1]) != 0
    assert principal.fast_data()["positive_fast_rate_fourth_power"].subs(c.r, 0) == 0
    assert principal.fast_data()["whole_shear_time_connection"] != s.zeros(8)


def test_current_physical_Hbar0_binding_retains_full_source():
    parameters = background.coefficients()["whole_parameter_substitution"]
    assert parameters[c.cs[1]] == parameters[c.ws[1]] == 0
    assert parameters[c.ds[1]] != 0 and parameters[c.vs[1]] != 0
    assert parameters[c.m11] != 0
