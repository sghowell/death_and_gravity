"""Independent endpoint, diagonal, weak-approximation and TT projection tests."""

from itertools import product

import pytest
import sympy as s
from p8_vacuum_affine_continuous_logarithmic_coefficient import (
    audit,
    extension,
    kernels,
    source,
)
from p8_vacuum_affine_radiative_logarithmic_soft_coefficient import radiative
from p8_vacuum_affine_radiative_state_soft_index import recoil

ROWS, GATES = audit.residuals(), audit.gates()
E, STATES = radiative.calibration_states()
U = s.Matrix([0, s.Rational(4, 5), s.Rational(3, 5)])


@pytest.mark.parametrize("name", tuple(ROWS))
def test_every_exact_residual(name):
    value = ROWS[name]
    assert all(
        v == 0 for v in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )
    assert not s.sympify(value).has(s.Float)


@pytest.mark.parametrize("name", tuple(GATES))
def test_every_analytic_proof_gate(name):
    assert GATES[name] is True


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[r[0] for r in audit.bad_cases()]
)
def test_rejected_inputs_and_scope(name, call, args):
    assert name
    with pytest.raises((TypeError, ValueError)):
        call(*args)


@pytest.mark.parametrize(
    "bad",
    (
        None,
        True,
        False,
        "1",
        1.0,
        s.Float(1),
        [],
        {},
        s.Matrix([1]),
        s.I,
        s.oo,
        s.nan,
        s.Symbol("R"),
    ),
)
def test_exact_bound_api_rejects_inexact_or_unresolved_values(bad):
    with pytest.raises((TypeError, ValueError)):
        extension.coefficient_change_upper(bad)


@pytest.mark.parametrize("bad", (-1, s.Rational(1, 7), 1))
def test_radiated_energy_domain(bad):
    with pytest.raises(ValueError):
        extension.coefficient_change_upper(bad)


@pytest.mark.parametrize(
    "R", (0, s.Rational(1, 10**50), s.Rational(1, 128), s.Rational(1, 8))
)
def test_original_parameter_bounds_and_zero_radiation(R):
    assert extension.coefficient_upper() == 300 / s.Integer(10) ** 1200
    assert extension.coefficient_change_upper(R) == 24000 * R / s.Integer(10) ** 1200
    assert source.HEAVY_MASS2 == s.Integer(10) ** 200 / 512 + 2


@pytest.mark.parametrize(
    "name,pol,axis,h",
    tuple(
        product(
            ("one", "two", "four"),
            ("plus", "complex"),
            ("y", "z"),
            (s.Rational(1, 100), s.Rational(1, 10**6), s.Rational(1, 10**12)),
        )
    ),
)
def test_two_azimuthal_approaches_to_unique_complete_coefficient(name, pol, axis, h):
    n0, p0, c0 = kernels.frame_x(0, "y")
    n, p, c = kernels.frame_x(h, axis)
    A0 = p0 if pol == "plus" else (p0 + s.I * c0) / s.sqrt(2)
    A = p if pol == "plus" else (p + s.I * c) / s.sqrt(2)
    target = kernels.coefficient(E, STATES[name], U, n0, A0)
    approached = radiative.coefficient(E, STATES[name], U, n, A)
    assert abs(s.N(approached - target, 90)) < h
    with pytest.raises(ValueError, match="soft-collinear"):
        radiative.coefficient(E, STATES[name], U, n0, A0)


@pytest.mark.parametrize(
    "alpha,pol",
    tuple(
        product(
            (s.Rational(1, 10**20), s.Rational(1, 3), s.Rational(9, 10)),
            ("plus", "complex"),
        )
    ),
)
def test_arbitrary_collinear_atom_splitting_at_its_soft_direction(alpha, pol):
    n, p, c = kernels.frame_x(0, "y")
    A = p if pol == "plus" else (p + s.I * c) / s.sqrt(2)
    base = STATES["two"]
    split = [alpha * base[0], (1 - alpha) * base[0], base[1]]
    assert (
        kernels.clean(
            kernels.coefficient(E, split, U, n, A)
            - kernels.coefficient(E, base, U, n, A)
        )
        == 0
    )


@pytest.mark.parametrize("pol", ("plus", "complex"))
def test_self_diagonal_is_required_by_the_original_coefficient(pol):
    n = s.Matrix([s.Rational(3, 5), s.Rational(4, 5), 0])
    e = s.Matrix([0, -s.Rational(4, 5), s.Rational(3, 5), 0])
    f = s.Matrix([0, 0, 0, 1])
    p = (e * e.T - f * f.T) / s.sqrt(2)
    c = (e * f.T + f * e.T) / s.sqrt(2)
    A = p if pol == "plus" else (p + s.I * c) / s.sqrt(2)
    rays = STATES["one"]
    points, _, _ = recoil.momenta(E, rays, U)
    q = s.Matrix([1, *n])
    P = lambda a, b: (a.T * A * b)[0]
    D = [recoil.dot(k, q) for k in points]
    L = [P(k, k) / d for k, d in zip(points, D)]
    from p8_vacuum_affine_radiative_logarithmic_soft_coefficient import kernel

    without = sum(L) * sum(d * s.log(abs(d)) for d in D)
    for i in range(4):
        for j in range(i + 1, 4):
            S = L[i] * D[j] + L[j] * D[i] - 2 * P(points[i], points[j])
            without -= (
                kernel.massive_fprime(abs(recoil.dot(points[i], points[j]))) / 4
                - s.log(2)
                - 1
            ) * S
    w = rays[0][0]
    v = rays[0] / w
    d = recoil.dot(v, q)
    ell = P(v, v) / d
    for i, k in enumerate(points):
        alpha = abs(recoil.dot(k, v))
        without += w * (
            L[i] * d * s.log(d / alpha)
            + D[i] * ell * s.log(abs(D[i]) / alpha)
            + 2 * P(k, v) * s.log(alpha)
        )
    original = radiative.components(E, rays, U, n, A)["F"]
    diagonal = w * w * P(v, v) * s.log(d) / (4 * s.pi**2)
    assert diagonal != 0
    assert kernels.clean(original - without / (4 * s.pi**2) - diagonal) == 0


def axis_frame(axis):
    if axis == 0:
        n = s.Matrix([1, 0, 0])
        e = s.Matrix([0, 0, 1, 0])
        f = s.Matrix([0, 0, 0, 1])
    elif axis == 1:
        n = s.Matrix([0, 1, 0])
        e = s.Matrix([0, 1, 0, 0])
        f = s.Matrix([0, 0, 0, 1])
    else:
        n = s.Matrix([0, 0, 1])
        e = s.Matrix([0, 1, 0, 0])
        f = s.Matrix([0, 0, 1, 0])
    return n, (e * e.T - f * f.T) / s.sqrt(2), (e * f.T + f * e.T) / s.sqrt(2)


@pytest.mark.parametrize(
    "axis,pol,h",
    tuple(
        product(
            (0, 1, 2), ("plus", "complex"), (s.Rational(1, 100), s.Rational(1, 10000))
        )
    ),
)
def test_moving_atomic_energy_approximation_including_the_limit_atom(axis, pol, h):
    n, p, c = axis_frame(axis)
    A = p if pol == "plus" else (p + s.I * c) / s.sqrt(2)
    w = STATES["one"][0][0]
    a, b = (1 - h * h) / (1 + h * h), 2 * h / (1 + h * h)
    approximants = [w * s.Matrix([1, a, b, 0]) / 2, w * s.Matrix([1, a, -b, 0]) / 2]
    target = kernels.coefficient(E, STATES["one"], U, n, A)
    approached = kernels.coefficient(E, approximants, U, n, A)
    assert abs(s.N(approached - target, 90)) < h
    assert sum(v[0] for v in approximants) == w
    assert sum(approximants, s.zeros(4, 1)) != STATES["one"][0]


@pytest.mark.parametrize("axis", (0, 1, 2))
def test_global_TT_projector_without_a_global_polarization_choice(axis):
    n, _, _ = axis_frame(axis)
    T = s.Matrix([[1, 2 + s.I, 3], [2 + s.I, -2, 5], [3, 5, 4]])
    P = extension.tt_project(T, n)
    assert P == P.T and P * n == s.zeros(3, 1) and s.trace(P) == 0
    assert extension.tt_project(P, n) == P
    assert s.trace(P.conjugate().T * P) <= s.trace(T.conjugate().T * T)


@pytest.mark.parametrize(
    "mutation",
    (
        "energy_none",
        "wrong_u",
        "wrong_n",
        "missing_quanta",
        "zero_ray",
        "bad_ray",
        "bad_tensor",
        "traceful",
        "non_TT",
        "nonunit",
    ),
)
def test_complete_coefficient_api_preserves_physical_input_validation(mutation):
    n, p, _ = kernels.frame_x(0, "y")
    args = [E, STATES["one"], U, n, p]
    if mutation == "energy_none":
        args[0] = None
    elif mutation == "wrong_u":
        args[2] = [1, 0]
    elif mutation == "wrong_n":
        args[3] = [0, 0, 2]
    elif mutation == "missing_quanta":
        args[1] = None
    elif mutation == "zero_ray":
        args[1] = [s.zeros(4, 1)]
    elif mutation == "bad_ray":
        args[1] = [s.Matrix([s.Rational(1, 32), 0, 0, 0])]
    elif mutation == "bad_tensor":
        args[4] = None
    elif mutation == "traceful":
        args[4] = s.diag(0, 0, 1, 1) / s.sqrt(2)
    elif mutation == "non_TT":
        args[4] = s.diag(0, 1, -1, 0) / s.sqrt(2)
    else:
        args[4] = 2 * p
    with pytest.raises((TypeError, ValueError)):
        kernels.coefficient(*args)


def test_scope_preserves_original_frontier_and_individual_current_qualification():
    assert len(audit.frontier()) == 9 and len(audit.matching()) == 184
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len(audit.qualifications()) == 6
    assert "No unique collinear value" in audit.observable()["not_established"]
    assert "interacting quantum state" in audit.observable()["not_established"]
