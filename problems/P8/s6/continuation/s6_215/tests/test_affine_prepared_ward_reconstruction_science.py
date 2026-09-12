"""Independent covariance, action, endpoint, projection and norm checks."""

from functools import cache

import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_full_adm_vertices import chart
from p8_vacuum_affine_prepared_ward_reconstruction import kinematics as k
from p8_vacuum_affine_prepared_ward_reconstruction import norms, scalar, ward

t, x, y, z = k.COORDS
Z3 = s.zeros(3)
ZERO = s.zeros(3, 1)
g0 = s.diag(1, -(k.a**2), -(k.a**2), -(k.a**2))
gi = s.diag(1, -1 / k.a**2, -1 / k.a**2, -1 / k.a**2)
EV = k.a**3 * gi / 2


def zero(expr):
    entries = list(expr) if isinstance(expr, s.MatrixBase) else [expr]
    assert all(s.cancel(v) == 0 for v in entries)


def volume_tangent(h):
    return k.a**3 * (s.trace(gi * h) * gi / 2 - gi * h * gi) / 2


def volume_response(D, G):
    nD, _bD, QD = D
    nG, _bG, QG = G
    return k.a**3 * (
        s.trace(QD) * s.trace(QG) / 4 + (nD * s.trace(QG) + nG * s.trace(QD)) / 2
    )


def integral(expr):
    return s.integrate(
        s.integrate(s.cancel(expr), (x, -1, 1)),
        (t, -s.Rational(1, 2), s.Rational(1, 2)),
    )


DIRECTIONS = (
    (1 + t * x, ZERO, Z3),
    (0, s.Matrix([t + x, x * x, 1 + t]), Z3),
    (0, ZERO, (1 + t * x) * s.eye(3)),
    (0, ZERO, s.Matrix([[t, x, 1], [x, -t, t * x], [1, t * x, 0]])),
    (
        t + x,
        s.Matrix([x, t * x, t]),
        s.Matrix([[1, x, t], [x, 2, x * t], [t, x * t, 3]]),
    ),
)
XIS = (
    s.Matrix([k.a**2 * (t + s.Rational(1, 2)) ** 2 * x, 0, 0, 0]),
    s.Matrix([0, (1 + t) * x * x, (1 + t * t) * x, t * x]),
    s.Matrix([k.a**2 * (1 + t) * x * x, t * x, t * t * x, (1 + t) * x * x]),
)


@pytest.mark.parametrize("side", ("retarded", "advanced"))
@pytest.mark.parametrize("index", range(4))
def test_exact_source_detector_inverse_maps(side, index):
    end = -s.Rational(1, 2) if side == "retarded" else s.Rational(1, 2)
    f = (t - end) ** (index + 2)
    eta = k.a**2 * f * (x + index * x * x)
    chi = f * s.Matrix([x * x, (index + 1) * x, x**3])
    n, b, Qxi = k.gauge(eta, chi)
    Q = f * s.Matrix([[1, x, index], [x, 2, t], [index, t, 3]])
    syn, er, cr = k.synchronous(n, b, Q, side)
    zero(er - eta)
    zero(cr - chi)
    zero(syn - (Q - Qxi))
    zero(
        chart.first(n, b, Q, k.a)
        - chart.first(0, ZERO, syn, k.a)
        - k.lie_covariant(g0, s.Matrix([eta, *chi]))
    )


@pytest.mark.parametrize("xi", XIS)
def test_density_source_covariance_from_literal_inverse_and_determinant(xi):
    h = k.lie_covariant(g0, xi)
    zero(volume_tangent(h) - ward.lie_density(EV, xi))


@pytest.mark.parametrize("D", DIRECTIONS)
@pytest.mark.parametrize("xi", XIS)
def test_source_Ward_with_full_volume_action_chart(D, xi):
    G = k.gauge(xi[0], xi[1:, 0])
    zero(volume_response(D, G) - ward.source_integrand(D, xi, EV))


@pytest.mark.parametrize("G", DIRECTIONS)
@pytest.mark.parametrize("xi", XIS)
def test_detector_Ward_includes_exact_differentiated_flux(G, xi):
    h = chart.first(*G, k.a)
    dE = volume_tangent(h)
    D = k.gauge(xi[0], xi[1:, 0])
    flux = [
        2
        * sum(
            xi[alpha]
            * (
                sum(
                    dE[mu, nu] * g0[alpha, nu] + EV[mu, nu] * h[alpha, nu]
                    for nu in range(4)
                )
            )
            for alpha in range(4)
        )
        for mu in range(4)
    ]
    divergence = sum(s.diff(flux[j], q) for j, q in enumerate(k.COORDS))
    zero(volume_response(D, G) - ward.detector_integrand(xi, G, EV) - divergence)


@pytest.mark.parametrize("index", range(3))
def test_integrated_full_ordered_volume_reconstruction(index):
    def fixture(end, shift):
        f = (t - end) ** 3
        space = (1 - x * x) ** 2
        eta = k.a**2 * f * space * (1 + shift * x)
        chi = f * space * s.Matrix([1 + shift * x, x, x * x])
        n, beta, _ = k.gauge(eta, chi)
        Q = f * s.Matrix([[1 + x, x, shift], [x, 2 - x, t], [shift, t, 1 + t * x]])
        return (n, beta, Q)

    D = fixture(s.Rational(1, 2), index + 1)
    G = fixture(-s.Rational(1, 2), index + 2)
    ds, de, dc = k.synchronous(*D, "advanced")
    gs, ge, gc = k.synchronous(*G, "retarded")
    reconstructed = volume_response((0, ZERO, ds), (0, ZERO, gs))
    reconstructed += ward.source_integrand(D, s.Matrix([ge, *gc]), EV)
    reconstructed += ward.detector_integrand(s.Matrix([de, *dc]), (0, ZERO, gs), EV)
    zero(integral(volume_response(D, G) - reconstructed))


@pytest.mark.parametrize("side", ("retarded", "advanced"))
@pytest.mark.parametrize("degree", range(6))
def test_time_primitive_L2_bound_direct_integral(side, degree):
    f = (t + s.Rational(1, 2)) ** degree
    F = k.primitive(f, side)
    fn = s.integrate(f * f, (t, -s.Rational(1, 2), s.Rational(1, 2)))
    Fn = s.integrate(F * F, (t, -s.Rational(1, 2), s.Rational(1, 2)))
    assert 0 <= Fn <= fn
    zero(s.diff(F, t) - f)


@cache
def coefficient_jet(kind, order):
    f = k.H if kind == "H" else 1 / k.a**2
    return s.lambdify(t, s.diff(f, t, order), "numpy")


@pytest.mark.parametrize("kind,leading", (("H", 7), ("c", 28)))
@pytest.mark.parametrize("order", range(14))
def test_literal_coefficient_jets_inside_Cauchy_majorants(kind, leading, order):
    actual = np.asarray(coefficient_jet(kind, order)(np.linspace(-0.5, 0.5, 25)))
    bound = float(leading * s.factorial(order) * 4**order)
    assert np.max(np.abs(actual)) < bound


AXES = (
    s.Matrix([0, 0, 1]),
    s.Matrix([s.Rational(3, 5), s.Rational(4, 5), 0]),
    s.Matrix([s.Rational(2, 3), s.Rational(1, 3), s.Rational(2, 3)]),
    s.Matrix([s.Rational(1, 3), s.Rational(-2, 3), s.Rational(2, 3)]),
    s.Matrix([s.Rational(-4, 5), 0, s.Rational(3, 5)]),
)


@pytest.mark.parametrize("e", AXES)
@pytest.mark.parametrize("complex_part", (0, 1))
def test_normalized_scalar_projection_and_rotation_covariance(e, complex_part):
    Q = s.Matrix([[1, 2, 3], [2, 4, -1], [3, -1, 6]]) + complex_part * s.I * s.Matrix(
        [[2, 1, 0], [1, -3, 4], [0, 4, 1]]
    )
    E = scalar.scalar_axis(e)
    u, v = scalar.scalar_components(Q, e)
    remainder = Q - u * scalar.TRACE - v * E
    zero(s.trace(E * E) - 1)
    zero(s.trace(E))
    zero(
        s.trace(Q.conjugate().T * Q)
        - s.conjugate(u) * u
        - s.conjugate(v) * v
        - s.trace(remainder.conjugate().T * remainder)
    )
    R = s.Matrix([[0, -1, 0], [0, 0, 1], [-1, 0, 0]])
    assert R.det() == 1
    zero(scalar.scalar_axis(R * e) - R * E * R.T)
    zero(scalar.scalar_components(R * Q * R.T, R * e)[1] - v)


def test_wrong_detector_primitive_has_nonzero_final_flux():
    G = (0, ZERO, (t + s.Rational(1, 2)) ** 2 * s.eye(3))
    xi = s.Matrix([k.primitive(1, "retarded"), 0, 0, 0])
    D = k.gauge(xi[0], ZERO)
    h = chart.first(*G, k.a)
    flux = ward.boundary_flux(EV, h, volume_tangent(h), xi, g0)
    actual = s.integrate(
        s.cancel(volume_response(D, G) - ward.detector_integrand(xi, G, EV)),
        (t, -s.Rational(1, 2), s.Rational(1, 2)),
    )
    expected = s.cancel(
        flux.subs(t, s.Rational(1, 2)) - flux.subs(t, -s.Rational(1, 2))
    )
    zero(actual - expected)
    assert expected != 0


def test_deleting_lapse_lapse_one_point_chart_is_detected():
    D = (1, ZERO, Z3)
    h = chart.first(*D, k.a)
    first = ward.pair(volume_tangent(h), h)
    contact = ward.pair(EV, chart.second(*D, *D, k.a))
    zero(first + contact - volume_response(D, D))
    assert s.cancel(contact) != 0 and s.cancel(first) != 0
    zero(volume_response(D, D))


def test_ordered_scalar_cross_entries_not_identified_by_symmetry():
    # An algebraic symmetry counterexample, not a numerical Proca kernel.
    K = s.Matrix([[1, 2], [3, 4]])
    left, right = s.Matrix([1, 0]), s.Matrix([0, 1])
    assert (left.T * K * right)[0] != (right.T * K * left)[0]


@pytest.mark.parametrize("module", (k, ward, scalar, norms))
def test_core_packets_exact_and_gate_complete(module):
    data = module.data()
    for value in data["checks"].values():
        zero(value)
    assert all(bool(v) for v in data["gates"].values())


from p8_vacuum_affine_prepared_ward_reconstruction import audit


@pytest.mark.parametrize("name,value", list(audit.residuals().items()))
def test_every_named_exact_audit_residual(name, value):
    zero(value)


@pytest.mark.parametrize("name,call,args", audit.bad_cases())
def test_unsupported_scopes_and_frontier_mutations_rejected(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_all_prior_corrected_rows_and_nine_primitive_statuses_unchanged():
    assert audit.frontier() == audit.previous.frontier()
    assert len(audit.frontier()) == 9
    assert audit.matching()[:-1] == audit.previous.matching()
    assert len({r["id"] for r in audit.matching()}) == len(audit.matching())
    assert audit.validate_scope(audit.frontier(), audit.matching())


def test_exact_counts_and_conditional_not_full_scope():
    assert len(audit.residuals()) == 38 and audit.scalar_entry_count() == 100
    assert len(audit.gates()) == 38 and all(v is True for v in audit.gates().values())
    assert len(audit.controls()) == 9 and audit.rejected_inputs() == 147
    assert len(scalar.data()["three_ordered_kernels"]) == 3
    assert "conditional" in audit.observable()["known_piece"]
    assert "remain" in audit.observable()["remaining_kernels"]


def test_both_source_and_detector_norm_rounded_constants():
    assert norms.DETECTOR == 100 and norms.SOURCE == 10**20
    assert norms.WARD == 10**35 and norms.KNOWN == 10**118
    assert k.primitive(1, "retarded").subs(t, s.Rational(1, 2)) == 1
    assert k.primitive(1, "advanced").subs(t, -s.Rational(1, 2)) == -1
