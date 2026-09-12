"""Exact hemisphere moments and full tracefree angular coefficient."""

from functools import cache
from itertools import product

import sympy as s

from . import polarizations

N = s.Matrix(s.symbols("n0:3", real=True))
AXIS = s.Matrix([0, 0, 1])


def require_orders(orders):
    if (
        len(orders) != 3
        or any(
            isinstance(x, bool) or not isinstance(x, (int, s.Integer)) or x < 0
            for x in orders
        )
        or sum(orders) > 4
    ):
        raise ValueError(
            "Three exact nonnegative exponents of total degree at most four required"
        )
    return tuple(int(x) for x in orders)


def hemisphere_moment(i, j, k):
    return _moment(*require_orders((i, j, k)))


@cache
def _moment(i, j, k):
    if i % 2 or j % 2:
        return s.Integer(0)

    def beta(a, b):
        return s.gamma(a) * s.gamma(b) / s.gamma(a + b)

    return s.simplify(
        (-1) ** k
        * beta(s.Rational(i + 1, 2), s.Rational(j + 1, 2))
        * beta(s.Rational(k + 2, 2), s.Rational(i + j + 2, 2))
    )


def integrate_hemisphere(expr):
    poly = s.Poly(s.expand(expr), *N)
    return s.simplify(
        sum(coeff * hemisphere_moment(*orders) for orders, coeff in poly.terms())
    )


def bracket(D, G, p=AXIS):
    return (
        18 * s.trace(D * G)
        - 12 * (D * p).dot(G * p)
        - (p.T * D * p)[0] * (p.T * G * p)[0]
    )


def coefficient(D, G, p=AXIS, a=1):
    return bracket(D, G, p) / (512 * s.pi**2 * a)


@cache
def data():
    x = s.symbols("D0:5", real=True)
    y = s.symbols("G0:5", real=True)
    D = s.Matrix([[x[0], x[1], x[2]], [x[1], x[3], x[4]], [x[2], x[4], -x[0] - x[3]]])
    G = s.Matrix([[y[0], y[1], y[2]], [y[1], y[3], y[4]], [y[2], y[4], -y[0] - y[3]]])
    a = s.symbols("a", positive=True)
    actual = (
        integrate_hemisphere(polarizations.leading_density(D, G, N, a))
        / (2 * s.pi) ** 3
    )
    basis = [
        s.diag(1, -1, 0) / s.sqrt(2),
        s.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / s.sqrt(2),
        s.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / s.sqrt(2),
        s.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / s.sqrt(2),
        s.diag(-1, -1, 2) / s.sqrt(6),
    ]
    eigen = s.diag(18, 18, 12, 12, s.Rational(28, 3))
    I2 = s.Matrix(3, 3, lambda i, j: integrate_hemisphere(N[i] * N[j]))
    tensor = s.MutableDenseNDimArray.zeros(3, 3, 3, 3)
    delta = lambda i, j: int(i == j)
    for i, j, k, l in product(range(3), repeat=4):
        expected = (
            s.pi
            / 24
            * (
                delta(i, j) * delta(k, l)
                + delta(i, k) * delta(j, l)
                + delta(i, l) * delta(j, k)
                + delta(i, j) * AXIS[k] * AXIS[l]
                + delta(i, k) * AXIS[j] * AXIS[l]
                + delta(i, l) * AXIS[j] * AXIS[k]
                + delta(j, k) * AXIS[i] * AXIS[l]
                + delta(j, l) * AXIS[i] * AXIS[k]
                + delta(k, l) * AXIS[i] * AXIS[j]
                - AXIS[i] * AXIS[j] * AXIS[k] * AXIS[l]
            )
        )
        tensor[i, j, k, l] = s.simplify(
            integrate_hemisphere(N[i] * N[j] * N[k] * N[l]) - expected
        )
    return {
        "weighted_measure": "Integrate w(n)=(-n.phat)_+ over the unit sphere. Its zeroth, second and fourth moments are retained. For phat=e3, the monomial formula follows from separate beta integrals in azimuth and u=-n3.",
        "invariant_coefficient": "The full current coefficient is [18tr(DG)-12(D phat).(G phat)-(phat.D.phat)(phat.G.phat)]/(512 pi^2 a). It includes the Fourier measure (2pi)^-3 and the inverse pair-frequency normalization.",
        "complete_five_channels": "In the Frobenius-orthonormal tensor/tensor/vector/vector/scalar basis relative to phat, the bracket has eigenvalues18,18,12,12,28/3. It is positive definite on real tracefree spatial tensors, not just on TT.",
        "checks": {
            "hemisphere_zeroth_moment": hemisphere_moment(0, 0, 0) - s.pi,
            "hemisphere_complete_second_moment": I2
            - s.pi * (s.eye(3) + AXIS * AXIS.T) / 4,
            "hemisphere_complete_fourth_moment": s.Matrix(
                9, 9, list(tensor.reshape(81))
            ),
            "complete_weighted_actual_symbol": s.expand(
                actual - coefficient(D, G, a=a)
            ),
            "full_five_channel_orthonormality": s.Matrix(
                5, 5, lambda i, j: s.trace(basis[i] * basis[j])
            )
            - s.eye(5),
            "full_five_channel_bracket": s.Matrix(
                5, 5, lambda i, j: s.simplify(bracket(basis[i], basis[j]))
            )
            - eigen,
        },
        "gates": {
            "positive_all_five_channels": all(x > 0 for x in eigen.diagonal()),
            "longitudinal_changes_coefficient": bool(hemisphere_moment(0, 0, 4) > 0),
            "complete_original_Fourier_measure": s.simplify(
                s.pi / (64 * (2 * s.pi) ** 3) - 1 / (512 * s.pi**2)
            )
            == 0,
            "not_only_transverse_spatial_sector": len(basis) == 5,
        },
    }
